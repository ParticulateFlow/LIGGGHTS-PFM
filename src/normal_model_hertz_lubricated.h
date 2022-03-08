/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HERTZ_LUBRICATED,hertz/lubricated,9)
#else
#ifndef NORMAL_MODEL_HERTZ_LUBRICATED_H_
#define NORMAL_MODEL_HERTZ_LUBRICATED_H_
#include "global_properties.h"
#include <math.h>

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<HERTZ_LUBRICATED> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      hminSigma(NULL),
      hco(1.0),
      coeffMu(NULL),
      limitForce(false),
      displayedSettings(false)
    {
      history_offset = hsetup->add_history_value("deltav0", "0");
      hsetup->add_history_value("hmin", "1");

#ifdef NONSPHERICAL_ACTIVE_FLAG
      error->all(FLERR,"HERTZ/LUBRICATED not defined for non-spherical particles\n");
#endif
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      error->all(FLERR,"HERTZ/LUBRICATED not defined for supersquadric particles\n");
#endif      
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "HERTZ/LUBRICATED loaded\n");
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("correctYoungsModulus", correctYoungsModulus);
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model hertz/lubricated");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model hertz/lubricated");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz/lubricated");
      registry.registerProperty("hminSigma", &MODEL_PARAMS::createHminSigma,"model hertz/lubricated");
      registry.registerProperty("hco", &MODEL_PARAMS::createLubricationCutoff,"model hertz/lubricated");
      registry.registerProperty("coeffMu", &MODEL_PARAMS::createCoeffMu,"model hertz/lubricated");

      registry.connect("Yeff", Yeff,"model hertz/lubricated");
      registry.connect("Geff", Geff,"model hertz/lubricated");
      registry.connect("betaeff", betaeff,"model hertz/lubricated");
      registry.connect("hminSigma", hminSigma,"model hertz/lubricated");
      registry.connect("hco", hco, "model hertz/lubricated");
      registry.connect("coeffMu", coeffMu,"model hertz/lubricated");

      if (correctYoungsModulus) {
        registry.registerProperty("YeffOriginal", &MODEL_PARAMS::createYeffOriginal,"model hertz/lubricated");
        registry.connect("YeffOriginal", YeffOriginal,"model hertz/lubricated");
        /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "HERTZ/LUBRICATED using YoungsModulusOriginal\n");
      }
    }

    // effective exponent for stress-strain relationship
    //NP used for area correction of heat transfer
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      double reff = cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      double * const deltav0 = &cdata.contact_history[history_offset];
      double * const hmin  = &cdata.contact_history[history_offset+1];

      // gap height
      double hij;
      if(cdata.is_wall)
        hij = cdata.r - cdata.radi;
      else
        hij = cdata.r - cdata.radsum;
      
      double Fn = 0.; // total normal force
      double Fl = 0.; // lubrication force
      double Fc = 0.; // contact force

      if (hij <= hco*reff) 
      {
        // approach velocity
        const double oldDeltav0 = *deltav0;
        if (cdata.vn<0. || hij<*hmin)
          *deltav0 = MAX(*deltav0,-cdata.vn);
        else
          *deltav0 = 0.;
        
        const double mul = coeffMu[itype][jtype]; //FIXME: couple with CFD

        // approach distance
        if (*deltav0!=oldDeltav0 && *deltav0>0.) {

          double YoungsModulusEff;
          if (correctYoungsModulus)
            YoungsModulusEff = YeffOriginal[itype][jtype];
          else
            YoungsModulusEff = Yeff[itype][jtype];

          const double hmine = 0.37 * pow(mul**deltav0/YoungsModulusEff,0.4) * pow(reff,0.6);
          *hmin = MAX(hminSigma[itype][jtype],hmine);
        }

        // lubrication force
        if (hij>0.)
          Fl = -6.*M_PI*mul*cdata.vn*reff*reff/MAX(hij,*hmin);

        // contact force
        if (hij<=*hmin)
        {
          // overlap
          const double deltan = MAX(*hmin - hij, 0.);
          
          const double meff = cdata.meff;
          const double sqrtval = sqrt(reff*deltan);

          const double Sn=2.*Yeff[itype][jtype]*sqrtval;
          const double St=8.*Geff[itype][jtype]*sqrtval;

          double kn=4./3.*Yeff[itype][jtype]*sqrtval;
          double kt=St;
          const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
          const double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
          double gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
          /*NL*/ //if (screen) fprintf(screen,"tangential_damping %s\n",tangential_damping?"yes":"no");
          if (!tangential_damping) gammat = 0.0;

          if(!displayedSettings)
          {
            displayedSettings = true;
          }
          // convert Kn and Kt from pressure units to force/distance^2
          kn /= force->nktv2p;
          kt /= force->nktv2p;

          const double Fc_damping = -gamman*cdata.vn;
          const double Fc_contact = kn*deltan;
          Fc = Fc_damping + Fc_contact;

          //limit force to avoid the artefact of negative repulsion force
          if(limitForce && (Fc<0.0) )
          {
            Fc = 0.;
          }

          cdata.kn = kn;
          cdata.kt = kt;
          cdata.gamman = gamman;
          cdata.gammat = gammat;
        }

        // add to total force
        if (hij>*hmin)
          Fn = Fl;
        else if (hij<=0.)
          Fn = Fc;
        else 
        {
          double frac = hij / *hmin;
          Fn = frac*Fl + (1-frac)*Fc;
        }
        cdata.Fn = Fc; // store only contact force for friction calculation in tangential model
      }  

      fprintf(screen,"h: %g\t",hij);
      fprintf(screen,"vn: %g\t",cdata.vn);
      fprintf(screen,"vn0: %g\t",*deltav0);
      fprintf(screen,"hmin: %g\t",*hmin);
      fprintf(screen,"Fl: %g\t",Fl); 
      fprintf(screen,"Fc: %g\t",Fc);
      fprintf(screen,"Fn: %g\t",Fn);

      // apply total normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
      } else {
        i_forces.delta_F[0] = Fn * cdata.en[0];
        i_forces.delta_F[1] = Fn * cdata.en[1];
        i_forces.delta_F[2] = Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** YeffOriginal;
    double ** Geff;
    double ** betaeff;
    double ** hminSigma;
    double hco;
    double ** coeffMu;

    int history_offset;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    bool correctYoungsModulus;
  };

}

}
#endif
#endif
