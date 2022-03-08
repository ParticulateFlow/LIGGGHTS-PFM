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
   Tim M.J. Nijssen (TU Delft)
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

      // gap height
      double hij;
      if(cdata.is_wall)
        hij = cdata.r - cdata.radi;
      else
        hij = cdata.r - cdata.radsum;
      fprintf(screen,"h: %g\t",hij);  
      
      double Fn = 0.; // total normal force

      if (hij <= hco*reff) 
      {
        fprintf(screen,"vn: %g\t",cdata.vn);

        // minimum approach distance
        const double hmin = compute_minimum_approach_distance(cdata, itype, jtype, cdata.vn, hij, reff);

        // lubrication force
        double Fl = compute_lubrication_force(itype, jtype, cdata.vn, reff, hij, hmin);

        // contact force
        double Fc = 0.;
        if (hij<=hmin)
        {
          // overlap
          const double deltan = MAX(hmin - hij, 0.);
          
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
          if (limitForce && (Fc<0.0))
            Fc = 0.;

          fprintf(screen,"Fc: %g\t",Fc);

          cdata.kn = kn;
          cdata.kt = kt;
          cdata.gamman = gamman;
          cdata.gammat = gammat;
        }
        cdata.Fn = Fc; // store only contact force for friction calculation in tangential model

        // total normal force
        Fn = compute_total_normal_force(Fl, Fc, hij, hmin);

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
      else 
      {
        // apply no force
        for (int ii=0; ii<3; ii++)
        {
          i_forces.delta_F[ii] = 0.;
          if(!cdata.is_wall)
            j_forces.delta_F[ii] = 0.;
        }
      }
      fprintf(screen,"\n");  
    }

    void noCollision(ContactData& cdata, ForceData & i_forces, ForceData & j_forces)
    {
      // compute approach properties
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double reff = cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      const double r = sqrt(cdata.rsq);

      // gap height
      double hij;
      if(cdata.is_wall)
        hij = r - cdata.radi;
      else
        hij = r - cdata.radsum;
      fprintf(screen,"h: %g\t",hij);

      double Fn = 0.;

      if (hij <= hco*reff)
      {
        int *type = atom->type;
        const int i = cdata.i;
        const int j = cdata.j;
        const int itype = type[i];
        const int jtype = type[j];

        // relative velocity
        double **v = atom->v;
        const double vr1 = v[i][0] - v[j][0];
        const double vr2 = v[i][1] - v[j][1];
        const double vr3 = v[i][2] - v[j][2];

        // normal vector
        double **x = atom->x;
        const double rinv = 1./r;
        const double enx = (x[i][0] - x[j][0])*rinv;
        const double eny = (x[i][1] - x[j][1])*rinv;
        const double enz = (x[i][2] - x[j][2])*rinv;

        // normal velocity
        const double vn = vr1 * enx + vr2 * eny + vr3 * enz;

        // minimum approach distance
        double hmin = compute_minimum_approach_distance(cdata, itype, jtype, vn, hij, reff);

        // lubrication force
        double Fl = compute_lubrication_force(itype, jtype, vn, reff, hij, hmin);

        // contact force
        // TODO contact force calculation
        double Fc = 0.;
        fprintf(screen,"Fc: %g\t",Fc);

        // total normal force
        Fn = compute_total_normal_force(Fl, Fc, hij, hmin);

        // apply total normal force
        if(cdata.is_wall) {
          const double Fn_ = Fn * cdata.area_ratio;
          i_forces.delta_F[0] = Fn_ * enx;
          i_forces.delta_F[1] = Fn_ * eny;
          i_forces.delta_F[2] = Fn_ * enz;
        } else {
          i_forces.delta_F[0] = Fn * enx;
          i_forces.delta_F[1] = Fn * eny;
          i_forces.delta_F[2] = Fn * enz;

          j_forces.delta_F[0] = -i_forces.delta_F[0];
          j_forces.delta_F[1] = -i_forces.delta_F[1];
          j_forces.delta_F[2] = -i_forces.delta_F[2];
        }
      }
      else
      {
        // apply no force
        for (int ii=0; ii<3; ii++)
        {
          i_forces.delta_F[ii] = 0.;
          if(!cdata.is_wall)
            j_forces.delta_F[ii] = 0.;
        }
      }
      cdata.has_force_update = true;
      fprintf(screen,"\n");
    }

    double compute_minimum_approach_distance(ContactData & cdata, int itype, int jtype, double vn, double hij, double reff)
    {
      double * const deltav0 = &cdata.contact_history[history_offset];
      double * const hmin  = &cdata.contact_history[history_offset+1];

      // approach velocity
      const double oldDeltav0 = *deltav0;
      if (vn<0. || hij<*hmin)
        *deltav0 = MAX(*deltav0,-vn);
      else
        *deltav0 = 0.;
        
      const double mul = coeffMu[itype][jtype]; //FIXME: couple with CFD

      // approach distance
      if (*deltav0!=oldDeltav0 && *deltav0>0.)
      {
        double YoungsModulusEff;
        if (correctYoungsModulus)
          YoungsModulusEff = YeffOriginal[itype][jtype];
        else
          YoungsModulusEff = Yeff[itype][jtype];

        const double hmine = 0.37 * pow(mul**deltav0/YoungsModulusEff,0.4) * pow(reff,0.6);
        *hmin = MAX(hminSigma[itype][jtype],hmine);
      }

      fprintf(screen,"vn: %g\t",vn);
      fprintf(screen,"vn0: %g\t",*deltav0);
      fprintf(screen,"hmin: %g\t",*hmin);

      return *hmin;
    }

    double compute_lubrication_force(int itype, int jtype, double vn, double reff, double hij, double hmin)
    {
      double Fl;
      if (hij>0.)
        Fl = -6.*M_PI*coeffMu[itype][jtype]*vn*reff*reff/MAX(hij,hmin);
      else
        Fl = 0.;

      fprintf(screen,"Fl: %g\t",Fl); 
      
      return Fl;
    }

    double compute_total_normal_force(double Fl, double Fc, double hij, double hmin)
    {
      double Fn;
      if (hij>hmin)
        Fn = Fl;
      else if (hij<=0.)
        Fn = Fc;
      else 
      {
        double frac = hij / hmin;
        Fn = frac*Fl + (1-frac)*Fc;
      }

      fprintf(screen,"Fn: %g\t",Fn);

      return Fn;
    }

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
