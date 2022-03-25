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
#include "fix_property_atom.h"
#include <math.h>

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<HERTZ_LUBRICATED> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      hminSigma(NULL),
      hco(1.0),
      limitForce(false),
      displayedSettings(false)
    {
      history_offset = hsetup->add_history_value("deltav0", "0");
      hsetup->add_history_value("hmin", "0");

      pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

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

      registry.connect("Yeff", Yeff,"model hertz/lubricated");
      registry.connect("Geff", Geff,"model hertz/lubricated");
      registry.connect("betaeff", betaeff,"model hertz/lubricated");
      registry.connect("hminSigma", hminSigma,"model hertz/lubricated");
      registry.connect("hco", hco, "model hertz/lubricated");

      if (correctYoungsModulus) {
        registry.registerProperty("YeffOriginal", &MODEL_PARAMS::createYeffOriginal,"model hertz/lubricated");
        registry.connect("YeffOriginal", YeffOriginal,"model hertz/lubricated");
        /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "HERTZ/LUBRICATED using YoungsModulusOriginal\n");
      }

      fix_visc = static_cast<FixPropertyAtom*>(modify->find_fix_property("fluidViscosity","property/atom","scalar",0,0,"normal_model hertz/lubricated"));
      visc = fix_visc->vector_atom;
    }

    // effective exponent for stress-strain relationship
    //NP used for area correction of heat transfer
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      double reff = cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));

      // gap height
      const double hij = cdata.is_wall ? (cdata.r-cdata.radi) : (cdata.r-cdata.radsum)

      // minimum approach distance
      const double hmin = compute_minimum_approach_distance(cdata, cdata.vn, hij, reff);

      // contact force
      double kn, kt, gamman, gammat;
      double Fc = compute_contact_force(cdata, cdata.vn, reff, cdata.meff, hij, hmin, kn, kt, gamman, gammat);
      
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;
      cdata.Fn = Fc;

      // during collision(), total force is the contact force
      double Fn = Fc;

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
      if(cdata.touch) *cdata.touch |= TOUCH_NORMAL_MODEL;
    }

    void noCollision(ContactData& cdata, ForceData & i_forces, ForceData & j_forces)
    {
      // compute approach properties
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double reff = cdata.is_wall ? ri : (ri*rj/(ri+rj));

      // contact radius
      const double rc = cdata.is_wall ? ri : (ri+rj);

      double cutoff = hco*reff + rc; // cutoff distance

      double Fn = 0.;

      if (cdata.rsq <= cutoff*cutoff)
      {
        const double r = sqrt(cdata.rsq);
        const double hij = r-rc;

        const int i = cdata.i;
        const int j = cdata.j;
        const int itype = cdata.itype;
        const int jtype = cdata.jtype;

        // relative velocity
        const double vr1 = cdata.v_i[0] - cdata.v_j[0];
        const double vr2 = cdata.v_i[1] - cdata.v_j[1];
        const double vr3 = cdata.v_i[2] - cdata.v_j[2];

        // normal vector
        const double rinv = 1./r;
        const double enx = cdata.delta[0]*rinv;
        const double eny = cdata.delta[1]*rinv;
        const double enz = cdata.delta[2]*rinv;

        // normal velocity
        const double vn = vr1 * enx + vr2 * eny + vr3 * enz;

        // minimum approach distance
        const double hmin = compute_minimum_approach_distance(cdata, vn, hij, reff);

        // viscosity
        visc = fix_visc->vector_atom;
        const double etaf = cdata.is_wall ? visc[i] : ((visc[i] + visc[j])/2.); 

        // lubrication force
        double Fl = -6.*M_PI*etaf*vn*reff*reff/MAX(hij,hmin);

        // contact force
        double Fc = 0;
        if (hij<=hmin)
        {
          // compute meff
          double *rmass = atom->rmass;
          double *mass = atom->mass;
          int *mask = atom->mask;

          double meff;
          if (cdata.is_wall)
          {
            if (rmass)
              meff = rmass[i];
            else
              meff = mass[itype];
          }
          else
          {
            double mi, mj;
            if (rmass) {
              mi = rmass[i];
              mj = rmass[j];
            } else {
              mi = mass[itype];
              mj = mass[jtype];
            }
            if (pair_gran->fr_pair()) {
              const double * mass_rigid = pair_gran->mr_pair();
              if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
              if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
            }

            meff = mi * mj / (mi + mj);
            if (mask[i] & pair_gran->freeze_group_bit())
              meff = mj;
            if (mask[j] & pair_gran->freeze_group_bit())
              meff = mi;
          }

          // contact force calculation
          double kn, kt, gamman, gammat;
          Fc = compute_contact_force(cdata, vn, reff, meff, hij, hmin, kn, kt, gamman, gammat);
        }

        // total normal force
        if (hij>hmin)
          Fn = Fl;
        else 
        {
          double frac = hij / hmin;
          Fn = frac*Fl + (1-frac)*Fc;
        }

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
        cdata.has_force_update = true;
        if(cdata.touch) *cdata.touch |= TOUCH_NORMAL_MODEL;
      }
      else
      {
        if(cdata.touch) *cdata.touch &= ~TOUCH_NORMAL_MODEL;
        if(!cdata.contact_history) return;
        double * const deltav0 = &cdata.contact_history[history_offset];
        double * const hmin  = &cdata.contact_history[history_offset+1];
        *deltav0 = 0.0;
        *hmin = 0.0;
      }
    }

    double compute_minimum_approach_distance(ContactData & cdata, double vn, double hij, double reff)
    {
      double * const deltav0 = &cdata.contact_history[history_offset];
      double * const hmin  = &cdata.contact_history[history_offset+1];

      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      // revaluate at first contact, or when approach velocity increased
      bool revaluate = true;
      if (cdata.touch && (*cdata.touch & TOUCH_NORMAL_MODEL))
        if ((hij<=*hmin) || (-vn<=*deltav0))
          revaluate = false;

      if (revaluate)
      {
        // approach velocity
        *deltav0 = MAX(*deltav0,MAX(-vn,0.));

        // effective youngs modulus
        const double YoungsModulusEff = correctYoungsModulus ? YeffOriginal[itype][jtype] : Yeff[itype][jtype];

        // viscosity
        visc = fix_visc->vector_atom;
        const double etaf = cdata.is_wall ? visc[cdata.i] : ((visc[cdata.i] + visc[cdata.j])/2.); 

        // elastic approach distance
        const double hmine = 0.37 * pow(etaf**deltav0/YoungsModulusEff,0.4) * pow(reff,0.6);
        *hmin = MAX(hminSigma[itype][jtype],hmine);
      }

      return *hmin;
    }

    double compute_contact_force(ContactData & cdata, double vn, double reff, double meff, double hij, double hmin, double &kn, double &kt, double &gamman, double &gammat)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      // overlap
      const double deltan = MAX(hmin - hij, 0.);
      
      const double sqrtval = sqrt(reff*deltan);

      const double Sn=2.*Yeff[itype][jtype]*sqrtval;
      const double St=8.*Geff[itype][jtype]*sqrtval;

      kn=4./3.*Yeff[itype][jtype]*sqrtval;
      kt=St;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
      gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
      /*NL*/ //if (screen) fprintf(screen,"tangential_damping %s\n",tangential_damping?"yes":"no");
      if (!tangential_damping) gammat = 0.0;

      if(!displayedSettings)
      {
        displayedSettings = true;
      }
      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fc_damping = -gamman*vn;
      const double Fc_contact = kn*deltan;
      double Fc = Fc_damping + Fc_contact;

      //limit force to avoid the artefact of negative repulsion force
      if (limitForce && (Fc<0.0))
        Fc = 0.;
        
      return Fc;
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

    double * visc;
    FixPropertyAtom* fix_visc;

    int history_offset;
    PairGran *pair_gran;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
    bool correctYoungsModulus;
  };

}

}
#endif
#endif
