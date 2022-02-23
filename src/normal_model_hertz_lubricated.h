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

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      limitForce(false),
      displayedSettings(false)
    {
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "HERTZ/LUBRICATED loaded\n");
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff,"model hertz/lubricated");
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff,"model hertz/lubricated");
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff,"model hertz/lubricated");

      registry.connect("Yeff", Yeff,"model hertz/lubricated");
      registry.connect("Geff", Geff,"model hertz/lubricated");
      registry.connect("betaeff", betaeff,"model hertz/lubricated");
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
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(cdata.is_non_spherical && atom->superquadric_flag)
        reff = cdata.reff;
#endif
      const double meff = cdata.meff;
      const double deltan = cdata.deltan;

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

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*deltan;
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if(limitForce && (Fn<0.0) )
      {
          Fn = 0.0;
      }

      cdata.Fn = Fn;
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;

#ifdef NONSPHERICAL_ACTIVE_FLAG
      double torque_i[3] = {0., 0., 0.};
      double Fn_i[3] = { Fn * cdata.en[0], Fn * cdata.en[1], Fn * cdata.en[2] };
      if(cdata.is_non_spherical) {
        double xci[3];
        vectorSubtract3D(cdata.contact_point, atom->x[cdata.i], xci);
        vectorCross3D(xci, Fn_i, torque_i);
      }
#endif

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
#ifdef NONSPHERICAL_ACTIVE_FLAG
        if(cdata.is_non_spherical) {
          // for non-spherical particles normal force can produce torque!
          i_forces.delta_torque[0] += torque_i[0];
          i_forces.delta_torque[1] += torque_i[1];
          i_forces.delta_torque[2] += torque_i[2];
        }
#endif
      } else {
        i_forces.delta_F[0] = cdata.Fn * cdata.en[0];
        i_forces.delta_F[1] = cdata.Fn * cdata.en[1];
        i_forces.delta_F[2] = cdata.Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
#ifdef NONSPHERICAL_ACTIVE_FLAG
        if(cdata.is_non_spherical) {
          // for non-spherical particles normal force can produce torque!
          double xcj[3], torque_j[3];
          double Fn_j[3] = { -Fn_i[0], -Fn_i[1], -Fn_i[2]};
          vectorSubtract3D(cdata.contact_point, atom->x[cdata.j], xcj);
          vectorCross3D(xcj, Fn_j, torque_j);

          i_forces.delta_torque[0] += torque_i[0];
          i_forces.delta_torque[1] += torque_i[1];
          i_forces.delta_torque[2] += torque_i[2];

          j_forces.delta_torque[0] += torque_j[0];
          j_forces.delta_torque[1] += torque_j[1];
          j_forces.delta_torque[2] += torque_j[2];
        }
#endif
      }
    }

    void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** betaeff;

    bool tangential_damping;
    bool limitForce;
    bool displayedSettings;
  };

}

}
#endif
#endif
