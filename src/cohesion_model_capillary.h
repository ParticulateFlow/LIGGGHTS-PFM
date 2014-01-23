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
#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_CAPILLARY,capillary,4)
#else
#ifndef COHESION_MODEL_CAPILLARY_H_
#define COHESION_MODEL_CAPILLARY_H_
#include "contact_models.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "global_properties.h"

namespace ContactModels {
  template<typename Style>
  class CohesionModel<COHESION_CAPILLARY, Style> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), liquidVolume(0.0), surfaceTension(0.0)
    {
      history_offset = hsetup->add_value("contflag", "0");
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("liquidVolume", &MODEL_PARAMS::createLiquidVolume);
      registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
      registry.connect("liquidVolume", liquidVolume);
      registry.connect("surfaceTension", surfaceTension);
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      //r is the distance between the sphere's centeres
      const double r = cdata.r;
      const double radi = cdata.radi;
      const double radj = cdata.radj;

      if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
      double * const shear = &cdata.contact_history[history_offset];
      const double dist = r - (radi + radj);

      const double volLi = (4./3.)*M_PI*radi*radi*radi*liquidVolume;   //Modified
      const double volLj = (4./3.)*M_PI*radj*radj*radj*liquidVolume;   //Modified
      const double volBond = (volLi+volLj)/8.;                          //Modified
      //volBond = 0.00000000001;
      //const double deltaMax = pow(volBond,(1./3.));                    //Modified

      const double R2 = (radi >=radj) ? radi : radj;
      const double lBond = volBond/(R2*R2*R2);
      const double Aparam = -1.1*pow((lBond),-0.53);
      const double Bparam = -0.0082*log(lBond)+0.48;
      const double Cparam = 0.0018*log(lBond)+0.078;
      const double delta = dist/R2;
      double Fn_coh;

      if (delta > 0)
        Fn_coh = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Aparam*delta+Bparam)+Cparam);
      else if (delta <= 0)
        Fn_coh = - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Bparam)+Cparam);
      cdata.Fn += Fn_coh;

      // apply normal force
      const double fx = Fn_coh * cdata.en[0];
      const double fy = Fn_coh * cdata.en[1];
      const double fz = Fn_coh * cdata.en[2];

      i_forces.delta_F[0] += fx;
      i_forces.delta_F[1] += fy;
      i_forces.delta_F[2] += fz;

      j_forces.delta_F[0] -= fx;
      j_forces.delta_F[1] -= fy;
      j_forces.delta_F[2] -= fz;

      // store for noCollision
      shear[0] = 1.0;
    }

    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) {
      const double r = sqrt(cdata.rsq);
      const double rinv = 1/r;
      const double radi = cdata.radi;
      const double radj = cdata.radj;
      const double dist = r - (radi + radj);

      const double volLi = (4./3.)*M_PI*radi*radi*radi*liquidVolume;
      const double volLj = (4./3.)*M_PI*radj*radj*radj*liquidVolume;
      const double volBond = (volLi+volLj)/8.;
      const double deltaMax = pow(volBond,(1./3.));
      if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;
      double * const shear = &cdata.contact_history[history_offset];

      if (dist < deltaMax && MathExtraLiggghts::compDouble(shear[0],1.0,1e-6))
      {
        const double R2 = (radi >=radj) ? radi : radj;
        const double lBond = volBond/(R2*R2*R2);
        const double Aparam = -1.1*pow((lBond),-0.53);
        const double Bparam = -0.0082*log(lBond)+0.48;
        const double Cparam = 0.0018*log(lBond)+0.078;
        const double delta = dist/R2;
        double Fn_coh;

        if (delta > 0)
          Fn_coh = -M_PI*surfaceTension*sqrt(radi*radj)*(exp(Aparam*delta+Bparam)+Cparam);
        else if (delta <= 0)
          Fn_coh = -M_PI*surfaceTension*sqrt(radi*radj)*(exp(Bparam)+Cparam);

        // apply normal force
        const double fx = Fn_coh * cdata.delta[0] * rinv;
        const double fy = Fn_coh * cdata.delta[1] * rinv;
        const double fz = Fn_coh * cdata.delta[2] * rinv;

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;

        cdata.has_force_update = true;
      }
      else
      {
        shear[0]=0.0;
      }
    }

  private:
    double liquidVolume;
    double surfaceTension;
    int history_offset;
  };
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
