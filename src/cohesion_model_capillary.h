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
   Stefan Radl (TU Graz)
   Bhageshvar Mohan (TU Graz)
------------------------------------------------------------------------- */
#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_CAPILLARY,capillary,4)
#else
#ifndef COHESION_MODEL_CAPILLARY_H_
#define COHESION_MODEL_CAPILLARY_H_
#include "contact_models.h"
#include <math.h>
#include "math_extra_liggghts.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_CAPILLARY> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), 
      liquidVolume(0.0), 
      surfaceTension(0.0), 
      capillaryModel(0)
    {
      // "contflag" is required to indicate bridge
      history_offset = hsetup->add_history_value("contflag", "0");
    }

    void registerSettings(Settings& settings) {
      NamedIntegerSetting* capillaryModelSetting = settings.registerNamedIntegerSetting("capillaryModel", capillaryModel, 0);
      capillaryModelSetting->addOption("Mikami", 0);
      capillaryModelSetting->addOption("Willett", 1);
      capillaryModelSetting->setDefault("Mikami");
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("liquidVolume", &MODEL_PARAMS::createLiquidVolume);
      registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
      registry.connect("liquidVolume", liquidVolume,"cohesion_model capilary");
      registry.connect("surfaceTension", surfaceTension,"cohesion_model capilary");
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
      double * const hist = &cdata.contact_history[history_offset];
      const double dist = r - (radi + radj);
      const double R2 = (radi >=radj) ? radi : radj;
      const double delta = dist/R2;

      const double volLi = (4./3.)*M_PI*radi*radi*radi*liquidVolume;  
      const double volLj = (4./3.)*M_PI*radj*radj*radj*liquidVolume; 

      const double volBond = (volLi+volLj);
      const double lBond = volBond/(R2*R2*R2);

      double Fn_coh = 0.0;
    
      //Section - Model Selection
      switch(capillaryModel) {
        case 0:
          Fn_coh = compute_force_mikami(lBond, radi, radj, delta);
          break;

        case 1:
          Fn_coh = compute_force_willet(lBond, radi, radj, delta);
          break;
      }

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
      hist[0] = 1.0;
    }

    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) {
      const double r = sqrt(cdata.rsq);
      const double rinv = 1/r;
      const double radi = cdata.radi;
      const double radj = cdata.radj;
      const double dist = r - (radi + radj);
      double * const hist = &cdata.contact_history[history_offset];
      const double R2 = (radi >=radj) ? radi : radj;
      double Fn_coh(0.0);
      
      const double volLi = (4./3.)*M_PI*radi*radi*radi*liquidVolume;
      const double volLj = (4./3.)*M_PI*radj*radj*radj*liquidVolume;

      const double volBond = (volLi+volLj);
      const double lBond = volBond/(R2*R2*R2);
      const double delta = dist/R2;

      const double deltaMax = pow(volBond,(1./3.));

      //if liquid bridge exists "no_collision" is true. Force calculation only if liquid bridge exists. 
      if((dist < deltaMax) && (MathExtraLiggghts::compDouble(hist[0],1.0,1e-6)))
      {
        switch(capillaryModel) {
          case 0:
            Fn_coh = compute_force_mikami(lBond, radi, radj, delta);
            break;

          case 1:
            Fn_coh = compute_force_willet(lBond, radi, radj, delta);
            break;
        }

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
        if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
      }
      else
      {
        hist[0] = 0.0;

        if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;
      }
    }

    /**
     * Mikami's Capillary Force Model
     * Mikami et al., CES 1998, 50(16)
     */
    double compute_force_mikami(double lBond, double radi, double radj, double delta) {
      const double Aparam = -1.1*pow((std::max(1e-64,lBond)),-0.53);
      const double Bparam = -0.0082*log(std::max(1e-64,lBond))+0.48;
      const double Cparam = 0.0018*log(std::max(1e-64,lBond))+0.078;

      if (delta > 0) {
        return - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Aparam*delta+Bparam)+Cparam);
      }
      return - M_PI*surfaceTension*sqrt(radi*radj)*(exp(Bparam)+Cparam);
    }

    /*
     * Willett's Capillary Force Model
     * Willett et al., Langmuir 2000,16
     */
    double compute_force_willet(double lBond, double radi, double radj, double delta) {
      double Fn_coh = 0.0;
      double tmpVar  = 2.3*log(std::max(1e-64,lBond));
      double tmpVar2 = tmpVar*tmpVar;
      double tmpVar3 = tmpVar*tmpVar2;

      double f1 = (-0.44507)
         + (-0.1119)    * tmpVar
         + (-0.0121010) * tmpVar2
         + (-0.000500)  * tmpVar3;

      double f2 = (1.9222)
         + (-0.0668)    * tmpVar
         + (-0.0013375) * tmpVar2;

      double f3 = (1.268)
         + (0.19800)    * tmpVar
         + (0.02232000) * tmpVar2
         + (0.0008585)  * tmpVar3;

      double f4 = (-0.010703)
         + (0.03345)    * tmpVar
         + (0.00185740) * tmpVar2;

      if (delta > 0)
      {
        double tmpVar_1= 2.3 *log(0.5*delta);
        Fn_coh = - 2.0 * M_PI
                       * surfaceTension
                       * sqrt(radi *radj)
                       * (
                          exp(
                                f1
                               -(
                                    f2 * exp(
                                               f3 * tmpVar_1
                                              +f4 * tmpVar_1 * tmpVar_1
                                            )
                                )
                             )
                         );
      }
      else
      {
        Fn_coh = - 2.0 * M_PI
                       * surfaceTension
                       * sqrt(radi *radj)
                       * exp(f1);
      }

      return Fn_coh;
    }


  private:
    double liquidVolume, surfaceTension;
    int history_offset;
    int capillaryModel;
  };
}
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
