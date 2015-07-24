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
#include "math.h"
#include "math_extra_liggghts.h"
#include "global_properties.h"

namespace LIGGGHTS {
namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_CAPILLARY> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp), 
      liquidVolume(0.0), 
      surfaceTension(0.0), 
      switchModel(0.0), 
      historyIndex(0.0)
    {
      //Check if liquidtracking fix is used
      #if 1
      int i = modify->find_fix("liquidtransfer");
      if (i < 0) liquidtracking = false;
      else liquidtracking = true;
      #endif

      //If no liquid tracking fix used "contflag" is required to indicate bridge 
      if(!liquidtracking) {      
      history_offset = hsetup->add_history_value("contflag", "0");
      }
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("liquidVolume", &MODEL_PARAMS::createLiquidVolume);
      registry.registerProperty("surfaceTension", &MODEL_PARAMS::createSurfaceTension);
      registry.registerProperty("switchModel", &MODEL_PARAMS::createSwitchModel);
      registry.registerProperty("historyIndex", &MODEL_PARAMS::createHistoryIndex);
      registry.connect("liquidVolume", liquidVolume,"cohesion_model capilary");
      registry.connect("surfaceTension", surfaceTension,"cohesion_model capilary");
      registry.connect("switchModel", switchModel,"cohesion_model capilary");
      registry.connect("historyIndex", historyIndex,"cohesion_model capilary");
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {

      if (liquidtracking) history_offset = 0; //if liquidtracking fix is used no "contflag" required
      Index = static_cast<int> (historyIndex);

      //r is the distance between the sphere's centeres
      const double r = cdata.r;
      const double radi = cdata.radi;
      const double radj = cdata.radj;
      double Fn_coh(0.0);

      if(cdata.touch) *cdata.touch |= TOUCH_COHESION_MODEL;
      double * const hist = &cdata.contact_history[history_offset];
      const double dist = r - (radi + radj);
      const double R2 = (radi >=radj) ? radi : radj;
      const double delta = dist/R2;

      const double volLi = (4./3.)*M_PI*radi*radi*radi*liquidVolume;  
      const double volLj = (4./3.)*M_PI*radj*radj*radj*liquidVolume; 

      double volBond;
      if (liquidtracking) volBond = hist[Index];
      else volBond = (volLi+volLj);   
      const double lBond = volBond/(R2*R2*R2);
    
      //Section - Model Selection
     if (switchModel ==0) 
     {
      #   include "cohesion_model_capillary_model_Mikami.h"
     }

     if (switchModel ==1) 
     {
      #   include "cohesion_model_capillary_model_Willett.h"
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
    
      // store for noCollision (when no fix liquidtracking is used)
      if (!liquidtracking) hist[0] = 1.0;
    }

    void noCollision(ContactData & cdata, ForceData & i_forces, ForceData & j_forces) {

      if (liquidtracking) history_offset = 0; //if liquidtracking fix is used no "contflag" required
      Index = static_cast<int> (historyIndex);

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

      double volBond;
      if (liquidtracking) volBond = hist[Index]; //Index '3' to be used only if tangential history model is used.
      else volBond = (volLi+volLj);
      const double lBond = volBond/(R2*R2*R2);
      const double delta = dist/R2;

      no_collision = false; 

      const double deltaMax = pow(volBond,(1./3.));

      //if liquid bridge exists "no_collision" is true. Force calculation only if liquid bridge exists. 
      if (((dist < deltaMax) && (MathExtraLiggghts::compDouble(hist[0],1.0,1e-6)) && (!liquidtracking)) || ((hist[Index]>0) && (liquidtracking)))
        no_collision = true;

      if (no_collision)
      {

       if (switchModel ==0) 
       {
        #   include "cohesion_model_capillary_model_Mikami.h"
       }

       if (switchModel ==1) 
       {
        #   include "cohesion_model_capillary_model_Willett.h"
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
        if (!liquidtracking) hist[0]=0.0; //(when no fix liquidtracking is used)

        if(cdata.touch) *cdata.touch &= ~TOUCH_COHESION_MODEL;
      }
    }

  private:
    double liquidVolume, surfaceTension, switchModel;
    double historyIndex;
    int history_offset, Index;
    bool liquidtracking, no_collision;
  };
}
}
#endif // COHESION_MODEL_CAPILLARY_H_
#endif
