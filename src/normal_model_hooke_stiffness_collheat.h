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
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HOOKE_STIFFNESS_COLLHEAT,hooke/stiffness/collheat,6)
#else
#ifndef NORMAL_MODEL_HOOKE_STIFFNESS_COLLHEAT_H_
#define NORMAL_MODEL_HOOKE_STIFFNESS_COLLHEAT_H_
#include "fix_property_atom.h"
#include "math.h"

namespace MODEL_PARAMS
{
  
  VectorProperty * createThermalConductivity(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    LAMMPS * lmp = registry.getLAMMPS();
    const int max_type = registry.max_type();

    VectorProperty * vec = new VectorProperty(max_type+1);
    FixPropertyGlobal * kappa = registry.getGlobalProperty("thermalConductivity","property/global","peratomtype",max_type,0,caller);

    for(int i=1; i < max_type+1; i++)
    {
      const double kappa_i = kappa->compute_vector(i-1);

      // error checks on kappa
      if(sanity_checks)
      {
        if(kappa_i < 0.0)
          lmp->error->all(FLERR,"thermal conductivity >= 0 required");
      }

      vec->data[i] = kappa_i;
    }

    return vec;
  }
  
  VectorProperty * createThermalCapacity(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    LAMMPS * lmp = registry.getLAMMPS();
    const int max_type = registry.max_type();

    VectorProperty * vec = new VectorProperty(max_type+1);
    FixPropertyGlobal * c = registry.getGlobalProperty("thermalCapacity","property/global","peratomtype",max_type,0,caller);

    for(int i=1; i < max_type+1; i++)
    {
      const double c_i = c->compute_vector(i-1);

      // error checks on c
      if(sanity_checks)
      {
        if(c_i < 0.0)
          lmp->error->all(FLERR,"thermal capacity >= 0 required");
      }

      vec->data[i] = c_i;
    }

    return vec;
  }
  
}

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_STIFFNESS_COLLHEAT> : protected NormalModel<HOOKE_STIFFNESS>
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup* hsetup) : NormalModel<HOOKE_STIFFNESS>(lmp, hsetup),
      history_offset(0)
    {
      history_offset = hsetup->add_history_value("contflag", "0");
      /*NL*/ if(comm->me == 0) fprintf(screen, "HOOKE/STIFFNESS/COLLHEAT loaded\n");
    }

    void registerSettings(Settings & settings)
    {
      NormalModel<HOOKE_STIFFNESS>::registerSettings(settings);
    }

    void connectToProperties(PropertyRegistry & registry) {
      NormalModel<HOOKE_STIFFNESS>::connectToProperties(registry);
      
      registry.registerProperty("thermalConductivity", &MODEL_PARAMS::createThermalConductivity);
      registry.registerProperty("thermalCapacity", &MODEL_PARAMS::createThermalCapacity);
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      
      registry.connect("thermalConductivity", thermalConductivity, "normal_model hooke stiffness collheat");
      registry.connect("thermalCapacity", thermalCapacity, "normal_model hooke stiffness collheat");
      registry.connect("Yeff", Yeff,"normal_model hooke stiffness collheat");

      fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,"normal_model hooke stiffness collheat"));
      Temp = fix_temp->vector_atom;

    }
    
    // effective exponent for stress-strain relationship
    //NP used for area correction of heat transfer
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      NormalModel<HOOKE_STIFFNESS>::collision(cdata, i_forces, j_forces);      
      
      // first-contact detection could be implemented more elegantly
      // using *cdata.touch |= TOUCH_NORMAL_MODEL ...
      
      double * const contflag = &cdata.contact_history[history_offset];
      
      // check how contflag[0] is initialized!!!
      if (contflag[0]<0.5)
      {
	const int i = cdata.i;
        const int j = cdata.j;
        const int itype = cdata.itype;
        const int jtype = cdata.jtype;

        const double radi = cdata.radi;
        const double radj = cdata.radj;
        const double reff = radi*radj/(radi+radj);
	const double mi = cdata.mi;
	const double mj = cdata.mj;
        const double meff = cdata.meff;
	const double vrel = fabs(cdata.vn);
	
	const double densityi = atom->density[i];
	const double densityj = atom->density[j];	
	const double kappai = thermalConductivity[itype];
	const double kappaj = thermalConductivity[jtype];
	const double ci = thermalCapacity[itype];
	const double cj = thermalCapacity[jtype];
	
	const double Ac = M_PI*0.9745*pow(meff/Yeff[itype][jtype],0.4)*pow(reff*vrel,0.8);
	const double tc = 2.87*pow(meff/Yeff[itype][jtype],0.4)*pow(reff*vrel,-0.2);
	
	const double numerator = 0.87 * (Temp[i] - Temp[j]) * Ac * sqrt(tc);
	const double denominator = 1.0/sqrt(densityi*ci*kappai) + 1.0/sqrt(densityj*cj*kappaj);
	const double qc = numerator/denominator;
	
	const double dTempi = - qc / (ci * mi);
	const double dTempj = qc / (cj * mj);
	
	Temp[i] += dTempi;
	Temp[j] += dTempj;	
      }
      // store for noCollision
      contflag[0] = 1.0;    
    }

    void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      double * const contflag = &cdata.contact_history[history_offset];
      // store for collision
      contflag[0] = 0.0;
    }
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    
    int history_offset;
    
    double ** Yeff;
    double * thermalConductivity;
    double * thermalCapacity;
    double * Temp;
    
    FixPropertyAtom* fix_temp;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_STIFFNESS_COLLHEAT_H_
#endif
