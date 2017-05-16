/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations
   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com
   Copyright 2015-     JKU Linz
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
   M.Efe Kinaci (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(chem/shrink/core,FixChemShrinkCore)

#else

#ifndef LMP_FIX_CHEM_SHRINKCORE_H
#define LMP_FIX_CHEM_SHRINKCORE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChemShrinkCore : public Fix  {

public:
  FixChemShrinkCore(class LAMMPS *, int, char **);
  ~FixChemShrinkCore();

  void pre_delete(bool unfixflag);
  int setmask();
  void post_create();
  virtual void updatePtrs();
  virtual void init();
  virtual void post_force(int);

  double partSurfArea(double);
  int active_layers(int);
  void calcMassLayer(int);
  double K_eq(int, double);
  void reaction(int, double *, double *, double *, double *, double *);   //, double *, double *
  void getA(int, double *, double *);
  void diffcoeff(int, double *);
  void getB(int, double *);
  void getXi(int, double *, double *);
  void layerRad(int);
  void update_atom_properties(int, double *);
  void update_gas_properties(int, double *);

  virtual void updatePtrsCreate();

 protected:
   int iarg_;

   // timestep
   double TimeStep;
   // modified strings of species concentrations
   char* massA, *massC;
   // name of diffusant species
   char *diffA;

  // maximum number of layers to be used for chemical reactions, currently 3
  const int nmaxlayers_;
  // number of active layers starts with 3, and reduces if a layer is depleted
  int layers_;

  // relative radius below which layers are neglected
  const double rmin_;

  // gas-phase properties
  char *speciesA, *speciesC;
  double Runiv;

  class FixPropertyAtom *fix_concA_, *fix_concC_;               //  concentration of species A and C [as mass fractions]
  class FixPropertyAtom *fix_changeOfA_, *fix_changeOfC_;       //  change of concentration of species A and C [as mass per volume and time]
  class FixPropertyAtom *fix_rhogas_;                           //  density of gas
  class FixPropertyAtom *fix_tgas_;                             //  temperature of gas
  class FixPropertyAtom *fix_reactionHeat_;                     //  DeltaQ
  class FixPropertyAtom *fix_totalmole_;                        // total mole in volume
  class FixPropertyAtom *fix_pressure_;
  class FixPropertyAtom *fix_diffcoeff_;

  // molar masses of gas species
  double molMass_A_, molMass_C_;

  // handle names
  double *changeOfA_, *changeOfC_, *rhogas_, *concA_, *concC_, *T_, *reactionHeat_, *pgas_, *N_, *dCoeff_;

  // particle properties
  // these are defined as vectors with the number of components corresponding to the number of active layers
  class FixPropertyAtom *fix_layerRelRad_;
  class FixPropertyAtom *fix_layerMass_;
  // molar masses and densities do not differ from particle to particle within a species, hence they are global properties
  class FixPropertyGlobal *fix_dens_, *fix_molMass_;
  // reaction properties
  // for each reaction type (e.g. CO + ore particle), global vectors containing reaction parameters have to be defined
  class FixPropertyGlobal *fix_k0_;
  class FixPropertyGlobal *fix_Ea_;

  double *radius_;                                  // radius of particle
  double **relRadii_;                               // relative radii
  double **radLayer_;
  double **massLayer_;
  double *pmass_;                                   // particle mass
  double *pdensity_;
  const double *layerDensities_, *layerMolMasses_;
  const double *k0_, *Ea_;
};
}

#endif
#endif

