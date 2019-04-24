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
#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixChemShrinkCore : public Fix  {

public:
  FixChemShrinkCore(class LAMMPS *, int, char **);
  ~FixChemShrinkCore();

  void post_create();
  void pre_delete(bool unfixflag);
  int setmask();

  virtual void updatePtrs();
  virtual void init();
  void init_defaults();
  virtual void post_force(int);

 protected:

  // functions declared in this class
  int active_layers(int);   // calculate number of active layers per-particle
  void calcMassLayer(int);  // calculate mass of layers per-particle
  void FractionalReduction(int, double *); // calculate fractional reduction per-layer depending on layer radius
  void getXi(int, double *);    // calculate molar equilibrium constant of reacting gas
  double K_eq(int, int); // calculate equilibrium constant based on the work of Valipour 2009
  void getA(int);   // calculate chemical reaction resistance term
  void getB(int);   // calculate diffusion resistance term
  void getMassT(int);   // calculate gas film mass transfer resistance term
  void reaction(int, double *, double *);   // calculate chemical reaction rate
  void update_atom_properties(int, double *, double *, double *, double *);   // update particle layers with depending on chemical reaction rate - per-particle
  void update_gas_properties(int, double *);    // update reactant and product gas masses depending on chemical reaction rate
  void heat_of_reaction(int, double *, double *, double *, double *);
  double conv_enthalpy(double *, double , int);
  double K_eq_low(int, int);
  void reaction_low(int, double *, double *);
  void FR_low(int, double *);
  void getXi_low(int, double *);
  void getA_low(int);

  // pre-defined variables for reduction process
  static int const nmaxlayers_ = 3;

  // variables
  bool screenflag_;
  double TimeStep;
  char* massA, *massC;

  double molMass_A_, molMass_C_;
  double scale_reduction_rate;

  char *diffA;
  char *moleFracA, *moleFracC;

  int layers_;          // current active layers
  double minMolarFrac_;
  const double rmin_;   // relative radius below which layers are neglected
  char *speciesA, *speciesC;

  // particle-layer variable values
  double **rhoeff_;
  double **porosity_;
  double pore_diameter_;
  double tortuosity_;
  double **relRadii_;
  double **massLayer_;
  double **Ea_, **k0_;
  //const double *k0_, *Ea_;

  // particle propertis
  double *radius_;
  double *pmass_;
  double *pdensity_;

  // handles of fixes
  double *changeOfA_, *changeOfC_, *T_, *molecularDiffusion_, *nuf_, *Rep_, *partP_, *Massterm, *reactionHeat_;
  double **Aterm, **Bterm, **effDiffBinary, **effDiffKnud, **fracRed_;

  // coarse_graining factor
  double cg_;

  class FixPropertyAtom *fix_changeOfA_, *fix_changeOfC_;
  class FixPropertyAtom *fix_tgas_;
  class FixPropertyAtom *fix_reactionHeat_;
  class FixPropertyAtom *fix_diffcoeff_;
  class FixPropertyAtom *fix_nuField_;
  class FixPropertyAtom *fix_partRe_;

  class FixPropertyAtom *fix_moleFractionA_, *fix_moleFractionC_;
  double *xA_, *xC_;

  class FixPropertyAtom *fix_fracRed;
  class FixPropertyAtom *fix_Aterm;
  class FixPropertyAtom *fix_Bterm;
  class FixPropertyAtom *fix_Massterm;
  class FixPropertyAtom *fix_effDiffBinary;
  class FixPropertyAtom *fix_effDiffKnud;
  class FixPropertyAtom *fix_partPressure_;


  // particle properties
  class FixPropertyAtom *fix_layerRelRad_;
  class FixPropertyAtom *fix_layerMass_;
  class FixPropertyAtom *fix_k0_;
  class FixPropertyAtom *fix_Ea_;
  class FixPropertyAtom *fix_porosity_;
  class FixPropertyAtom *fix_rhoeff_;
  class FixPropertyGlobal *fix_tortuosity_;
  class FixPropertyGlobal *fix_pore_diameter_;

  class FixPropertyAtom *fix_totalMole_;
  double *molarConc_;

  class FixPropertyAtom *fix_dY_;
  double **dY;

  class FixPropertyAtom *fix_dmA_;
  double **dmA_f_;

};
}

#endif
#endif

