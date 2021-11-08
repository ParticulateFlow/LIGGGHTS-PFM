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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(chem/shrink/core/single,FixChemShrinkCoreSingle)

#else

#ifndef LMP_FIX_CHEM_SHRINKCORESINGLE_H
#define LMP_FIX_CHEM_SHRINKCORESINGLE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChemShrinkCoreSingle : public Fix  {

public:
  FixChemShrinkCoreSingle(class LAMMPS *, int, char **);
  ~FixChemShrinkCoreSingle();

  void post_create();
  void pre_delete(bool unfixflag);
  int setmask();

  virtual void updatePtrs();
  virtual void init();
  void init_defaults();
  virtual void setup(int);
  virtual void post_force(int);
  virtual void update_fix(int, char **);

 protected:

  int active_layers(int);   // calculate number of active layers per-particle
  void calcMassLayer(int);  // calculate mass of layers per-particle
  void FractionalReduction(int); // calculate fractional reduction per-layer depending on layer radius
  void getXi(int, double &);    // calculate molar equilibrium constant of reacting gas
  double K_eq(int); // calculate equilibrium constant
  void getA(int);   // calculate chemical reaction resistance term
  void getB(int);   // calculate diffusion resistance term
  void getMassT(int);   // calculate gas film mass transfer resistance term
  void reaction(int, double &, const double);   // calculate chemical reaction rate
  // update particle layers (relative radii, mass) depending on chemical reaction rate
  void update_atom_properties(int, const double);
  // update reactant and product gas masses depending on chemical reaction rate
  void update_gas_properties(int, const double);
  void heat_of_reaction(int, const double);

  // pre-defined variables for reduction process
//  int const nmaxlayers_;
  static double const Runiv; // universal gas constant

  // variables
  bool screenflag_;
  double TimeStep;
  char* massA, *massC;

  double molMass_A_, molMass_B_, molMass_C_;
  int nu_A_, nu_B_, nu_C_;   
  double scale_reduction_rate;

  char *diffA;
  char *moleFracA, *moleFracC;

  int layers_;          // current active layers
  double minMolarFrac_;
  double rmin_;   // radius below which layers are neglected
  char *speciesA, *speciesC;
  bool layerDiffusion_;
  bool heatToFluid_;
  bool heatToParticle_;
  bool shrink_;

  int reactionHeatIndex_;
  int KeqIndex_;

  // particle-layer variable values
  double **rhoeff_;
  double *porosity_;
  double pore_diameter_;
  double tortuosity_;
  double **relRadii_; // relative radii of individual layers
  double **massLayer_; // mass of individual layers
  double *effvolfactors_;
  double k0_; // frequency factor
  double T0_; // activation temperature
  double Tmin_; // minimum temperature for reaction
  int nPreFactor_;

  // particle propertis
  double *radius_;
  double *pmass_;
  double *pdensity_;

  // handles of fixes
  double *changeOfA_, *changeOfC_;
  double *T_;
  double *Tpart_;
  double *molecularDiffusion_;
  double *nuf_;
  double *Rep_;
  double *partP_;
  double *heatFlux_;       // heat flux to/from grain
  double *Massterm;         // mass transfer resistance
  double *Aterm;           // reaction resistance
  double *Bterm;           // diffusion resistance
  double *effDiffBinary;   // effective binary diffusion coefficient
  double *effDiffKnud;     // effective Knudsen diffusion coefficient
  double **fracRed_;        // fractional reduction
#ifdef PER_ATOM_LAYER_DENSITIES
  double **layerDensities_;
#else
  double *layerDensities_;
#endif

  // coarse_graining factor
  double cg_;

  class FixPropertyAtom *fix_changeOfA_;    // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_changeOfC_;    // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_tgas_;         // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_tpart_;        // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_heatFlux_;     // [cfd/coupling/convection]
  class FixPropertyAtom *fix_diffcoeff_;    // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_nuField_;      // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_partRe_;       // [cfd/coupling/chemistry]

  class FixPropertyAtom *fix_moleFractionA_; // [cfd/coupling/chemistry]
  class FixPropertyAtom *fix_moleFractionC_; // [cfd/coupling/chemistry]
  double *xA_, *xC_;

  class FixPropertyAtom *fix_fracRed;       // [script]
  class FixPropertyAtom *fix_Aterm;         // [internal]
  class FixPropertyAtom *fix_Bterm;         // [internal]
  class FixPropertyAtom *fix_Massterm;      // [internal]
  class FixPropertyAtom *fix_effDiffBinary; // [internal]
  class FixPropertyAtom *fix_effDiffKnud;   // [internal]
  class FixPropertyAtom *fix_partPressure_; // [cfd/coupling/chemistry]


  // particle properties
  class FixPropertyAtom *fix_layerRelRad_;  // [script]
  class FixPropertyAtom *fix_layerMass_;    // [script]

#ifdef PER_ATOM_LAYER_DENSITIES
  class FixPropertyAtom *fix_layerDens_;
#else
  class FixPropertyGlobal *fix_layerDens_;  // [script]
#endif

  class FixPropertyAtomPolydispParcel *fix_polydisp_;

  class FixPropertyAtom *fix_rhoeff_;       // [script]
  class FixPropertyGlobal *fix_porosity_;     // [script/internal]
  class FixPropertyGlobal *fix_tortuosity_; // [script]
  class FixPropertyGlobal *fix_pore_diameter_; // [script]

  class FixPropertyAtom *fix_dY_; // [internal]
  double *dY;

  class FixPropertyAtom *fix_dmA_; // [internal]
  double *dmA_f_;

  class FixPropertyAtom *fix_reactionheat;
  double *reactionheat_;

  class FixPropertyAtom *fix_reactantPerParticle_;
  double *reactantPerParticle_;
  bool limit_reactant_consumption_;
  double maxReactantConsumptionFrac_;

  // constant parameters for reactions
  const double Cp_coke_; // heat capacity coke in J/(mol K)
  const double T_room_;

};
}

#endif
#endif

