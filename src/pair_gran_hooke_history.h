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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gran/hooke/history,PairGranHookeHistory)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_H

#include "pair_gran.h"

namespace LAMMPS_NS {

class PairGranHookeHistory : public PairGran {

 friend class FixWallGranHookeHistory;
 friend class FixCheckTimestepGran;

 public:

  PairGranHookeHistory(class LAMMPS *);
  ~PairGranHookeHistory();

  virtual void settings(int, char **);
  virtual void init_granular(); //NP makes inits specific to certain gran style     //NP modified C.K.

  virtual void compute_force(int eflag, int vflag, int addflag);

  template <int ROLLINGFRICTION>
  void compute_force_eval(int eflag, int vflag, int addflag);

  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

 protected:

  virtual void history_args(char**);
  void allocate_properties(int);

  bool forceoff()
  { return force_off; }

  //NP modified C.K.

  //NP define the mechanical properties for all models here, since other model classes are derived from
  //NP the PairGranHookeHistory class

  //NP these are properties defined for each atom type (each material) - thus vectors
  class FixPropertyGlobal* Y1; //Youngs Modulus
  class FixPropertyGlobal* v1; //Poisson's ratio
  class FixPropertyGlobal* cohEnergyDens1; //Cohesion energy density

  class FixPropertyGlobal* coeffMu1; // Fluid viscosity
  class FixPropertyGlobal* coeffRestMax1;  // Maximum restitution coefficient (for mu=0)
  class FixPropertyGlobal* coeffStc1; // Critical Stokes number (10-30 for glass beads)

  //NP these are properties defined for each pair of atom types (pair of materials) - thus matrices
  class FixPropertyGlobal* coeffRest1; //coefficient of restitution
  class FixPropertyGlobal* coeffFrict1; //coefficient of (static) friction
  class FixPropertyGlobal* coeffRollFrict1; //characteristic velocity needed for Linear Spring Model
  class FixPropertyGlobal* coeffRollVisc1; //coefficient of rolling viscous damping (epsd model)

  //NP these are general properties defined only once - thus scalars
  int charVelflag;
  class FixPropertyGlobal* charVel1; //characteristic velocity needed for Linear Spring Model

  //NP here, pre-calculated contact parameters for all possible material combinations
  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict;
  double charVel, **coeffRollFrict,**coeffRollVisc,**coeffMu,**coeffRestMax,**coeffStc;

  //NP these values are now calculated from the material properties for each contact, they are thus not constant!
  //NP double kn,kt,gamman,gammat,xmu;

  virtual void deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu,double &vnnr);
  virtual void addCohesionForce(int &, int &,double &,double &);

  int cohesionflag; //NP indicates if linear cohesion model is used
  int dampflag,rollingflag,viscousflag; //NP indicates if tang damping or rolling friction is used

  // option to turn off all force computations
  bool force_off;

  // option to do some sanity checks
  bool sanity_checks;
};

}

#endif
#endif
