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

FixStyle(chem/shrink,FixChemShrink)

#else

#ifndef LMP_FIX_CHEM_SHRINK_H
#define LMP_FIX_CHEM_SHRINK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixChemShrink : public Fix  {

public:
  FixChemShrink(class LAMMPS *, int, char **);
  ~FixChemShrink();
  void pre_delete(bool unfixflag);


  void reaction();
  double partSurfArea(double);
  void updatePtrs();

  int setmask();
  virtual void init();
  virtual void post_force(int);
  void delete_atoms();
  void default_radius();

 protected:
  char* speciesA, *speciesC;
  char* massA, *massC;

  class FixPropertyAtom *fix_concA_, *fix_concC_;               //  concentration of species A and C
  class FixPropertyAtom *fix_changeOfA_, *fix_changeOfC_;       //  change of concentration of species A and C
  class FixPropertyAtom *fix_rhogas;                           //  density of gas
  class FixPropertyAtom *fix_tgas;                             //  temperature of gas
  class FixPropertyAtom *fix_reactionheat_;                     //  DeltaQ
  class FixPropertyAtom *fix_totalMole_;


  int iarg_;

  // values from user
  double k;                                         // reaction rate coefficient
  double molMass_A_, molMass_B_, molMass_C_;        // Molecular mass of species A, B and C
  double hertzpct;

  // particle properties
  double *radius_;                                  // radius of particle
  double *pmass_;                                   // particle mass
  double *pdensity_;

  // minimum radius value -rmin input from user
  double rmin;
  double rdefault;
  double radius_origin;
  bool rdef;

  // pointer updaters
  double *changeOfA_;
  double *changeOfC_;
  double *rhogas_;
  double *concA_;
  double *concC_;
  double *tgas_;
  double *N_;

  // while loop integers
  int spcA;
  int spcC;

  // TimeStep
  double TimeStep;

  // delete atoms lsit & shrink list
  int *dlist;

  // get values for default radius minimum
  double vmax_;
  double **Yeff_;
  const double *Y, *nu;


};

}

#endif
#endif


