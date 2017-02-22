/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
contributing authors
Gerhard Holzinger (K1MET)
---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(speedControl,FixSpeedControl)

#else

#ifndef LMP_FIX_SPEEDCONTROL_H
#define LMP_FIX_SPEEDCONTROL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSpeedControl : public Fix {
 public:
  FixSpeedControl(class LAMMPS *, int, char **);
  ~FixSpeedControl();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  double memory_usage();

 private:
  double xvalue,yvalue,zvalue;
  int varflag,iregion;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle;
  double foriginal[4],foriginal_all[4];
  int force_flag;
  int nlevels_respa;

  int maxatom;
  double **sforce;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix speedControl does not exist

Self-explanatory.

E: Variable name for fix speedControl does not exist

Self-explanatory.

E: Variable for fix speedControl is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix speedControl

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix speedControl

Must define an energy vartiable when applyting a dynamic
force during minimization.

*/
