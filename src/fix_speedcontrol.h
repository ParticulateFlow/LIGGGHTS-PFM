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

this fix controls the velocity of the particles
based on the fix addforce
may need some cleaning

contributing authors
Gerhard Holzinger (K1MET)
---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(speedcontrol,FixSpeedControl)

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
  double memory_usage();

 private:
  double xvalue,yvalue,zvalue;
  int varflag,iregion;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle;
  int nlevels_respa;

  int maxatom;
  double **sforce;

  double K;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix speedcontrol does not exist

Self-explanatory.

E: Variable name for fix speedcontrol does not exist

Self-explanatory.

E: Variable for fix speedcontrol is invalid style

Self-explanatory.

*/
