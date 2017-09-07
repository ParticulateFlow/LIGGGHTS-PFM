/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   K1MET Metallurgical Competence Center
   Copyright 2017- K1MET

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Gerhard Holzinger (K1MET)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(speedcontrol,FixSpeedControl)

#else

#ifndef LMP_FIX_SPEEDCONTROL_H
#define LMP_FIX_SPEEDCONTROL_H

#include "fix.h"

namespace LAMMPS_NS {

/**
 * @brief FixSpeedControl
 *        control the velocity of particles.
 *
 * This class controls the velocity of the particles in a given region
 * and is based on FixAddForce.
 */
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
