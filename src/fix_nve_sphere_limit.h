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

#ifdef FIX_CLASS

FixStyle(nve/sphere/limit,FixNVESphereLimit)

#else

#ifndef LMP_FIX_NVE_SPHERE_LIMIT_H
#define LMP_FIX_NVE_SPHERE_LIMIT_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVESphereLimit : public FixNVE {
 public:
  FixNVESphereLimit(class LAMMPS *, int, char **);
  ~FixNVESphereLimit() {}
  void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void reset_dt();
  double compute_scalar();

 private:
  int extra, ncount;
  double vlimit,vlimitsq, omegalimit,omegalimitsq;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere requires atom style sphere

Self-explanatory.

E: Fix nve/sphere requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere requires extended particles

This fix can only be used for particles of a finite size.

*/
