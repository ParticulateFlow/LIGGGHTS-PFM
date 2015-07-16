/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2013-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(gravity/partitioned/omp,FixGravityPartitionedOMP)

#else

#ifndef LMP_FIX_GRAVITY_PARTITIONED_OMP_H
#define LMP_FIX_GRAVITY_PARTITIONED_OMP_H

#include "fix_gravity.h"

namespace LAMMPS_NS {

class FixGravityPartitionedOMP : public FixGravity {

 public:
  FixGravityPartitionedOMP(class LAMMPS *, int, char **);

  virtual int setmask();
  virtual void post_force(int);
};

}

#endif
#endif
