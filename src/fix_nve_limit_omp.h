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
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/limit/omp,FixNVELimitOMP)

#else

#ifndef LMP_FIX_NVE_LIMIT_OMP_H
#define LMP_FIX_NVE_LIMIT_OMP_H

#include "fix_nve_limit.h"

namespace LAMMPS_NS {

class FixNVELimitOMP : public FixNVELimit {
 public:
  FixNVELimitOMP(class LAMMPS *lmp, int narg, char **arg) :
    FixNVELimit(lmp, narg, arg) {};

  virtual void initial_integrate(int);
  virtual void final_integrate();
};

}

#endif
#endif
