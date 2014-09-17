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

FixStyle(nve/sphere/aos/omp,FixNVESphereAOSOMP)

#else

#ifndef LMP_FIX_NVE_SPHERE_AOS_OMP_H
#define LMP_FIX_NVE_SPHERE_AOS_OMP_H

#include "fix_nve_sphere.h"

namespace LAMMPS_NS {

class FixNVESphereAOSOMP : public FixNVESphere {
 public:
  FixNVESphereAOSOMP(class LAMMPS *lmp, int narg, char **arg) :
    FixNVESphere(lmp, narg, arg) {};

  virtual void initial_integrate(int);
  virtual void final_integrate();
};

}

#endif
#endif
