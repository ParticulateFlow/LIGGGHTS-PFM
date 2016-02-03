/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2015- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Richard Berger <richard.berger@jku.at> (JKU Linz)
------------------------------------------------------------------------- */
#ifdef COMPUTE_CLASS
ComputeStyle(wall/ghosts,ComputeWallGhosts)
#else

#ifndef LMP_COMPUTE_WALL_GHOSTS_H
#define LMP_COMPUTE_WALL_GHOSTS_H

#include "compute.h"
#include "fix_wall_gran.h"
#include <string>

namespace LAMMPS_NS {

class ComputeWallGhosts : public Compute {
  std::string wall_fix_id;
  FixWallGran * wall;

 public:
  ComputeWallGhosts(class LAMMPS *, int, char **);
  ~ComputeWallGhosts();
  void init();
  virtual double compute_scalar();
};

}

#endif
#endif
