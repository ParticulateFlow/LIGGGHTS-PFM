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
#include "compute_wall_ghosts.h"
#include "lammps.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeWallGhosts::ComputeWallGhosts(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute wall/ghosts command");

  scalar_flag = 1;

  wall_fix_id = arg[3];
  wall = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeWallGhosts::~ComputeWallGhosts()
{
}

/* ---------------------------------------------------------------------- */

void ComputeWallGhosts::init()
{
  wall = dynamic_cast<FixWallGran*>(modify->find_fix_id_style(wall_fix_id.c_str(), "wall/gran"));

  if(!wall) {
    error->all(FLERR, "Could not find wall fix for compute wall/ghosts");
  }
}

/* ---------------------------------------------------------------------- */

double ComputeWallGhosts::compute_scalar()
{
  scalar = wall->n_ghosts_all();
  return scalar;
}
