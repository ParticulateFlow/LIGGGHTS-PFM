/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2017- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_forcecontrol_region_universe.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_extra_liggghts.h"
#include "mpi_liggghts.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;

/* ---------------------------------------------------------------------- */

FixForceControlRegionUniverse::FixForceControlRegionUniverse(LAMMPS *lmp, int narg, char **arg) :
  FixForceControlRegion(lmp, narg, arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition for fix forcecontrol/region/universe");

  target_stress_ = NULL;
  target_v_ave_ = NULL;
  target_v_min_ = NULL;
  target_v_max_ = NULL;

  int iarg = 3;
  while (iarg < narg) {
    if(strcmp(arg[iarg],"target_val") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'target_val'");
      idtarget_hash = bitwiseHash(arg[iarg+1], mpi_tag_upper_bound(universe->uworld));
      iarg += 2;
    } else if (strcmp(arg[iarg],"receive_from_partition") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      receive_from_world_ = atoi(arg[iarg+1])-1;
      if(receive_from_world_ < 0 || receive_from_world_ >= universe->nworlds)
        error->fix_error(FLERR,this,"receive_from_world_ must be a valid communicator");
      iarg += 2;
    } else if (strcmp(arg[iarg],"couple_every") == 0 || strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      couple_every_ = atoi(arg[iarg+1]);
      if(couple_every_ < 0) error->fix_error(FLERR,this,"couple_every must be >= 0");
      nevery = couple_every_;
      iarg += 2;
    } else {
      ++iarg;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixForceControlRegionUniverse::~FixForceControlRegionUniverse()
{
  memory->destroy(target_stress_);
  memory->destroy(target_v_ave_);
  memory->destroy(target_v_min_);
  memory->destroy(target_v_max_);
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegionUniverse::receive_post_create_data()
{
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&cg_target_, 1, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&cg_target_, 1, MPI_DOUBLE, 0, world);

  int ncells = 0;
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&ncells, 1, MPI_INT, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&ncells, 1, MPI_INT, 0, world);

  std::vector<int> cellids(ncells,0);
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&cellids[0], ncells, MPI_INT, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&cellids[0], ncells, MPI_INT, 0, world);

  for(int i=0; i<ncells; ++i)
    target_cellid2index_[cellids[i]] = i;

  if (ncells > 0) {
    target_count_.resize(ncells, 0);
    target_mass_.resize(ncells, 0.);
    target_vol_fr_.resize(ncells, 0.);
    memory->grow(target_stress_, ncells, 7, "forcecontrol:target_stress_");
    std::fill_n(&target_stress_[0][0],   7*ncells, 0.);
    memory->grow(target_v_ave_,  ncells, 3, "forcecontrol:target_v_ave_");
    std::fill_n(&target_v_ave_[0][0],    3*ncells, 0.);
    memory->grow(target_v_min_,  ncells, 3, "forcecontrol:target_v_min_");
    std::fill_n(&target_v_min_[0][0],    3*ncells, 0.);
    memory->grow(target_v_max_,  ncells, 3, "forcecontrol:target_v_max_");
    std::fill_n(&target_v_max_[0][0],    3*ncells, 0.);
  }
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegionUniverse::post_create_stress_part()
{
  double maxrd,minrd;
  modify->max_min_rad(maxrd,minrd);
  const_part_ = (2.*maxrd)*1.2;
  used_part_ = (2.*maxrd)*1.4;
  sinesq_part_ = used_part_ - const_part_;
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegionUniverse::receive_coupling_data()
{
  if((couple_every_ > 0) && (receive_from_world_ >= 0) && (update->ntimestep % couple_every_ == 0)) {
    int ncells = target_cellid2index_.size();
    if (comm->me == 0) {
      MPI_Recv(&target_count_[0],       ncells, MPI_INT,    universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_mass_[0],        ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_vol_fr_[0],      ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_stress_[0][0], 7*ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_v_ave_[0][0],  3*ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_v_min_[0][0],  3*ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&target_v_max_[0][0],  3*ncells, MPI_DOUBLE, universe->root_proc[receive_from_world_], idtarget_hash, universe->uworld, MPI_STATUS_IGNORE);
    }

    MPI_Bcast(&target_count_[0],       ncells, MPI_INT,    0, world);
    MPI_Bcast(&target_mass_[0],        ncells, MPI_DOUBLE, 0, world);
    MPI_Bcast(&target_vol_fr_[0],      ncells, MPI_DOUBLE, 0, world);
    MPI_Bcast(&target_stress_[0][0], 7*ncells, MPI_DOUBLE, 0, world);
    MPI_Bcast(&target_v_ave_[0][0],  3*ncells, MPI_DOUBLE, 0, world);
    MPI_Bcast(&target_v_min_[0][0],  3*ncells, MPI_DOUBLE, 0, world);
    MPI_Bcast(&target_v_max_[0][0],  3*ncells, MPI_DOUBLE, 0, world);
  }
}
