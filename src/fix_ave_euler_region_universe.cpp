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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "comm.h"
#include "mpi_liggghts.h"
#include "fix_ave_euler_region_universe.h"
#include "math_extra_liggghts.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;

/* ---------------------------------------------------------------------- */

FixAveEulerRegionUniverse::FixAveEulerRegionUniverse(LAMMPS *lmp, int narg, char **arg) :
  FixAveEulerRegion(lmp, narg, arg),
  send_to_world_(-1),
  synchronize_(false)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition for fix ave/euler/region/universe");

  id_hash_ = JSHash(id, mpi_tag_upper_bound(universe->uworld));

  int iarg = 3;

  // parse args
  while(iarg < narg) {
    if (strcmp(arg[iarg],"send_to_partition") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      send_to_world_ = atoi(arg[iarg+1])-1;
      if(send_to_world_ >= universe->nworlds)
        error->fix_error(FLERR,this,"send_to_world must be a valid communicator");
      iarg += 2;
    } else if (strcmp(arg[iarg],"sync") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      if (strcmp(arg[iarg+1],"yes") == 0) synchronize_ = true;
      iarg += 2;
    } else {
      ++iarg;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixAveEulerRegionUniverse::~FixAveEulerRegionUniverse()
{
}

/* ---------------------------------------------------------------------- */

template <typename K, typename V> V get_key(const std::pair<K,V>& p) { return p.first; }

template <typename K, typename V> void mapkeys2vector(const std::map<K, V> &m, std::vector<K>& vkey)
{
  vkey.clear();
  vkey.reserve(m.size());
  std::transform(m.begin(), m.end(), std::back_inserter(vkey), get_key<K,V>);
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegionUniverse::send_post_create_data()
{
  // send information to child
  if (comm->me == 0 && send_to_world_ >= 0) {
    double cg = force->cg();
    MPI_Send(&cg, 1, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
    std::vector<int> keys;
    int ncellids = cellid2index_.size();
    mapkeys2vector(cellid2index_, keys);
    MPI_Send(&ncellids, 1, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
    MPI_Send(&keys[0], ncellids, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
  }
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegionUniverse::send_coupling_data()
{
  // check for timestep is necessary because end_of_step() is also called from setup function
  if(comm->me == 0 && send_to_world_ >= 0 && update->ntimestep % nevery == 0) {
    if (synchronize_) {
      MPI_Ssend(&ncount_[0],      ncells_, MPI_INT,    universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&mass_[0],        ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&vol_fr_[0],      ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&stress_[0][0], 7*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&v_av_[0][0],   3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&v_min_[0][0],  3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Ssend(&v_max_[0][0],  3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
    } else {
      MPI_Send(&ncount_[0],      ncells_, MPI_INT,    universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&mass_[0],        ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&vol_fr_[0],      ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&stress_[0][0], 7*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&v_av_[0][0],   3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&v_min_[0][0],  3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
      MPI_Send(&v_max_[0][0],  3*ncells_, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
    }
  }
}
