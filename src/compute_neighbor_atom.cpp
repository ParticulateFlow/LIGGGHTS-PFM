/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-2014 DCS Computing GmbH, Linz
                     Copyright 2015-     JKU Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_neighbor_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeNeighborAtom::ComputeNeighborAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute neighbor/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 1;

  nmax = 0;
  neighs = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeNeighborAtom::~ComputeNeighborAtom()
{
  memory->destroy(neighs);
}

/* ---------------------------------------------------------------------- */

void ComputeNeighborAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute neighbor/atom requires a pair style be defined");

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"neighbor/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute neighbor/atom");

  // need an occasional neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeNeighborAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeNeighborAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow neighbor array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(neighs);
    nmax = atom->nmax;
    memory->create(neighs,nmax,"neighbor/atom:neighs");
    vector_atom = neighs;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  int * const numneigh = list->numneigh;
  int * const mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  for (int i = 0; i < nall; ++i) {
    if (mask[i] & groupbit) {
      neighs[i] = numneigh[i];
    } else {
      neighs[i] = 0.0;
    }
  }

  // communicate ghost atom counts between neighbor procs if necessary

  if (force->newton_pair) comm->reverse_comm_compute(this);
}

/* ---------------------------------------------------------------------- */

int ComputeNeighborAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; ++i)
    buf[m++] = neighs[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeNeighborAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; ++i) {
    int j = list[i];
    neighs[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeNeighborAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
