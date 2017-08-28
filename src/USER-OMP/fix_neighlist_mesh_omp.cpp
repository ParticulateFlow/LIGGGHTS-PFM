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
   Contributing authors:
   Richard Berger (JKU Linz)
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh_omp.h"
#include "fix_mesh_surface_omp.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "atom.h"
#include "comm.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>
#include "thr_omp.h"
#include <assert.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

FixNeighlistMeshOMP::FixNeighlistMeshOMP(LAMMPS *lmp, int narg, char **arg)
: FixNeighlistMesh(lmp,narg,arg)
{
}

/* ---------------------------------------------------------------------- */

FixNeighlistMeshOMP::~FixNeighlistMeshOMP()
{
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMeshOMP::initializeNeighlist()
{
  FixNeighlistMesh::initializeNeighlist();
  int nall = mesh_->sizeLocal()+mesh_->sizeGhost();

  trilist.clear();
  trilist.reserve(nall);

  for(int iTri = 0; iTri < nall; iTri++) {
    trilist.push_back(iTri);
  }

  sortfreq = 10000;
  nextsort = 1000;
}

/* ---------------------------------------------------------------------- */

int * FixNeighlistMeshOMP::partition_begin(int tid) {
  return &partition_local_indices[0] + thread_offsets[tid];
}

int * FixNeighlistMeshOMP::partition_end(int tid) {
  return &partition_local_indices[0] + thread_offsets[tid+1];
}

bool FixNeighlistMeshOMP::in_thread_partition(int tid, int i) {
  return partition_global_thread[i] == tid;
}

/* ---------------------------------------------------------------------- */

struct TriangleComparator {
  std::vector<int> & triangles;

  TriangleComparator(std::vector<int> & triangles) : triangles(triangles)
  {
  }

  bool operator() (int i,int j) {
    return triangles[i] > triangles[j];
  }
};

struct ThreadComparator {
  std::vector<int> & thread;

  ThreadComparator(std::vector<int> & thread) : thread(thread)
  {
  }

  bool operator() (int i, int j) {
    return thread[i] < thread[j];
  }
};

/* ---------------------------------------------------------------------- */

void FixNeighlistMeshOMP::pre_force(int vflag)
{
    if(!buildNeighList) return;

    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    buildNeighList = false;

    const int nthreads = comm->nthreads;

    numAllContacts_ = 0;

    // copy current to old # of neighbors
    memset(fix_nneighs_->vector_atom, 0, sizeof(double)*atom->nlocal);

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    double prev_skin = skin;
    double prev_distmax = distmax;

    if(changingMesh)
    {
      skin = neighbor->skin;
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;

    const size_t nall = mesh_->sizeLocal() + mesh_->sizeGhost();

    // update cache if necessary
    if (triangles.size() != nall) {
      initializeNeighlist();
    }

    // update precomputed bins if necessary
    if((skin != prev_skin) || (distmax != prev_distmax) || (neighbor->last_setup_bins_timestep > last_bin_update)) {
      generate_bin_list(nall);
    }

  const bool use_parallel = nall > static_cast<size_t>(nthreads*nthreads);
  const bool load_balance = true;
  const bool sort = true;

  if(use_parallel) {
    if(load_balance) {
      #if defined(_OPENMP)
      #pragma omp parallel default(none)
      #endif
      {
        int twidth, tstride, ioffset, tid;
        loop_setup_thr(twidth, tstride, ioffset, tid, nall, nthreads);

        size_t i = ioffset;
        int offset = ioffset;
        int col = 0;
        while(i < nall) {
          const int iTri = trilist[i];
          handleTriangle(iTri);

          col++;
          if (col == twidth) {
            offset += tstride;
            col = 0;
          }
          i = offset + col;
        }
      }
    } else {
      #if defined(_OPENMP)
      #pragma omp parallel default(none)
      #endif
      {
        int ifrom, ito, tid;
        loop_setup_thr(ifrom, ito, tid, nall, nthreads);

        for(int i = ifrom; i < ito; i++) {
          const int iTri = trilist[i];
          handleTriangle(iTri);
        }
      }
    }

    if(sort && update->ntimestep >= nextsort) {
      sort_triangles_by_nchecked();
      nextsort = (update->ntimestep/sortfreq)*sortfreq + sortfreq;
    }
  } else {
    for(size_t i = 0; i < nall; i++) {
      const int iTri = trilist[i];
      handleTriangle(iTri);
    }
  }

  // prepare memory for partition generation
  const int nlocal = atom->nlocal;

  thread_offsets.resize(nthreads+1);
  partition_local_triangles.clear();
  partition_global_indices.clear();
  partition_global_triangles.resize(nlocal);
  partition_global_thread.resize(nlocal);

  std::fill_n(partition_global_triangles.begin(), nlocal, 0);
  std::fill_n(partition_global_thread.begin(), nlocal, -1);

  std::fill(thread_offsets.begin(), thread_offsets.end(), 0);

  // update nneighs
  for(size_t iTri = 0; iTri < nall; ++iTri) {
    std::vector<int> & neighbors = triangles[iTri].contacts;
    for(std::vector<int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
      const int i =  *it;
      ++partition_global_triangles[i];
    }
    numAllContacts_ += neighbors.size();
  }

  for(int i = 0; i < nlocal; ++i) {
    const int ntriangles = partition_global_triangles[i];
    if(ntriangles > 0) {
      partition_global_indices.push_back(i);
      partition_local_triangles.push_back(ntriangles);
      fix_nneighs_->set_vector_atom_int(i, ntriangles);
    }
  }

  const int ncontacts = partition_global_indices.size();

  if(ncontacts == 0) {
    partition_local_indices.push_back(0);

    if(globalNumAllContacts_) {
      MPI_Sum_Scalar(numAllContacts_,world);
    }

    return;
  }

  partition_local_indices.resize(ncontacts);
  partition_local_thread.resize(ncontacts);

  for(int i = 0; i < ncontacts; ++i) {
    partition_local_indices[i] = i;
  }

  // ======================================================================================================
  // generate particle partitions for mesh (inefficiently for a start)

  // 1. sort particles based on triangle and existing thread assignment
  std::sort(partition_local_indices.begin(), partition_local_indices.end(), TriangleComparator(partition_local_triangles));

  // 2. assign each index to a thread
  for(int i = 0; i < ncontacts; ++i) {
    const int tid = i % nthreads;
    partition_local_thread[partition_local_indices[i]] = tid;
    ++thread_offsets[tid+1]; // count number of particles per thread
  }

  // 3. collect partition indices
  std::stable_sort(partition_local_indices.begin(), partition_local_indices.end(), ThreadComparator(partition_local_thread));

  // 4. set thread boundaries
  for(int tid = 1; tid <= nthreads; ++tid) {
    thread_offsets[tid] += thread_offsets[tid-1];
  }

  //  thread_mapping.clear();

  // store final permutation in partition_local_indices
  for(int i = 0; i < ncontacts; ++i) {
    int local_idx = partition_local_indices[i];
    const int global_idx = partition_global_indices[local_idx];
    partition_global_thread[global_idx] = partition_local_thread[local_idx];
    partition_local_indices[i] = global_idx;
  }

//#define DEBUG
#ifdef DEBUG
  for(int i = 0; i < ncontacts; ++i) {
    const int global_idx = partition_local_indices[i];
    partition_local_thread[i] = partition_global_thread[global_idx];
  }

  for(int tid = 0; tid < nthreads; ++tid) {
    int ifrom = thread_offsets[tid];
    int ito = thread_offsets[tid+1];

    for(int i = ifrom; i < ito; ++i) {
      const int global_idx = partition_local_indices[i];
      assert(partition_global_thread[global_idx] == tid);
    }
  }
#endif

  // ======================================================================================================

  if(globalNumAllContacts_) {
    MPI_Sum_Scalar(numAllContacts_,world);
  }
}

/* ---------------------------------------------------------------------- */

// Replace FixNeighlistMesh::handleTriangle (no override, static linking!)
// Removes increments of nneighs, this is instead done serially at the end of pre_force

void FixNeighlistMeshOMP::handleTriangle(int iTri)
{
    TriangleNeighlist & triangle = triangles[iTri];
    std::vector<int> & neighbors = triangle.contacts;
    int & nchecked = triangle.nchecked;
    int *mask = atom->mask;
    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    const double contactDistanceFactor = neighbor->contactDistanceFactor;

    neighbors.clear();

    nchecked = 0;

    // only do this if I own particles
    if(nlocal)
    {
      if(changingMesh || changingDomain)
      {
        getBinBoundariesForTriangle(iTri,ixMin,ixMax,iyMin,iyMax,izMin,izMax);

        for(int ix=ixMin;ix<=ixMax;ix++) {
          for(int iy=iyMin;iy<=iyMax;iy++) {
            for(int iz=izMin;iz<=izMax;iz++) {
              int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;

              int iAtom = binhead[iBin];
              while(iAtom != -1 && iAtom < nlocal)
              {
                if(! (mask[iAtom] & groupbit))
                {
                    if(bins) iAtom = bins[iAtom];
                    else iAtom = -1;
                    continue;
                }
                nchecked++;

                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  neighbors.push_back(iAtom);
                  // Incrementing fix_nneighs is a data race, fix this after parallel region!
                }
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
              }
            }
          }
        }
      } else {
        const std::vector<int> & triangleBins = triangle.bins;
        const int bincount = triangleBins.size();
        for(int i = 0; i < bincount; i++) {
          const int iBin = triangleBins[i];

          int iAtom = binhead[iBin];
          while(iAtom != -1 && iAtom < nlocal)
          {
            if(! (mask[iAtom] & groupbit))
            {
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
                continue;
            }
            nchecked++;

            if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
            {
              neighbors.push_back(iAtom);
              // Incrementing fix_nneighs is a data race, fix this after parallel region!
            }
            if(bins) iAtom = bins[iAtom];
            else iAtom = -1;
          }
        }
      }
    }
}

/* ---------------------------------------------------------------------- */

template<typename X>
inline void swap(X& a, X& b) {
  X tmp = a;
  a = b;
  b = tmp;
}

int FixNeighlistMeshOMP::quicksort_partition(int p, int r) {
  int x = trilist[r];
  int i = p - 1;

  for (int j = p; j < r; j++) {
    TriangleNeighlist & triangle_j = triangles[trilist[j]];
    TriangleNeighlist & triangle_x = triangles[x];
    if (triangle_j.nchecked > triangle_x.nchecked) {
      i++;

      // exchange
      swap(trilist[i], trilist[j]);
    }
  }
  swap(trilist[i+1], trilist[r]);
  return i+1;
}

void FixNeighlistMeshOMP::quicksort(int p, int r) {
  if (p < r) {
    int q = quicksort_partition(p, r);
    quicksort(p, q-1);
    quicksort(q+1, r);
  }
}

void FixNeighlistMeshOMP::sort_triangles_by_nchecked() {
  const int nall = mesh_->sizeLocal()+mesh_->sizeGhost();
  quicksort(0, nall-1);
}
