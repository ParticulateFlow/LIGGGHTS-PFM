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
#define NDEBUG
#include <assert.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

/*NL*/ #define DEBUGMODE_LMP_FIX_NEIGHLIST_MESH false //(update->ntimestep>15400 && comm->me ==1)
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID 0
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID 60

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

//NP this is called before FixContactHistoryMesh::pre_force()
//NP because of the order of creation in FixWallGran::post_create()
//NP this is important so building new neigh list before refreshing contact hist

void FixNeighlistMeshOMP::pre_force(int vflag)
{
    if(!buildNeighList) return;

    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    /*NL*/ //fprintf(screen,"***building neighbor list at timestep "BIGINT_FORMAT"\n",update->ntimestep);

    buildNeighList = false;

    const int nthreads = comm->nthreads;

    numAllContacts_ = 0;

    // copy current to old # of neighbors
    //NP this is important for correct contact history copy
    memset(fix_nneighs_->vector_atom, 0, sizeof(double)*atom->nlocal);

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    //NP cutneighmax includes contactDistanceFactor, thus rmax includes this as well
    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    double prev_skin = skin;
    double prev_distmax = distmax;

    if(changingMesh)
    {
      skin = neighbor->skin;
      //NP cutneighmax includes contactDistanceFactor, thus distmax includes this as well
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      //NP cutneighmax includes contactDistanceFactor, thus distmax includes this as well
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    //NP from neighbor; get the binning
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

    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID <= atom->get_map_size() && update->ntimestep > 0 &&
    /*NL*/      atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID) >= 0)
    /*NL*/ {
    /*NL*/          if(!r) error->one(FLERR,"debugmode not made for SPH");
    /*NL*/          int iTriDeb = mesh_->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID);
    /*NL*/          int iAtomDeb = atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID);
    /*NL*/          int ixDeb, iyDeb, izDeb;
    /*NL*/          int iBinDeb = neighbor->coord2bin(atom->x[iAtomDeb],ixDeb, iyDeb, izDeb);
    /*NL*/          fprintf(screen, "**step "BIGINT_FORMAT", particle id %d at bin %d (indixes %d %d %d) on proc %d, within skin to target tri %s\n",
    /*NL*/                      update->ntimestep,DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID,
    /*NL*/                      iBinDeb,ixDeb, iyDeb, izDeb,comm->me,
    /*NL*/                      mesh_->resolveTriSphereNeighbuild(iTriDeb,atom->radius[iAtomDeb]*neighbor->contactDistanceFactor,atom->x[iAtomDeb],skin) ? "true" : "false" );
    /*NL*/ }

  const bool use_parallel = nall > static_cast<size_t>(nthreads*nthreads);
  const bool load_balance = true;
  const bool sort = true;

  if(use_parallel) {
    if(load_balance) {
      // DEPRECATED: strided access version, TODO: replace with Zoltan partitioning
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
      printf("%ld - %s, %s: sorting %lu triangles...\n", update->ntimestep, id, style, nall);
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

    // NL printf("[%d] %d / %d particles in contact with this wall\n", comm->me, ncontacts, nlocal);

    if(globalNumAllContacts_) {
      MPI_Sum_Scalar(numAllContacts_,world);
    }

    return;
  }

  // NL printf("[%d] %d / %d particles in contact with this wall\n", comm->me, ncontacts, nlocal);

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
    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
    /*NL*/ {
    /*NL*/          double lo[3],hi[3];
    /*NL*/          //b.getBoxBounds(lo,hi);
    /*NL*/          vectorScalarSubtract3D(lo,distmax);
    /*NL*/          vectorScalarAdd3D(hi,distmax);
    /*NL*/          int ixDebLo, iyDebLo, izDebLo, ixDebHi, iyDebHi, izDebHi;
    /*NL*/          neighbor->coord2bin(lo,ixDebLo, iyDebLo, izDebLo);
    /*NL*/          neighbor->coord2bin(hi,ixDebHi, iyDebHi, izDebHi);
    /*NL*/          fprintf(screen, "handleTriangle for tri id %d on proc %d, indices from here  are %d %d //  %d %d //  %d %d , bbox is %f %f //  %f %f //  %f %f  \n",
    /*NL*/                      mesh_->id(iTri),comm->me,ixMin,ixMax,iyMin,iyMax,izMin,izMax,lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
    /*NL*/          fprintf(screen, "handleTriangle for tri id %d on proc %d, indices from neigh are %d %d //  %d %d //  %d %d , bbox is %f %f //  %f %f //  %f %f  \n",
    /*NL*/                      mesh_->id(iTri),comm->me,ixDebLo,ixDebHi,iyDebLo,iyDebHi,izDebLo,izDebHi,lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
    /*NL*/ }

        for(int ix=ixMin;ix<=ixMax;ix++) {
          for(int iy=iyMin;iy<=iyMax;iy++) {
            for(int iz=izMin;iz<=izMax;iz++) {
              int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;

              /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
              /*NL*/          fprintf(screen, "       handleTriangle tri id %d on proc %d - checking bin %d\n",
              /*NL*/                      mesh_->id(iTri),comm->me, iBin);

              int iAtom = binhead[iBin];
              //NP only handle local atoms
              while(iAtom != -1 && iAtom < nlocal)
              {
                if(! (mask[iAtom] & groupbit))
                {
                    if(bins) iAtom = bins[iAtom];
                    else iAtom = -1;
                    continue;
                }
                nchecked++;

                /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri)
                /*NL*/                                     && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID == atom->map(iAtom) )
                /*NL*/          fprintf(screen, "   handleTriangle atom %d for tri id %d on proc %d\n",
                /*NL*/                      atom->map(iAtom),mesh_->id(iTri),comm->me);
                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  //NP include iAtom in neighbor list
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
          /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
          /*NL*/          fprintf(screen, "       handleTriangle tri id %d on proc %d - checking bin %d\n",
          /*NL*/                      mesh_->id(iTri),comm->me, iBin);

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

            /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri)
            /*NL*/                                     && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID == atom->map(iAtom) )
            /*NL*/          fprintf(screen, "   handleTriangle atom %d for tri id %d on proc %d\n",
            /*NL*/                      atom->map(iAtom),mesh_->id(iTri),comm->me);
            if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
            {
              //NP include iAtom in neighbor list
              neighbors.push_back(iAtom);
              // Incrementing fix_nneighs is a data race, fix this after parallel region!
            }
            if(bins) iAtom = bins[iAtom];
            else iAtom = -1;
          }
        }
      }
    }

    /*NL*/// fprintf(screen,"iTri %d numContacts %d\n",iTri, neighbors.size());
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
