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

#include <cstdlib>
#include "partitioner_zoltan.h"
#include "comm.h"
#include "atom.h"
#include "error.h"
#include "domain.h"
#include <mpi.h>
#include <assert.h>
#include "force.h"
#include "update.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include <vector>
#include <algorithm>

#include <omp.h>

using namespace LAMMPS_NS;
using namespace std;

int ZoltanPartitioner::me = 0;
int ZoltanPartitioner::box_change = 0;
bool ZoltanPartitioner::neighbor_list_exists = false;

ZoltanPartitioner::ZoltanPartitioner(class LAMMPS * lmp, int argc, const char * const * argv) : Partitioner(lmp),
    zinstance(NULL),
    numGidEntries(1),
    numLidEntries(1),
    numImport(-1),
    numExport(-1),
    importGlobalIds(NULL),
    importLocalIds(NULL),
    exportGlobalIds(NULL),
    exportLocalIds(NULL),
    importProcs(NULL),
    importToPart(NULL),
    exportProcs(NULL),
    exportToPart(NULL),
    use_spatial_sort(true),
    nevery(1000),
    last_partitioning(-1)
{
  // Initialize the Zoltan library with a C language call
  float version;
  Zoltan_Initialize(0, NULL, &version);

  // Dynamically create Zoltan object.
  zinstance = new Zoltan(MPI::COMM_SELF);
  part_lists = new std::vector<int>[comm->nthreads];

  me = comm->me;
  box_change = 0;

  char num_parts[10];
  sprintf(num_parts, "%d", comm->nthreads);

  // Configure instance
  zinstance->Set_Param("DEBUG_LEVEL", "0");
  zinstance->Set_Param("LB_METHOD", "RCB");
  zinstance->Set_Param("NUM_GID_ENTRIES", "1");
  zinstance->Set_Param("NUM_LID_ENTRIES", "1");
  zinstance->Set_Param("NUM_LOCAL_PARTS", num_parts);
  zinstance->Set_Param("RETURN_LISTS", "PARTS");
  zinstance->Set_Param("KEEP_CUTS", "1");
  zinstance->Set_Param("RCB_OUTPUT_LEVEL", "0");
  zinstance->Set_Param("RCB_RECTILINEAR_BLOCKS", "0");

  if(argc % 2) error->all(FLERR,"Bad partitioner parameters");

  for(int a = 0; a < argc; a += 2) {
    if(strcmp(argv[a], "USE_SPATIAL_SORT") == 0) {
      if(strcmp(argv[a+1], "yes") == 0) {
        use_spatial_sort = true;
      } else if(strcmp(argv[a+1], "no") == 0) {
        use_spatial_sort = false;
      } else {
        error->all(FLERR, "Bad USE_SPATIAL_SORT value");
      }
    } else if(strcmp(argv[a], "EVERY") == 0) {
      nevery = atoi(argv[a+1]);
    } else {
      zinstance->Set_Param(argv[a], argv[a+1]);
    }
  }

  // Setup callbacks
  zinstance->Set_Num_Obj_Fn(ZoltanPartitioner::getNumberOfAssignedObjects, this);
  zinstance->Set_Obj_List_Fn(ZoltanPartitioner::getObjectList, this);
  zinstance->Set_Num_Geom_Fn(ZoltanPartitioner::getObjectSize, this);
  zinstance->Set_Geom_Multi_Fn(ZoltanPartitioner::getObjects,  this);
}

ZoltanPartitioner::~ZoltanPartitioner()
{
  delete [] part_lists;
  delete zinstance;
}

// number of items on processor
int ZoltanPartitioner::getNumberOfAssignedObjects(void *data, int *ierr)
{
  ZoltanPartitioner * p = static_cast<ZoltanPartitioner*>(data);
  *ierr = 0;
  return p->atom->nlocal;
}

void ZoltanPartitioner::getObjectList(void *data, int, int,
    ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim,
    float *obj_wgts, int *ierr){

  ZoltanPartitioner * p = static_cast<ZoltanPartitioner*>(data);
  int nlocal = p->atom->nlocal;
  int * tag = p->atom->tag;

  for(int i = 0; i < nlocal; i++) {
    global_ids[i] = tag[i];
    local_ids[i]  = i;
  }

  if(wgt_dim) {
    PairGran * pg = static_cast<PairGran*>(p->force->pair);

    if(neighbor_list_exists) {
      int * numneigh = pg->list->numneigh;

      // this will produce garbage if box has changed
      int max_neigh  = *max_element(&numneigh[0], &numneigh[nlocal-1]);
      //printf("max neigh: %d\n", max_neigh);

      for(int i = 0; i < nlocal; i++) {
        if(box_change || max_neigh == 0) {
          // neighbor list data outdated, use this for the time being
          obj_wgts[i] = 1.0f;
        } else {
          obj_wgts[i] = ((float)numneigh[i]) / max_neigh;
        }
        if(obj_wgts[i] < 0.0) printf("[%d] WARNING at %d: %f = %d / %d\n", me, i, obj_wgts[i], numneigh[i], max_neigh);
      }
    } else {
      std::fill_n(obj_wgts, nlocal, 1.0f);
    }
  }

  *ierr = 0;
}


int ZoltanPartitioner::getObjectSize(void *, int *ierr) {
  *ierr = 0;
  return 3;
}

void ZoltanPartitioner::getObjects(void *data, int, int, int num_obj, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR local_ids, int, double *geom_vec, int *ierr) {
  int next = 0;

  ZoltanPartitioner * p = static_cast<ZoltanPartitioner*>(data);
  double ** x = p->atom->x;

  for(int i = 0; i < num_obj; i++) {
    int id = local_ids[i];
    geom_vec[next++] = x[id][0];
    geom_vec[next++] = x[id][1];
    geom_vec[next++] = x[id][2];
  }

  *ierr = 0;
}

bool ZoltanPartitioner::is_cost_effective() const
{
  return atom->natoms > 1000;
}


/**
 * @brief Generate partitions periodically using Zoltan and amend them if necessary
 *
 * @param permute - permutation vector of particle indices
 * @param thread_offsets - list of thread region boundaries
 *
 * @return Partitioner::NEW_PARTITIONS if partitioning was successful and particles should be permuted,
 *         Partitioner::NO_CHANGE if existing partitions are sufficient,
 *         Partitioner::FAILED otherwise
 */
Partitioner::Result ZoltanPartitioner::generate_partitions(int * permute, std::vector<int> & thread_offsets)
{
  PairGran* pg = (PairGran*)force->pair_match("gran",1);
  if(!pg)   pg = (PairGran*)force->pair_match("gran/omp",1);
  if(!pg) {
    if(comm->me == 0) printf("ZoltanPartitioner: Granular Pair Style required!\n");
    return FAILED;
  }

  const int nlocal = atom->nlocal;
  const int nthreads = comm->nthreads;

  if(nlocal == 0 || nthreads == 1) return FAILED;

  // update only nevery times
  bool no_update_needed = (update->ntimestep - last_partitioning) < nevery;

  // force full partitioning if nlocal change is larger than first partition
  bool small_change = atom->thread_offsets.size() > 0 && (std::abs(atom->thread_offsets.back() - nlocal) < atom->thread_offsets[1]);

  if(last_partitioning > 0 && no_update_needed && small_change){
    return amend_partitions(permute, thread_offsets);
  }

  return update_partitions(permute, thread_offsets);
}

/**
 * @brief Amend existing partitions with new particles and move relocated particles
 *        back into a valid thread range.
 *
 * @param permute - permutation vector of particle indices
 * @param thread_offsets - list of thread region boundaries
 *
 * @return Partitioner::NEW_PARTITIONS if amending was successful and particles should be permuted,
 *         Partitioner::FAILED otherwise
 */
Partitioner::Result ZoltanPartitioner::amend_partitions(int * permute, std::vector<int> & thread_offsets) {
  const int nlocal = atom->nlocal;
  const int nthreads = comm->nthreads;
  int * thread = atom->thread;

  // partitioning happened before, avoid full repartitioning and only assign
  // new particles to existing partitions and do sort using permute
  //
  // NOTE: be aware that because of exchange, particles might have been deleted.
  // this will cause particles of the last partitions to move forward to fill
  // in the blanks.
#ifdef ZOLTAN_DEBUG
  if(comm->me == 0) printf("[%lu] Updating partitioning!\n", update->ntimestep);
#endif
  double ** x = atom->x;

  for(int tid = 0; tid < nthreads; ++tid) {
    part_lists[tid].clear();
  }

  for(int i = 0; i < nlocal; ++i) {
    // new particles or particles exchanged will be set to thread[i] = -1
    int part = thread[i];

    if(part < 0) {
      // figure out into which partition a new particle should go
      int proc = 0;
      if(zinstance->LB_Point_PP_Assign(x[i], proc, part) == ZOLTAN_OK) {
        thread[i] = part;
      } else {
        printf("ERROR in ZoltanPartitioner: Could not assign particle #%d to a partition!\n", i);
        thread[i] = -1;
        return FAILED;
      }
    }
    part_lists[part].push_back(i);
  }

  // update thread boundaries
  thread_offsets.clear();
  int n = 0;
  for(int tid = 0; tid < nthreads; ++tid) {
    std::vector<int> & part_list = part_lists[tid];

    // fill permute vector
    std::copy(part_list.begin(), part_list.end(), &permute[n]);
    thread_offsets.push_back(n);
    n += part_list.size();
  }
  thread_offsets.push_back(n);
  return NEW_PARTITIONS;
}

/**
 * @brief Partition particles using Zoltan algorithms
 *
 * @param permute - permutation vector of particle indices
 * @param thread_offsets - list of thread region boundaries
 *
 * @return Partitioner::NEW_PARTITIONS if partitioning was successful and particles should be permuted,
 *         Partitioner::NO_CHANGE if existing partitions are sufficient,
 *         Partitioner::FAILED otherwise
 */
Partitioner::Result ZoltanPartitioner::update_partitions(int * permute, std::vector<int> & thread_offsets) {
  const int nthreads = comm->nthreads;
  int * thread = atom->thread;

  box_change = domain->box_change;
  neighbor_list_exists = neighbor->ncalls > 0;

  // full partitioning using Zoltan
  int changes;
  int result = zinstance->LB_Partition(changes, numGidEntries, numLidEntries, numImport, importGlobalIds, importLocalIds, importProcs, importToPart, numExport, exportGlobalIds, exportLocalIds, exportProcs, exportToPart);

  if(result == ZOLTAN_OK) {
      if(changes) {
#ifdef ZOLTAN_DEBUG
        if(comm->me == 0) printf("[%lu] Zoltan partitioning!\n", update->ntimestep);
#endif

        for(int tid = 0; tid < nthreads; ++tid) {
          part_lists[tid].clear();
        }

        // something changed => create permutation
        for(int i = 0; i < numExport; i++) {
          const int tid = exportToPart[i];
          const int pidx = exportLocalIds[i];
          part_lists[tid].push_back(pidx);
          thread[pidx] = tid;
        }

        // update thread boundaries
        thread_offsets.clear();
        int n = 0;
        for(int tid = 0; tid < nthreads; ++tid) {
          std::vector<int> & part_list = part_lists[tid];
          thread_offsets.push_back(n);
          n += part_list.size();
        }
        thread_offsets.push_back(n);

        // create permute vector for partitioning (optional: use spatial sorting)
        if(use_spatial_sort) {
          // TODO: pre-sorting really necessary?
          #pragma omp parallel
          {
            const int tid = omp_get_thread_num();
            std::vector<int> & part_list = part_lists[tid];
            std::sort(part_list.begin(), part_list.end());
          }

          // create permute in serial
          // duplicating memory buffers does not scale
          for(int tid = 0; tid < nthreads; ++tid) {
            std::vector<int> & part_list = part_lists[tid];
            const int offset = thread_offsets[tid];
            atom->fill_permute_by_spatial_sorted_bins(part_list, &permute[offset]);
          }
        } else {
          for(int tid = 0; tid < nthreads; ++tid) {
            std::vector<int> & part_list = part_lists[tid];
            const int offset = thread_offsets[tid];
            std::copy(part_list.begin(), part_list.end(), &permute[offset]);
          }
        }

        last_partitioning = update->ntimestep;

        //if(comm->me == 0) printf("offset: %d\n", n);
        zinstance->LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs, &importToPart);
        zinstance->LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs, &exportToPart);
        return NEW_PARTITIONS;
      }
      else {
        // no changes
        return NO_CHANGE;
      }
  } else {
    // something went wrong
    return FAILED;
  }

  // TODO: should this really happen here, or can they be reused?
  zinstance->LB_Free_Part(&importGlobalIds, &importLocalIds, &importProcs, &importToPart);
  zinstance->LB_Free_Part(&exportGlobalIds, &exportLocalIds, &exportProcs, &exportToPart);
  return FAILED;
}
