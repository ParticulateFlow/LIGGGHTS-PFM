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

#ifdef PARTITIONER_CLASS

PartitionerStyle(zoltan,ZoltanPartitioner)

#else

#ifndef PARTITIONER_ZOLTAN_H_
#define PARTITIONER_ZOLTAN_H_

#include "atom.h"
#include <vector>
#include <zoltan_cpp.h>

namespace LAMMPS_NS {

class ZoltanPartitioner: public Partitioner {
  Zoltan * zinstance;
  std::vector<int> * part_lists;

  static int me;
  static int box_change;
  static bool neighbor_list_exists;

  int numGidEntries, numLidEntries;
  int numImport, numExport;
  ZOLTAN_ID_PTR importGlobalIds;
  ZOLTAN_ID_PTR importLocalIds;
  ZOLTAN_ID_PTR exportGlobalIds;
  ZOLTAN_ID_PTR exportLocalIds;
  int * importProcs;
  int * importToPart;
  int * exportProcs;
  int * exportToPart;

  bool use_spatial_sort;

  int nevery;             // frequency of full partitioning
  bigint last_partitioning;  // last timestep which full partitioning occured

  // callbacks for zoltan
  static int  getNumberOfAssignedObjects(void *data, int *ierr);
  static void getObjectList(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr);
  static int  getObjectSize(void *data, int *ierr);
  static void getObjects(void *data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr);

public:
  ZoltanPartitioner(class LAMMPS * lmp, int argc, const char * const * argv);
  virtual ~ZoltanPartitioner();

  virtual bool is_cost_effective() const;
  virtual Result generate_partitions(int * permute, std::vector<int> & thread_offsets);

private:
  Result amend_partitions(int * permute, std::vector<int> & thread_offsets);

  Result update_partitions(int * permute, std::vector<int> & thread_offsets);
};

} /* namespace LAMMPS_NS */
#endif /* PARTITIONER_ZOLTAN_H_ */

#endif
