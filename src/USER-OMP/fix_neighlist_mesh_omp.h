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

#ifdef FIX_CLASS

FixStyle(neighlist/mesh/omp,FixNeighlistMeshOMP)

#else

#ifndef LMP_FIX_NEIGHLIST_MESH_OMP_H
#define LMP_FIX_NEIGHLIST_MESH_OMP_H

#include "fix.h"
#include "container.h"
#include <vector>
#include <iterator>
//#include <map>
#include "fix_neighlist_mesh.h"

namespace LAMMPS_NS
{


class FixNeighlistMeshOMP : public FixNeighlistMesh
{
  public:
    FixNeighlistMeshOMP(LAMMPS *lmp, int narg, char **arg);
    virtual ~FixNeighlistMeshOMP();

    virtual void initializeNeighlist();
    virtual void pre_force(int vflag); //NP builds tri neigh list if flag

    std::vector<int> partition_global_indices;
    std::vector<int> partition_global_triangles;
    std::vector<int> partition_global_thread;

    std::vector<int> partition_local_indices;
    std::vector<int> partition_local_triangles;
    std::vector<int> partition_local_thread;
    //std::map<int, int> thread_mapping;

    std::vector<int> thread_offsets;
    int * partition_begin(int tid);
    int * partition_end(int tid);
    bool in_thread_partition(int tid, int i);

  protected:
    void handleTriangle(int iTri);

    int quicksort_partition(int p, int r);
    void quicksort(int p, int r);
    void sort_triangles_by_nchecked();

    std::vector<int> trilist;
    int sortfreq;
    bigint nextsort;
};

} /* namespace LAMMPS_NS */


#endif /* FIX_MESH_NEIGHLIST_OMP_H_ */
#endif /* FIX_CLASS */
