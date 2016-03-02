/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2016- JKU Linz

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

#ifdef FIX_CLASS

FixStyle(ave/euler/region,FixAveEulerRegion)

#else

#ifndef LMP_FIX_AVE_EULER_REGION_H
#define LMP_FIX_AVE_EULER_REGION_H

#include "fix_ave_euler.h"
#include <map>
#include <vector>

namespace LAMMPS_NS {

class FixAveEulerRegion : public FixAveEuler {

 public:

  FixAveEulerRegion(class LAMMPS *, int, char **);
  ~FixAveEulerRegion();

  void post_create();
  int setmask();
  void init();

  double compute_array(int i, int j);

  double compute_array_by_id(int cell_id, int j);

  double cell_volume(int i);

  double cell_center(int i, int j);

  inline int cell_id(int i)
  { return cellid_[i]; }

  inline int cell(int cell_id)
  { return cellid2index_[cell_id]; }

  inline double cell_v_min(int i, int j)
  { return v_min_[i][j]; }

  inline double cell_v_max(int i, int j)
  { return v_max_[i][j]; }

  void cell_bounds(int i, double bounds[6]);

 private:

  void setup_bins();
  void bin_atoms();
  void calculate_eu();

  // cell-based min/max velocity
  double **v_min_;
  double **v_max_;

  char *idregion_grid_;
  class Region *region_grid_;
  class RegHexMesh *region_grid_mesh_hex_;

  std::map<int, int> cellid2index_;
  std::vector<int> cellid_;
};

}

#endif
#endif
