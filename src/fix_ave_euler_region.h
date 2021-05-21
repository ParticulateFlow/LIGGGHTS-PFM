/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
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

#if defined(LAMMPS_VTK)

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
  virtual ~FixAveEulerRegion();

  virtual void post_create();
  virtual int setmask();
  virtual void init();

  virtual double compute_array(int i, int j);

  virtual double compute_array_by_id(int cell_id, int j);

  virtual double cell_volume(int i) const;

  virtual double cell_center(int i, int j) const;

  virtual bool has_cell_id(int cell_id) const
  { return cellid2index_.find(cell_id) != cellid2index_.end(); }

  virtual int cell_id(int i) const
  { return cellid_[i]; }

  virtual int cell(int cell_id)
  { return cellid2index_[cell_id]; }

  virtual double cell_v_min(int i, int j) const
  { return v_min_[i][j]; }

  virtual double cell_v_max(int i, int j) const
  { return v_max_[i][j]; }

  virtual void cell_bounds(int i, double bounds[6]) const;

  virtual void cell_points(int i, double points[24]) const;

  virtual double* cell_vector_property(int i, const char* property);

 private:

  virtual void setup_bins();
  virtual void bin_atoms();
  virtual void lazy_bin_atoms(int i);
  virtual void calculate_eu();
  virtual bool is_inside_bounds(double bounds[6], double *pos) const;

  char *idregion_grid_;
  class Region *region_grid_;
  class RegHexMesh *region_grid_mesh_hex_;

 protected:
  // cell-based min/max velocity
  double **v_min_;
  double **v_max_;

  std::map<int, int> cellid2index_;
  std::vector<int> cellid_;

 private:
  virtual void send_post_create_data() {}
  virtual void send_coupling_data() {}
};

}

#endif
#endif
#endif
