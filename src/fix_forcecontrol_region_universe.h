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

#ifdef FIX_CLASS

FixStyle(forcecontrol/region/universe,FixForceControlRegionUniverse)

#else

#ifndef LMP_FIX_FORCECONTROL_REGION_UNIVERSE_H
#define LMP_FIX_FORCECONTROL_REGION_UNIVERSE_H

#include "fix_forcecontrol_region.h"
#include <map>
#include <set>
#include <vector>

namespace LAMMPS_NS {

class FixForceControlRegionUniverse : public FixForceControlRegion {
 public:
  FixForceControlRegionUniverse(class LAMMPS *, int, char **);
  ~FixForceControlRegionUniverse();

 private:
  int couple_every_;
  int receive_from_world_;
  unsigned int idtarget_hash;
  std::map<int, int> target_cellid2index_;
  std::vector<int> target_count_;
  std::vector<double> target_mass_;
  std::vector<double> target_vol_fr_;
  double **target_stress_;
  double **target_v_ave_;
  double **target_v_min_;
  double **target_v_max_;

  virtual int target_couple_every() { return couple_every_; }

  virtual int target_has_cell_id(int cell_id)
  { return (target_cellid2index_.find(cell_id) != target_cellid2index_.end()); }

  virtual int target_cell_index(int cell_id) { return target_cellid2index_[cell_id]; }
  virtual int target_cell_count(int cell_index) { return target_count_[cell_index]; }
  virtual double target_cell_mass(int cell_index) { return target_mass_[cell_index]; }
  virtual double target_cell_vol_fr(int cell_index) { return target_vol_fr_[cell_index]; }
  virtual double target_cell_stress_xx(int cell_index) { return target_stress_[cell_index][1]; }
  virtual double target_cell_stress_yy(int cell_index) { return target_stress_[cell_index][2]; }
  virtual double target_cell_stress_zz(int cell_index) { return target_stress_[cell_index][3]; }
  virtual double target_cell_v_ave_x(int cell_index) { return target_v_ave_[cell_index][0]; }
  virtual double target_cell_v_ave_y(int cell_index) { return target_v_ave_[cell_index][1]; }
  virtual double target_cell_v_ave_z(int cell_index) { return target_v_ave_[cell_index][2]; }
  virtual double target_cell_v_min_x(int cell_index) { return target_v_min_[cell_index][0]; }
  virtual double target_cell_v_min_y(int cell_index) { return target_v_min_[cell_index][1]; }
  virtual double target_cell_v_min_z(int cell_index) { return target_v_min_[cell_index][2]; }
  virtual double target_cell_v_max_x(int cell_index) { return target_v_max_[cell_index][0]; }
  virtual double target_cell_v_max_y(int cell_index) { return target_v_max_[cell_index][1]; }
  virtual double target_cell_v_max_z(int cell_index) { return target_v_max_[cell_index][2]; }

  virtual void post_create_stress_part();

  virtual void receive_post_create_data();
  virtual void receive_coupling_data();
};

}

#endif
#endif
