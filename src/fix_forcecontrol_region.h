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

#if defined(LAMMPS_VTK)

#ifdef FIX_CLASS

FixStyle(forcecontrol/region,FixForceControlRegion)

#else

#ifndef LMP_FIX_FORCECONTROL_REGION_H
#define LMP_FIX_FORCECONTROL_REGION_H

#include "fix.h"
#include <map>
#include <set>
#include <vector>

namespace LAMMPS_NS {

class FixForceControlRegion : public Fix {
 public:
  FixForceControlRegion(class LAMMPS *, int, char **);
  ~FixForceControlRegion();
  virtual void post_create();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void min_setup(int);
  virtual void post_force(int);
  virtual void reset_dt();
  virtual double compute_scalar();
  virtual double compute_vector(int);
  virtual int modify_param(int narg, char **arg);
  virtual void end_of_step();

  virtual void write_restart(FILE *fp);
  virtual void restart(char *buf);

 protected:
  double cg_actual_, cg3_actual_, cg_target_, cg3_target_, cg_ratio_, cg3_ratio_;
  double cg_max_; // input only if single-partition, else equal to actual_cg
  double const_part_, sinesq_part_, used_part_;

 private:
  double *xvalue;
  double *yvalue;
  double *zvalue;
  double *const_part_cell_; // part of cell where control force is applied constantly
  double *used_part_cell_; // part of cell where control force is applied
  double *sinesq_part_cell_; // part of cell where correction force is faded out
  int *ncontrolled_cell_; // # of controlled particles in each cell
  double foriginal[4],foriginal_all[4];
  int force_flag;

  double axis_[3];
  double vel_max_[3];
  double vel_min_[3];
  double ctrl_op_[3];
  double sp_vec_[3]; // set point
  double pv_vec_[3]; // process value
  double err_[3];
  double kp_,ki_,kd_;

  int ctrl_style_; // stress or velocity control style
  // timesteps
  double dtf_,dtv_,dtv_inverse_;
  double fadex_,fadey_,fadez_;

  class FixAveEulerRegion *actual_;
  int ncells_max_;
  double **old_pv_vec_;
  double **sum_err_;

  std::set<int> active_; // set of active cells
  std::vector<bool> modifier_; // if true apply massflow corrections
  std::map<class FixScaleDiameter*, std::set<int> > modifier_scale_;
  double acceptable_deviation_min;
  double acceptable_deviation_max;
  char limit_velocity_;
  double limit[3];
  double force_limit_percentage_threshold_;

 private:
  class FixAveEulerRegion *target_;

  virtual int target_couple_every();
  virtual int target_has_cell_id(int cell_id);
  virtual int target_cell_index(int cell_id);
  virtual int target_cell_count(int cell_index);
  virtual double target_cell_mass(int cell_index);
  virtual double target_cell_vol_fr(int cell_index);
  virtual double target_cell_stress_xx(int cell_index);
  virtual double target_cell_stress_yy(int cell_index);
  virtual double target_cell_stress_zz(int cell_index);
  virtual double target_cell_v_ave_x(int cell_index);
  virtual double target_cell_v_ave_y(int cell_index);
  virtual double target_cell_v_ave_z(int cell_index);
  virtual double target_cell_v_min_x(int cell_index);
  virtual double target_cell_v_min_y(int cell_index);
  virtual double target_cell_v_min_z(int cell_index);
  virtual double target_cell_v_max_x(int cell_index);
  virtual double target_cell_v_max_y(int cell_index);
  virtual double target_cell_v_max_z(int cell_index);

  virtual void post_create_stress_part();

  virtual void receive_post_create_data() {}
  virtual void receive_coupling_data() {}

  void grow_cell_arrays(int ncells);
};

}

#endif
#endif
#endif
