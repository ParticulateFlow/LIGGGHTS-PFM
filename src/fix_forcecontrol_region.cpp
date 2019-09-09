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

#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#include "fix_ave_euler_region.h"
#include "fix_forcecontrol_region.h"
#include "fix_scale_diameter.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,STRESS,VELOCITY};

extern MPI_Op MPI_ABSMIN_OP;
extern MPI_Op MPI_ABSMAX_OP;

/* ---------------------------------------------------------------------- */

FixForceControlRegion::FixForceControlRegion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  cg_target_(1.0),
  cg3_target_(1.0),
  cg_ratio_(1.0),
  const_part_(1.0),
  sinesq_part_(0.0),
  used_part_(1.0),
  xvalue(NULL),
  yvalue(NULL),
  zvalue(NULL),
  const_part_cell_(NULL),
  used_part_cell_(NULL),
  sinesq_part_cell_(NULL),
  ncontrolled_cell_(NULL),
  kp_(1.),
  ki_(0.),
  kd_(0.),
  ctrl_style_(NONE),
  dtf_(1.0),
  dtv_(1.0),
  dtv_inverse_(1.0),
  fadex_ (1.0),
  fadey_ (1.0),
  fadez_ (1.0),
  ncells_max_(0),
  old_pv_vec_(NULL),
  sum_err_(NULL),
  acceptable_deviation_min(0.97),
  acceptable_deviation_max(1.03),
  limit_velocity_(false)
{
  if (narg < 6) error->all(FLERR,"Illegal fix forcecontrol/region command");

  restart_global = 1;
  virial_flag = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  vectorZeroize3D(axis_);
  vectorZeroize3D(ctrl_op_);

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ctrlPV") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'ctrlPV'");
      if (strcmp(arg[iarg+1],"stress") == 0) {
        ctrl_style_ = STRESS;
      } else if (strcmp(arg[iarg+1],"velocity") == 0) {
        ctrl_style_ = VELOCITY;
        limit_velocity_ = true;
      }
      else error->fix_error(FLERR,this,"only 'stress' and 'velocity' are valid arguments for ctrlPV");
      iarg = iarg + 2;
    } else if (strcmp(arg[iarg],"actual_val") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'actual_val'");
      ++iarg;

      int ifix = modify->find_fix(arg[iarg]);

      if(ifix < 0)
          error->all(FLERR,"Illegal fix forcecontrol/region command, invalid ID for fix ave/euler/region provided");

      if(strncmp(modify->fix[ifix]->style,"ave/euler/region",16))
          error->all(FLERR,"Illegal fix forcecontrol/region command, fix is not of type ave/euler/region");

      actual_ = static_cast<FixAveEulerRegion*>(modify->fix[ifix]);

      ++iarg;
    } else if(strcmp(arg[iarg],"target_val") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'target_val'");
      ++iarg;
      if (strcmp(style,"forcecontrol/region") == 0) {
        int ifix = modify->find_fix(arg[iarg]);

        if(ifix < 0)
            error->all(FLERR,"Illegal fix forcecontrol/region command, invalid ID for fix ave/euler/region provided");

        if(strncmp(modify->fix[ifix]->style,"ave/euler/region",16))
            error->all(FLERR,"Illegal fix forcecontrol/region command, fix is not of type ave/euler/region");

        target_ = static_cast<FixAveEulerRegion*>(modify->fix[ifix]);
      }

      ++iarg;
    } else if(strcmp(arg[iarg],"kp") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments");
      kp_ = force->numeric(FLERR,arg[iarg+1]);
      iarg = iarg+2;
    } else if(strcmp(arg[iarg],"ki") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments");
      ki_ = force->numeric(FLERR,arg[iarg+1]);
      iarg = iarg+2;
    } else if(strcmp(arg[iarg],"kd") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments");
      kd_ = force->numeric(FLERR,arg[iarg+1]);
      iarg = iarg+2;
    } else if (strcmp(arg[iarg],"velocity_limit") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      if (strcmp(arg[iarg+1],"on") == 0)
        limit_velocity_ = true;
      iarg += 2;
    } else if (strcmp(arg[iarg],"cg") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      cg_target_ = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(style,"forcecontrol/region") == 0) {
      error->all(FLERR,"Illegal fix forcecontrol/region command");
    } else {
      ++iarg;
    }
  }

  if(!actual_)
    error->fix_error(FLERR, this, "missing ave/euler/region for actual_val");
  if((strcmp(style,"forcecontrol/region") == 0) && !target_)
    error->fix_error(FLERR, this, "missing ave/euler/region for target_val");

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixForceControlRegion::~FixForceControlRegion()
{
  memory->destroy(xvalue);
  memory->destroy(yvalue);
  memory->destroy(zvalue);
  memory->destroy(const_part_cell_);
  memory->destroy(used_part_cell_);
  memory->destroy(sinesq_part_cell_);
  memory->destroy(ncontrolled_cell_);
  memory->destroy(old_pv_vec_);
  memory->destroy(sum_err_);
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::post_create()
{
  receive_post_create_data();
  cg3_target_ = cg_target_*cg_target_*cg_target_;
  cg_ratio_ = force->cg()/cg_target_;

  // all cells active by default
  int ncells = actual_->ncells();
  for (int icell=0; icell<ncells; ++icell) {
    active_.insert(icell);
  }
  modifier_.resize(ncells, false);

  if (ctrl_style_ == STRESS) {
    post_create_stress_part();
  }
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::post_create_stress_part()
{
  double maxrd,minrd;
  modify->max_min_rad(maxrd,minrd);
  const_part_ = (2.*maxrd/cg_target_)*1.2;
  used_part_ = (2.*maxrd/cg_target_)*1.4;
  sinesq_part_ = used_part_ - const_part_;
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::init()
{
  // get timestep
  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::min_setup(int vflag)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::post_force(int vflag)
{
  if (vflag) v_setup(vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *rmass = atom->rmass;
  int ncells = actual_->ncells();
  double fx, fy, fz;
  double vx, vy, vz;

  // reallocate old_pv_vec_ array if necessary
  if (ncells > ncells_max_) {
    memory->grow(xvalue,ncells,"forcecontrol/region:xvalue");
    memory->grow(yvalue,ncells,"forcecontrol/region:yvalue");
    memory->grow(zvalue,ncells,"forcecontrol/region:zvalue");
    memory->grow(const_part_cell_,ncells,"forcecontrol/region:const_part_cell_");
    memory->grow(used_part_cell_,ncells,"forcecontrol/region:used_part_cell_");
    memory->grow(sinesq_part_cell_,ncells,"forcecontrol/region:sinesq_part_cell_");
    memory->grow(ncontrolled_cell_,ncells,"forcecontrol/region:ncontrolled_cell_");
    memory->grow(old_pv_vec_,ncells,3,"forcecontrol/region:old_pv_mag_");
    memory->grow(sum_err_,ncells,3,"forcecontrol/region:sum_err_");

    for (; ncells_max_ < ncells; ++ncells_max_) {
      xvalue[ncells_max_] = 0.0;
      yvalue[ncells_max_] = 0.0;
      zvalue[ncells_max_] = 0.0;
      const_part_cell_[ncells_max_] = const_part_;
      used_part_cell_[ncells_max_] = used_part_;
      sinesq_part_cell_[ncells_max_] = sinesq_part_;
      ncontrolled_cell_[ncells_max_] = -1;

      for (int i = 0; i < 3; ++i) {
          old_pv_vec_[ncells_max_][i] = 0.0;
          sum_err_[ncells_max_][i] = 0.0;
      }
    }
  }

  std::set<int>::iterator it_cell = active_.begin();
  for (; it_cell!=active_.end(); ++it_cell) { // active actual_ cell indices

    int cell_id = actual_->cell_id(*it_cell);

    int tcell = target_cell_index(cell_id);

    // skip cell if target_ cell has less than half of equivalent cg particle
    // mostly this happens when a stream of particles enters an empty cell
    if (target_cell_vol_fr(tcell) <= actual_->cell_vol_fr(*it_cell)*0.5) {
      continue;
    }

    // clear controller history for empty cells
    if (actual_->cell_count(*it_cell) < 1 || target_cell_count(tcell) < 1) {
      old_pv_vec_[*it_cell][0] = 0.;
      old_pv_vec_[*it_cell][1] = 0.;
      old_pv_vec_[*it_cell][2] = 0.;

      sum_err_[*it_cell][0] = 0.;
      sum_err_[*it_cell][1] = 0.;
      sum_err_[*it_cell][2] = 0.;

      continue;
    }

    switch (ctrl_style_) {
    case STRESS:
      {
        const double vol_fr_mismatch_correction = actual_->cell_vol_fr(*it_cell)/target_cell_vol_fr(tcell);
        actual_->set_cell_weight(*it_cell, vol_fr_mismatch_correction);
        // calculate half average normal force per particle (half, because same force acting from 2 sides)
        const double r_ave = actual_->cell_radius(*it_cell); // average fg radius (assume that in cg we get the same r_ave when scaled with 1/cg)
        const double vol_cell = actual_->cell_volume(*it_cell); // equal for fg and cg
        const double pre_t = 0.5 * vol_cell / r_ave; // (along positive and negative axis -> * 2)
        const double pre_a = pre_t * vol_fr_mismatch_correction; // (along positive and negative axis -> * 2)

        pv_vec_[0] = -(pre_a*actual_->cell_stress(*it_cell, 0)/actual_->cell_count(*it_cell)); // fg
        pv_vec_[1] = -(pre_a*actual_->cell_stress(*it_cell, 1)/actual_->cell_count(*it_cell)); // fg
        pv_vec_[2] = -(pre_a*actual_->cell_stress(*it_cell, 2)/actual_->cell_count(*it_cell)); // fg
        sp_vec_[0] = -(pre_t*target_cell_stress_xx(tcell)/(target_cell_count(tcell)*cg3_target_));// cg
        sp_vec_[1] = -(pre_t*target_cell_stress_yy(tcell)/(target_cell_count(tcell)*cg3_target_));// cg
        sp_vec_[2] = -(pre_t*target_cell_stress_zz(tcell)/(target_cell_count(tcell)*cg3_target_));// cg

        // do nothing in case of negative pressure
        if(pv_vec_[0] < 0.) pv_vec_[0] = 0.;
        if(pv_vec_[1] < 0.) pv_vec_[1] = 0.;
        if(pv_vec_[2] < 0.) pv_vec_[2] = 0.;

        // also clear error sum if there is no target pressure
        if(sp_vec_[0] <= 0.) { sp_vec_[0] = 0.; sum_err_[*it_cell][0] = 0.; }
        if(sp_vec_[1] <= 0.) { sp_vec_[1] = 0.; sum_err_[*it_cell][1] = 0.; }
        if(sp_vec_[2] <= 0.) { sp_vec_[2] = 0.; sum_err_[*it_cell][2] = 0.; }

        // get stress directions from region mesh/hex (vector properties stored in vtk file)
        double *stress_ctrl_dir = actual_->cell_vector_property(*it_cell, "stress_ctrl_dir");
        axis_[0] = stress_ctrl_dir[0];
        axis_[1] = stress_ctrl_dir[1];
        axis_[2] = stress_ctrl_dir[2];
        break;
      }
    case VELOCITY:
      {
        pv_vec_[0] = actual_->cell_v_av(*it_cell, 0);
        pv_vec_[1] = actual_->cell_v_av(*it_cell, 1);
        pv_vec_[2] = actual_->cell_v_av(*it_cell, 2);

        double massflow_correction = 1.;

        if (modifier_[tcell] && actual_->cell_vol_fr(*it_cell) > 0.)
          massflow_correction = target_cell_vol_fr(tcell)/actual_->cell_vol_fr(*it_cell);

        sp_vec_[0] = massflow_correction*target_cell_v_ave_x(tcell);
        sp_vec_[1] = massflow_correction*target_cell_v_ave_y(tcell);
        sp_vec_[2] = massflow_correction*target_cell_v_ave_z(tcell);

        axis_[0] = axis_[1] = axis_[2] = 1.;

        break;
      }
    default:
      error->fix_error(FLERR, this, "unknown ctrl style");
      break;
    }

    // simple PID-controller

    // calc error and sum of the errors
    err_[0] = sp_vec_[0] - pv_vec_[0];
    err_[1] = sp_vec_[1] - pv_vec_[1];
    err_[2] = sp_vec_[2] - pv_vec_[2];

    if (ctrl_style_ == STRESS) {
      // axis value (one of {-1,0,1}) turns on/off control
      err_[0] *= fabs(axis_[0]);
      err_[1] *= fabs(axis_[1]);
      err_[2] *= fabs(axis_[2]);
    }

    sum_err_[*it_cell][0] += err_[0] * dtv_;
    sum_err_[*it_cell][1] += err_[1] * dtv_;
    sum_err_[*it_cell][2] += err_[2] * dtv_;

    // derivative term
    // force used instead of the error in order to avoid signal spikes in case of change of the set point
    double dfdt[3];
    dfdt[0] = -(pv_vec_[0] - old_pv_vec_[*it_cell][0])/dtv_;
    dfdt[1] = -(pv_vec_[1] - old_pv_vec_[*it_cell][1])/dtv_;
    dfdt[2] = -(pv_vec_[2] - old_pv_vec_[*it_cell][2])/dtv_;

    // vel points opposite to force vector
    ctrl_op_[0] = err_[0] * kp_ + sum_err_[*it_cell][0] * ki_ + dfdt[0] * kd_;
    ctrl_op_[1] = err_[1] * kp_ + sum_err_[*it_cell][1] * ki_ + dfdt[1] * kd_;
    ctrl_op_[2] = err_[2] * kp_ + sum_err_[*it_cell][2] * ki_ + dfdt[2] * kd_;

    // save process value for next timestep
    old_pv_vec_[*it_cell][0] = pv_vec_[0];
    old_pv_vec_[*it_cell][1] = pv_vec_[1];
    old_pv_vec_[*it_cell][2] = pv_vec_[2];

    if (ctrl_style_ == STRESS) {
      // do not apply an outward force; particles should relax in that direction automatically
      xvalue[*it_cell] = (ctrl_op_[0] > 0.) ? (axis_[0] * ctrl_op_[0]) : 0.;
      yvalue[*it_cell] = (ctrl_op_[1] > 0.) ? (axis_[1] * ctrl_op_[1]) : 0.;
      zvalue[*it_cell] = (ctrl_op_[2] > 0.) ? (axis_[2] * ctrl_op_[2]) : 0.;
    } else {
      xvalue[*it_cell] = ctrl_op_[0];
      yvalue[*it_cell] = ctrl_op_[1];
      zvalue[*it_cell] = ctrl_op_[2];
    }

    // TODO: remove
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
    force_flag = 0;
    ncontrolled_cell_[*it_cell] = 0;

    if (xvalue[*it_cell]!=0. || yvalue[*it_cell]!=0. || zvalue[*it_cell]!=0.) {
      double bounds[6];
      actual_->cell_bounds(*it_cell, bounds);

      vel_min_[0] = target_cell_v_min_x(tcell);
      vel_min_[1] = target_cell_v_min_y(tcell);
      vel_min_[2] = target_cell_v_min_z(tcell);
      vel_max_[0] = target_cell_v_max_x(tcell);
      vel_max_[1] = target_cell_v_max_y(tcell);
      vel_max_[2] = target_cell_v_max_z(tcell);

      // allow small deviation from min/max velocity
      vel_min_[0] *= (vel_min_[0] > 0.) ? acceptable_deviation_min : acceptable_deviation_max;
      vel_min_[1] *= (vel_min_[1] > 0.) ? acceptable_deviation_min : acceptable_deviation_max;
      vel_min_[2] *= (vel_min_[2] > 0.) ? acceptable_deviation_min : acceptable_deviation_max;
      vel_max_[0] *= (vel_max_[0] > 0.) ? acceptable_deviation_max : acceptable_deviation_min;
      vel_max_[1] *= (vel_max_[1] > 0.) ? acceptable_deviation_max : acceptable_deviation_min;
      vel_max_[2] *= (vel_max_[2] > 0.) ? acceptable_deviation_max : acceptable_deviation_min;

      for (int i = actual_->cell_head(*it_cell); i >= 0; i = actual_->cell_ptr(i)) {
        if (i >= nlocal)
          continue;
        if (mask[i] & groupbit) {
          fadex_ = fadey_ = fadez_ = 1.0;

          // fade out factors
          if (ctrl_style_ == VELOCITY) {

            fadex_ = fadey_ = fadez_ = rmass[i]*dtv_inverse_;

          } else if (ctrl_style_ == STRESS  && sinesq_part_cell_[*it_cell] > 0.0) {

            if (axis_[0] != 0.) {
              // leftward force and atom position to the left of const force part
              if (axis_[0] < 0. && x[i][0] < bounds[1] - const_part_cell_[*it_cell]) {
                if(x[i][0] < bounds[1] - used_part_cell_[*it_cell]) {
                  fadex_ = 0.;
                } else {
                  fadex_ = sin(M_PI*0.5*(x[i][0] - (bounds[1]-used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadex_ *= fadex_;
                }
              // rightward force and atom position to the right of const force part
              } else if (axis_[0] > 0. && x[i][0] > bounds[0] + const_part_cell_[*it_cell]) {
                if(x[i][0] > bounds[0] + used_part_cell_[*it_cell]) {
                  fadex_ = 0.;
                } else {
                  fadex_ = sin(M_PI*0.5*(x[i][0] - (bounds[0]+used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadex_ *= fadex_;
                }
              }
            }

            if (axis_[1] != 0.) {
              if (axis_[1] < 0. && x[i][1] < bounds[3] - const_part_cell_[*it_cell]) {
                if(x[i][1] < bounds[3] - used_part_cell_[*it_cell]) {
                  fadey_ = 0.;
                } else {
                  fadey_ = sin(M_PI*0.5*(x[i][1] - (bounds[3]-used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadey_ *= fadey_;
                }
              } else if (axis_[1] > 0. && x[i][1] > bounds[2] + const_part_cell_[*it_cell]) {
                if(x[i][1] > bounds[2] + used_part_cell_[*it_cell]) {
                  fadey_ = 0.;
                } else {
                  fadey_ = sin(M_PI*0.5*(x[i][1] - (bounds[2]+used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadey_ *= fadey_;
                }
              }
            }

            if (axis_[2] != 0.) {
              // downward force and atom position below const force part
              if (axis_[2] < 0. && x[i][2] < bounds[5] - const_part_cell_[*it_cell]) {
                // atom position below used part
                if(x[i][2] < bounds[5] - used_part_cell_[*it_cell]) {
                  fadez_ = 0.;
                // atom position in sine squared part
                } else {
                  fadez_ = sin(M_PI*0.5* (x[i][2] - (bounds[5]-used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadez_ *= fadez_;
                }
              // upward force and atom position above const force part
              } else if (axis_[2] > 0. && x[i][2] > bounds[4] + const_part_cell_[*it_cell]) {
                // atom position above used part
                if(x[i][2] > bounds[4] + used_part_cell_[*it_cell]) {
                  fadez_ = 0.;
                // atom position in sine squared part
                } else {
                  fadez_ = sin(M_PI*0.5* (x[i][2] - (bounds[4]+used_part_cell_[*it_cell]))/sinesq_part_cell_[*it_cell]);
                  fadez_ *= fadez_;
                }
              }
            }
          }


          if (limit_velocity_ && !modifier_[tcell]) {
            limit[0] = limit[1] = limit[2] = 1.0;
            double dtfm = dtf_ / rmass[i];

            fx = fadex_ * xvalue[*it_cell];

#define FORCE_LIMIT_PERCENTAGE_TRESHOLD 0.05

            if (fabs(fx) > fabs(f[i][0]*FORCE_LIMIT_PERCENTAGE_TRESHOLD)) {
              vx = v[i][0] + dtfm * (f[i][0] + fx);

              if (vx > vel_max_[0]) {
                if (vx - dtfm * fx < vel_max_[0])
                  limit[0] = (vel_max_[0] - v[i][0] - dtfm * f[i][0])/(dtfm * fx);
                else if (fx > 0.)
                  limit[0] = 0.;
              } else if (vx < vel_min_[0]) {
                if (vx - dtfm * fx > vel_min_[0])
                  limit[0] = (vel_min_[0] - v[i][0] - dtfm * f[i][0])/(dtfm * fx);
                else if (fx < 0)
                  limit[0] = 0.;
              }
            }

            fy = fadey_ * yvalue[*it_cell];

            if (fabs(fy) > fabs(f[i][1]*FORCE_LIMIT_PERCENTAGE_TRESHOLD)) {
              vy = v[i][1] + dtfm * (f[i][1] + fy);
              if (vy > vel_max_[1]) {
                if (vy - dtfm * fy < vel_max_[1])
                  limit[1] = (vel_max_[1] - v[i][1] - dtfm * f[i][1])/(dtfm * fy);
                else if (fy > 0.)
                  limit[1] = 0.;
              } else if (vy < vel_min_[1]) {
                if (vy - dtfm * fy > vel_min_[1])
                  limit[1] = (vel_min_[1] - v[i][1] - dtfm * f[i][1])/(dtfm * fy);
                else if (fy < 0.)
                  limit[1] = 0.;
              }
            }

            fz = fadez_ * zvalue[*it_cell];

            if (fabs(fz) > fabs(f[i][2]*FORCE_LIMIT_PERCENTAGE_TRESHOLD)) {
              vz = v[i][2] + dtfm * (f[i][2] + fz);
              if (vz > vel_max_[2]) {
                if (vz - dtfm * fz < vel_max_[2])
                  limit[2] = (vel_max_[2] - v[i][2] - dtfm * f[i][2])/(dtfm * fz);
                else if (fz > 0.)
                  limit[2] = 0.;
              } else if (vz < vel_min_[2]) {
                if (vz - dtfm * fz > vel_min_[2])
                  limit[2] = (vel_min_[2] - v[i][2] - dtfm * f[i][2])/(dtfm * fz);
                else if (fz < 0.)
                  limit[2] = 0.;
              }
            }

            double a = std::min(limit[0], std::min(limit[1], limit[2]));
            f[i][0] += a * fx;
            f[i][1] += a * fy;
            f[i][2] += a * fz;

#define DEBUG_EPS 0.0
            if (fabs(a*fx) > DEBUG_EPS || fabs(a*fy) > DEBUG_EPS || fabs(a*fz) > DEBUG_EPS ) {
              ++ncontrolled_cell_[*it_cell];
            }
            if (vflag) {
                const double delta = 2.0*actual_->cell_radius(*it_cell);
                const double pre = 0.5*a;
                vatom[i][0] += pre*fabs(fx)*delta;
                vatom[i][1] += pre*fabs(fy)*delta;
                vatom[i][2] += pre*fabs(fz)*delta;
            }
          } else {
            f[i][0] += fadex_ * xvalue[*it_cell];
            f[i][1] += fadey_ * yvalue[*it_cell];
            f[i][2] += fadez_ * zvalue[*it_cell];

            if (   fabs(fadex_ * xvalue[*it_cell]) > DEBUG_EPS
                || fabs(fadey_ * yvalue[*it_cell]) > DEBUG_EPS
                || fabs(fadez_ * zvalue[*it_cell]) > DEBUG_EPS ) {
              ++ncontrolled_cell_[*it_cell];
            }
            if (vflag) {
                const double delta = 2.0*actual_->cell_radius(*it_cell);
                vatom[i][0] += 0.5*fabs(fadex_ * xvalue[*it_cell])*delta;
                vatom[i][1] += 0.5*fabs(fadey_ * yvalue[*it_cell])*delta;
                vatom[i][2] += 0.5*fabs(fadez_ * zvalue[*it_cell])*delta;
            }
          }
        }
      }
    }
  }

  // apply scale modifier to correct volume fraction
  std::map<class FixScaleDiameter*, std::set<int> >::iterator it = modifier_scale_.begin();
  for (; it!=modifier_scale_.end(); ++it) {
    std::set<int>::iterator it_cell = it->second.begin();
    double mass_ratio = 0.;
    for (; it_cell!=it->second.end(); ++it_cell) {
      int cell_id = actual_->cell_id(*it_cell);
      int tcell = target_cell_index(cell_id);
      if(actual_->cell_vol_fr(*it_cell) > 0.)
        mass_ratio += target_cell_mass(tcell)/actual_->cell_mass(*it_cell); // fg/cg
      else
        mass_ratio += 1.0;
    }
    mass_ratio /= it->second.size();

    if (mass_ratio > 1.005) {
      it->first->set_scale(std::max(cg_ratio_, it->first->get_scale() * 0.9999));
    } else if (mass_ratio < 0.995) {
      it->first->set_scale(std::min(1.0,it->first->get_scale() * 1.0001));
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::end_of_step()
{
  receive_coupling_data();
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::modify_param(int narg, char **arg)
{
    if (narg < 2) error->fix_error(FLERR,this,"Illegal fix_modify command");
    int nusedarg = 2;

    if (strcmp(arg[0],"massflow_correction_scale") == 0) {
      if (narg < 3) error->fix_error(FLERR,this,"Illegal fix_modify command");
      FixScaleDiameter * fix_scale_ = static_cast<FixScaleDiameter*>(modify->find_fix_id_style(arg[1], "scale/diameter"));
      ++nusedarg;
      if(fix_scale_) {
        fix_scale_->set_scale_mass(false);
        int start_id = atoi(arg[2]);
        int end_id = start_id;
        if (narg > 3) {
          ++nusedarg;
          end_id = atoi(arg[3]);
        }
        while (start_id <= end_id) {
          if (actual_->has_cell_id(start_id))
            modifier_scale_[fix_scale_].insert(actual_->cell(start_id));
          ++start_id;
        }
      }
      return nusedarg;
    }

    int start_id = atoi(arg[1]);
    int end_id = start_id;
    if (narg > 2) {
      ++nusedarg;
      end_id = atoi(arg[2]);
    }

    if (strcmp(arg[0],"activate") == 0) {
      while (start_id <= end_id) {
        if (actual_->has_cell_id(start_id))
          active_.insert(actual_->cell(start_id));
        ++start_id;
      }
    } else if (strcmp(arg[0],"deactivate") == 0) {
      while (start_id <= end_id) {
        if (actual_->has_cell_id(start_id))
          active_.erase(actual_->cell(start_id));
        ++start_id;
      }
    } else if (strcmp(arg[0],"massflow_correction_on") == 0) {
      while (start_id <= end_id) {
        if (target_has_cell_id(start_id))
          modifier_[target_cell_index(start_id)] = true;
        ++start_id;
      }
    } else if (strcmp(arg[0],"massflow_correction_off") == 0) {
      while (start_id <= end_id) {
        if (target_has_cell_id(start_id))
          modifier_[target_cell_index(start_id)] = false;
        ++start_id;
      }
    }
    else error->all(FLERR,"Illegal fix_modify command");

    return nusedarg;
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::reset_dt()
{
  dtv_ = update->dt;
  dtv_inverse_ = 1.0/update->dt;
  dtf_ = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixForceControlRegion::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixForceControlRegion::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into restart file
------------------------------------------------------------------------- */

void FixForceControlRegion::write_restart(FILE *fp)
{
  std::vector<double> used_part_cell_all(ncells_max_,0.0);
  std::vector<double> old_pv_vec_all(3*ncells_max_,0.0);
  std::vector<double> sum_err_all(3*ncells_max_,0.0);

  MPI_Reduce(used_part_cell_,      &used_part_cell_all[0], ncells_max_, MPI_DOUBLE, MPI_MAX,       0, world);
  MPI_Reduce(&(old_pv_vec_[0][0]), &old_pv_vec_all[0],   3*ncells_max_, MPI_DOUBLE, MPI_ABSMAX_OP, 0, world);
  MPI_Reduce(&(sum_err_[0][0]),    &sum_err_all[0],      3*ncells_max_, MPI_DOUBLE, MPI_ABSMAX_OP, 0, world);

  if (comm->me == 0) {
    int size = (1 + ncells_max_ * (1+3+3)) * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    double ncells = static_cast<double>(ncells_max_);
    fwrite(&ncells,sizeof(double),1,fp);
    fwrite(&used_part_cell_all[0],sizeof(double),  ncells_max_,fp);
    fwrite(&old_pv_vec_all[0],    sizeof(double),3*ncells_max_,fp);
    fwrite(&sum_err_all[0],       sizeof(double),3*ncells_max_,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixForceControlRegion::restart(char *buf)
{
  double *list = (double *) buf;

  int ncells = static_cast<int>(list[0]);
  if (ncells != actual_->ncells())
    error->fix_error(FLERR, this, "number of cells has changed upon restart");

  ++list;

  if (ncells > ncells_max_) {
    memory->grow(xvalue,ncells,"forcecontrol/region:xvalue");
    memory->grow(yvalue,ncells,"forcecontrol/region:yvalue");
    memory->grow(zvalue,ncells,"forcecontrol/region:zvalue");
    memory->grow(const_part_cell_,ncells,"forcecontrol/region:const_part_cell_");
    memory->grow(used_part_cell_,ncells,"forcecontrol/region:used_part_cell_");
    memory->grow(sinesq_part_cell_,ncells,"forcecontrol/region:sinesq_part_cell_");
    memory->grow(ncontrolled_cell_,ncells,"forcecontrol/region:ncontrolled_cell_");
    memory->grow(old_pv_vec_,ncells,3,"forcecontrol/region:old_pv_mag_");
    memory->grow(sum_err_,ncells,3,"forcecontrol/region:sum_err_");

    for (; ncells_max_ < ncells; ++ncells_max_) {
      xvalue[ncells_max_] = 0.0;
      yvalue[ncells_max_] = 0.0;
      zvalue[ncells_max_] = 0.0;
      const_part_cell_[ncells_max_] = const_part_;
      used_part_cell_[ncells_max_] = used_part_;
      sinesq_part_cell_[ncells_max_] = sinesq_part_;
      ncontrolled_cell_[ncells_max_] = -1;

      for (int i = 0; i < 3; ++i) {
        old_pv_vec_[ncells_max_][i] = 0.0;
        sum_err_[ncells_max_][i] = 0.0;
      }
    }
  }

  for (int icell=0; icell<ncells_max_; ++icell) {
    used_part_cell_[icell] = list[icell];
    const_part_cell_[icell] = used_part_cell_[icell] * CONST_TO_USED_PART_RATIO;
    sinesq_part_cell_[icell] = used_part_cell_[icell] - const_part_cell_[icell];
    for (int i=0; i<3; ++i) {
      old_pv_vec_[icell][i] = list[  ncells_max_+3*icell+i];
      sum_err_   [icell][i] = list[4*ncells_max_+3*icell+i];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::target_couple_every()
{
  return target_->nevery;
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::target_has_cell_id(int cell_id)
{
  return target_->has_cell_id(cell_id);
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::target_cell_index(int cell_id)
{
  return target_->cell(cell_id);
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::target_cell_count(int cell_index)
{
  return target_->cell_count(cell_index);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_mass(int cell_index)
{
  return target_->cell_mass(cell_index);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_vol_fr(int cell_index)
{
  return target_->cell_vol_fr(cell_index);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_stress_xx(int cell_index)
{
  return target_->cell_stress(cell_index, 0);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_stress_yy(int cell_index)
{
  return target_->cell_stress(cell_index, 1);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_stress_zz(int cell_index)
{
  return target_->cell_stress(cell_index, 2);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_ave_x(int cell_index)
{
  return target_->cell_v_av(cell_index, 0);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_ave_y(int cell_index)
{
  return target_->cell_v_av(cell_index, 1);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_ave_z(int cell_index)
{
  return target_->cell_v_av(cell_index, 2);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_min_x(int cell_index)
{
  return target_->cell_v_min(cell_index, 0);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_min_y(int cell_index)
{
  return target_->cell_v_min(cell_index, 1);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_min_z(int cell_index)
{
  return target_->cell_v_min(cell_index, 2);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_max_x(int cell_index)
{
  return target_->cell_v_max(cell_index, 0);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_max_y(int cell_index)
{
  return target_->cell_v_max(cell_index, 1);
}

/* ---------------------------------------------------------------------- */

double FixForceControlRegion::target_cell_v_max_z(int cell_index)
{
  return target_->cell_v_max(cell_index, 2);
}

#endif
