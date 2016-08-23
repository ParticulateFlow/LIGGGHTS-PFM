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

#include "string.h"
#include "stdlib.h"
#include "fix_ave_euler_region.h"
#include "fix_forcecontrol_region.h"
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

/* ---------------------------------------------------------------------- */

FixForceControlRegion::FixForceControlRegion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xvalue(NULL),
  yvalue(NULL),
  zvalue(NULL),
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
  cg_(1),
  ncells_max_(0),
  old_pv_vec_(NULL),
  sum_err_(NULL),
  const_part_(1.0),
  sinesq_part_(0.0),
  used_part_(1.0),
  acceptable_deviation_min(0.97),
  acceptable_deviation_max(1.03),
  limit_velocity_(false)
{
  if (narg < 6) error->all(FLERR,"Illegal fix addforce command");

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
      int ifix = modify->find_fix(arg[iarg]);

      if(ifix < 0)
          error->all(FLERR,"Illegal fix forcecontrol/region command, invalid ID for fix ave/euler/region provided");

      if(strncmp(modify->fix[ifix]->style,"ave/euler/region",16))
          error->all(FLERR,"Illegal fix forcecontrol/region command, fix is not of type ave/euler/region");

      target_ = static_cast<FixAveEulerRegion*>(modify->fix[ifix]);

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
    } else if (strcmp(arg[iarg],"cg") == 0) { // TODO: remove when cg and fg are separate simulations -> cg_ = force->cg();
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      cg_ = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal fix forcecontrol/region command");
  }

  if(!actual_)
    error->fix_error(FLERR, this, "missing ave/euler/region for actual_val");
  if(!target_)
    error->fix_error(FLERR, this, "missing ave/euler/region for target_val");

  cg3_ = cg_*cg_*cg_;
  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixForceControlRegion::~FixForceControlRegion()
{
  memory->destroy(xvalue);
  memory->destroy(yvalue);
  memory->destroy(zvalue);
  memory->destroy(old_pv_vec_);
  memory->destroy(sum_err_);
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::post_create()
{
  // all cells active by default
  int ncells = actual_->ncells();
  for (int icell=0; icell<ncells; ++icell) {
    active_.insert(icell);
  }
  modifier_.resize(ncells, false);

  double maxrd,minrd;
  modify->max_min_rad(maxrd,minrd);
  if (ctrl_style_ == STRESS)
  {
    const_part_ = (2.*maxrd/cg_)*1.2; // TODO: change when cg and fg are separate simulations
    used_part_ = (2.*maxrd/cg_)*1.4;
    sinesq_part_ = used_part_ - const_part_;
  }
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
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
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::post_force(int vflag)
{
  UNUSED(vflag);
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
    memory->grow(old_pv_vec_,ncells,3,"forcecontrol/region:old_pv_mag_");
    memory->grow(sum_err_,ncells,3,"forcecontrol/region:sum_err_");

    for (; ncells_max_ < ncells; ++ncells_max_) {
      xvalue[ncells_max_] = 0.0;
      yvalue[ncells_max_] = 0.0;
      zvalue[ncells_max_] = 0.0;

      for (int i = 0; i < 3; ++i) {
          old_pv_vec_[ncells_max_][i] = 0.0;
          sum_err_[ncells_max_][i] = 0.0;
      }
    }
  }

  std::set<int>::iterator it_cell = active_.begin();
  for (; it_cell!=active_.end(); ++it_cell) { // active actual_ cell indices

    int cell_id = actual_->cell_id(*it_cell);

    // skip cell if target_ cell has less than half of equivalent cg particle
    // mostly this happens when a stream of particles enters an empty cell
    if (target_->compute_array_by_id(cell_id,3) <= actual_->cell_vol_fr(*it_cell)*0.5) {
      continue;
    }

    int tcell = target_->cell(cell_id);

    // clear controller history for empty cells
    if (actual_->cell_count(*it_cell) < 1 || target_->cell_count(tcell) < 1) {
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
        // get stress directions from region mesh/hex (vector properties stored in vtk file)
        double rsq = actual_->cell_radius(*it_cell)*actual_->cell_radius(*it_cell); // fg radius squared
        double pre = (2.*M_PI*rsq/3.);//*actual_->cell_volume(*it_cell);

        pv_vec_[0] = -(pre*actual_->cell_stress(*it_cell, 0)/actual_->cell_count(*it_cell)); // fg
        pv_vec_[1] = -(pre*actual_->cell_stress(*it_cell, 1)/actual_->cell_count(*it_cell)); // fg
        pv_vec_[2] = -(pre*actual_->cell_stress(*it_cell, 2)/actual_->cell_count(*it_cell)); // fg
        sp_vec_[0] = -(pre*target_->cell_stress(tcell, 0)/(target_->cell_count(tcell)*cg3_));// cg
        sp_vec_[1] = -(pre*target_->cell_stress(tcell, 1)/(target_->cell_count(tcell)*cg3_));// cg
        sp_vec_[2] = -(pre*target_->cell_stress(tcell, 2)/(target_->cell_count(tcell)*cg3_));// cg

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
          massflow_correction = target_->cell_vol_fr(tcell)/actual_->cell_vol_fr(*it_cell);

        sp_vec_[0] = massflow_correction*target_->cell_v_av(tcell, 0);
        sp_vec_[1] = massflow_correction*target_->cell_v_av(tcell, 1);
        sp_vec_[2] = massflow_correction*target_->cell_v_av(tcell, 2);

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
      err_[0] *= axis_[0];
      err_[1] *= axis_[1];
      err_[2] *= axis_[2];
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
      xvalue[*it_cell] = (axis_[0] < 0.) ? std::min(0., ctrl_op_[0]) : std::max(0., ctrl_op_[0]);
      yvalue[*it_cell] = (axis_[1] < 0.) ? std::min(0., ctrl_op_[1]) : std::max(0., ctrl_op_[1]);
      zvalue[*it_cell] = (axis_[2] < 0.) ? std::min(0., ctrl_op_[2]) : std::max(0., ctrl_op_[2]);
    } else {
      xvalue[*it_cell] = ctrl_op_[0];
      yvalue[*it_cell] = ctrl_op_[1];
      zvalue[*it_cell] = ctrl_op_[2];
    }

    // TODO: remove
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
    force_flag = 0;

    if (xvalue[*it_cell]!=0. || yvalue[*it_cell]!=0. || zvalue[*it_cell]!=0.) {
      double bounds[6];
      actual_->cell_bounds(*it_cell, bounds);

      vel_min_[0] = target_->cell_v_min(tcell, 0);
      vel_min_[1] = target_->cell_v_min(tcell, 1);
      vel_min_[2] = target_->cell_v_min(tcell, 2);
      vel_max_[0] = target_->cell_v_max(tcell, 0);
      vel_max_[1] = target_->cell_v_max(tcell, 1);
      vel_max_[2] = target_->cell_v_max(tcell, 2);

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

          if (ctrl_style_ == VELOCITY) {
            fadex_ = fadey_ = fadez_ = rmass[i]*dtv_inverse_;
          } else if (ctrl_style_ == STRESS  && sinesq_part_ > 0.0) {
            if (axis_[0] != 0.) {
              if (axis_[0] < 0. && x[i][0] < bounds[1] - const_part_) {
                if(x[i][0] < bounds[1] - used_part_) {
                  fadex_ = 0.;
                } else {
                  fadex_ = sin(M_PI*0.5*(x[i][0] - (bounds[1]-used_part_))/sinesq_part_);
                  fadex_ *= fadex_;
                }
              } else if (axis_[0] > 0. && x[i][0] > bounds[0] + const_part_) {
                if(x[i][0] > bounds[0] + used_part_) {
                  fadex_ = 0.;
                } else {
                  fadex_ = sin(M_PI*0.5*(x[i][0] - (bounds[0]+used_part_))/sinesq_part_);
                  fadex_ *= fadex_;
                }
              }
            }

            if (axis_[1] != 0.) {
              if (axis_[1] < 0. && x[i][1] < bounds[3] - const_part_) {
                if(x[i][1] < bounds[3] - used_part_) {
                  fadey_ = 0.;
                } else {
                  fadey_ = sin(M_PI*0.5*(x[i][1] - (bounds[3]-used_part_))/sinesq_part_);
                  fadey_ *= fadey_;
                }
              } else if (axis_[1] > 0. && x[i][1] > bounds[2] + const_part_) {
                if(x[i][1] > bounds[2] + used_part_) {
                  fadey_ = 0.;
                } else {
                  fadey_ = sin(M_PI*0.5*(x[i][1] - (bounds[2]+used_part_))/sinesq_part_);
                  fadey_ *= fadey_;
                }
              }
            }

            if (axis_[2] != 0.) {
              if (axis_[2] < 0. && x[i][2] < bounds[5] - const_part_) {
                if(x[i][2] < bounds[5] - used_part_) {
                  fadez_ = 0.;
                } else {
                  fadez_ = sin(M_PI*0.5* (x[i][2] - (bounds[5]-used_part_))/sinesq_part_);
                  fadez_ *= fadez_;
                }
              } else if (axis_[2] > 0. && x[i][2] > bounds[4] + const_part_) {
                if(x[i][2] > bounds[4] + used_part_) {
                  fadez_ = 0.;
                } else {
                  fadez_ = sin(M_PI*0.5* (x[i][2] - (bounds[4]+used_part_))/sinesq_part_);
                  fadez_ *= fadez_;
                }
              }
            }
          }


          if (limit_velocity_ && !modifier_[tcell]) {
            limit[0] = limit[1] = limit[2] = 1.0;
            double dtfm = dtf_ / rmass[i];

            fx = fadex_ * xvalue[*it_cell];

            if (fabs(fx) > fabs(f[i][0]*0.005)) {
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

            if (fabs(fy) > fabs(f[i][1]*0.005)) {
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

            if (fabs(fz) > fabs(f[i][2]*0.005)) {
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
          } else {
            f[i][0] += fadex_ * xvalue[*it_cell];
            f[i][1] += fadey_ * yvalue[*it_cell];
            f[i][2] += fadez_ * zvalue[*it_cell];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::modify_param(int narg, char **arg)
{
    if (narg < 2) error->fix_error(FLERR,this,"Illegal fix_modify command");
    int nusedarg = 2;
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
        if (target_->has_cell_id(start_id))
          modifier_[target_->cell(start_id)] = true;
        ++start_id;
      }
    } else if (strcmp(arg[0],"massflow_correction_off") == 0) {
      while (start_id <= end_id) {
        if (target_->has_cell_id(start_id))
          modifier_[target_->cell(start_id)] = false;
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
