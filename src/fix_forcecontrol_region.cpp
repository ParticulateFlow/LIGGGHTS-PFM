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

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,STRESS,VELOCITY};

/* ---------------------------------------------------------------------- */

FixForceControlRegion::FixForceControlRegion(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  kp_(0.01),
  ki_(0.),
  kd_(0.),
  ctrl_style_(NONE),
  cg_(1),
  ncells_max_(0),
  old_pv_vec_(NULL),
  const_part_(0.0)
{
  if (narg < 6) error->all(FLERR,"Illegal fix addforce command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  xvalue = yvalue = zvalue = 0.;

  vectorZeroize3D(axis_);
  vectorZeroize3D(ctrl_op_);

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ctrlPV") == 0) {
      if (narg < iarg+2) error->fix_error(FLERR,this,"not enough arguments for 'ctrlPV'");
      if (strcmp(arg[iarg+1],"stress") == 0) {
        ctrl_style_ = STRESS;
        const_part_ = 0.75;
      }
      else if (strcmp(arg[iarg+1],"velocity") == 0) {
        ctrl_style_ = VELOCITY;
        const_part_ = 1.0;
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
    } else if(strcmp(arg[iarg],"axis") == 0) {
      if (narg < iarg+4) error->fix_error(FLERR,this,"not enough arguments for 'axis'");
      axis_[0] = force->numeric(FLERR,arg[iarg+1]);
      axis_[1] = force->numeric(FLERR,arg[iarg+2]);
      axis_[2] = force->numeric(FLERR,arg[iarg+3]);
      // normalize axis
      vectorNormalize3D(axis_);
      iarg = iarg+4;
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

  sinesq_part_ = 1.0 - const_part_;
  cg3_ = cg_*cg_*cg_;
  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixForceControlRegion::~FixForceControlRegion()
{
  memory->destroy(old_pv_vec_);
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

  double vx, vy, vz;
  int ncells = actual_->ncells();

  // reallocate old_pv_vec_ array if necessary
  if (ncells > ncells_max_) {
    memory->grow(old_pv_vec_,ncells,3,"forcecontrol/region:old_pv_mag_");
    for (; ncells_max_ < ncells; ++ncells_max_) {
      for (int i = 0; i < 3; ++i)
        old_pv_vec_[ncells_max_][i] = 0.0;
    }
  }


  std::set<int>::iterator it_cell = active_.begin();
  for (; it_cell!=active_.end(); ++it_cell) { // active actual_ cell indices

    int cell_id = actual_->cell_id(*it_cell);
    // skip cell if target_ cell has less than half of equivalent cg particle
    // mostly this happens when a stream of particles enters an empty cell
    if (target_->compute_array_by_id(cell_id,3) < actual_->cell_vol_fr(*it_cell)*0.5)
      continue;

    int tcell = target_->cell(cell_id);

    switch (ctrl_style_) {
    case STRESS:
      {
        // stress target is supposed to come from higher cg level
        double rsq = actual_->cell_radius(*it_cell)*actual_->cell_radius(*it_cell);
        double pre = (2.*M_PI*rsq/3.)*actual_->cell_volume(*it_cell);
        pv_vec_[0] = pre*actual_->cell_stress(*it_cell, 0)/actual_->cell_count(*it_cell); // fg
        pv_vec_[1] = pre*actual_->cell_stress(*it_cell, 1)/actual_->cell_count(*it_cell); // fg
        pv_vec_[2] = pre*actual_->cell_stress(*it_cell, 2)/actual_->cell_count(*it_cell); // fg
        sp_vec_[0] = pre*target_->cell_stress(tcell, 0)/(target_->cell_count(tcell)*cg3_);// cg
        sp_vec_[1] = pre*target_->cell_stress(tcell, 1)/(target_->cell_count(tcell)*cg3_);// cg
        sp_vec_[2] = pre*target_->cell_stress(tcell, 2)/(target_->cell_count(tcell)*cg3_);// cg
        // get stress directions from region mesh/hex (vector properties stored in vtk file)
        double *stress_ctrl_dir = actual_->cell_vector_property(*it_cell, "stress_ctrl_dir");
        axis_[0] = stress_ctrl_dir[0];
        axis_[1] = stress_ctrl_dir[1];
        axis_[2] = stress_ctrl_dir[2];
        break;
      }
    case VELOCITY:
      {
        pv_vec_[0] = -actual_->cell_v_av(*it_cell, 0);
        pv_vec_[1] = -actual_->cell_v_av(*it_cell, 1);
        pv_vec_[2] = -actual_->cell_v_av(*it_cell, 2);
        double massflow_correction = 1.;
        if (modifier_[tcell] && actual_->cell_vol_fr(*it_cell) > 0.)
          massflow_correction = target_->cell_vol_fr(tcell)/actual_->cell_vol_fr(*it_cell);
        sp_vec_[0] = -massflow_correction*target_->cell_v_av(tcell, 0);
        sp_vec_[1] = -massflow_correction*target_->cell_v_av(tcell, 1);
        sp_vec_[2] = -massflow_correction*target_->cell_v_av(tcell, 2);
        axis_[0] = axis_[1] = axis_[2] = 1.;
        break;
      }
    }

    // simple PID-controller

    // calc error and sum of the errors
    err_[0] = sp_vec_[0] - pv_vec_[0];
    err_[1] = sp_vec_[1] - pv_vec_[1];
    err_[2] = sp_vec_[2] - pv_vec_[2];
    sum_err_[0] += err_[0] * dtv_;
    sum_err_[1] += err_[1] * dtv_;
    sum_err_[2] += err_[2] * dtv_;
    // derivative term
    // force used instead of the error in order to avoid signal spikes in case of change of the set point
    double dfdt[3];
    dfdt[0] = -( pv_vec_[0] - old_pv_vec_[*it_cell][0])/dtv_;
    dfdt[1] = -( pv_vec_[1] - old_pv_vec_[*it_cell][1])/dtv_;
    dfdt[2] = -( pv_vec_[2] - old_pv_vec_[*it_cell][2])/dtv_;

    // vel points opposite to force vector
    ctrl_op_[0] =  err_[0] * kp_ + sum_err_[0] * ki_ + dfdt[0] * kd_;
    ctrl_op_[1] =  err_[1] * kp_ + sum_err_[1] * ki_ + dfdt[1] * kd_;
    ctrl_op_[2] =  err_[2] * kp_ + sum_err_[2] * ki_ + dfdt[2] * kd_;
    // save process value for next timestep
    old_pv_vec_[*it_cell][0] = pv_vec_[0];
    old_pv_vec_[*it_cell][1] = pv_vec_[1];
    old_pv_vec_[*it_cell][2] = pv_vec_[2];

    xvalue+=ctrl_op_[0];
    yvalue+=ctrl_op_[1];
    zvalue+=ctrl_op_[2];

    // TODO: remove
    foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
    force_flag = 0;

    vel_min_[0] = target_->cell_v_min(tcell, 0);
    vel_min_[1] = target_->cell_v_min(tcell, 1);
    vel_min_[2] = target_->cell_v_min(tcell, 2);
    vel_max_[0] = target_->cell_v_max(tcell, 0);
    vel_max_[1] = target_->cell_v_max(tcell, 1);
    vel_max_[2] = target_->cell_v_max(tcell, 2);

    double bounds[6];
    actual_->cell_bounds(*it_cell, bounds);

    if (xvalue!=0. || yvalue!=0. || zvalue!=0.) {
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {

          double fadex = 1.0;
          double fadey = 1.0;
          double fadez = 1.0;


          if (axis_[0] != 0.) {
            double extent_x = bounds[1] - bounds[0];

            if (axis_[0] < 0. && x[i][0] < bounds[1] - const_part_*extent_x) {
              fadex = sin(M_PI*0.5*(x[i][0] - bounds[0])/(sinesq_part_*extent_x));
              fadex *= fadex;
            } else if (axis_[0] > 0. && x[i][0] > bounds[0] + const_part_*extent_x) {
              fadex = sin(M_PI*0.5*(x[i][0] - bounds[1])/(sinesq_part_*extent_x));
              fadex *= fadex;
            }
          }
          if (axis_[1] != 0.) {
            double extent_y = bounds[3] - bounds[2];

            if (axis_[1] < 0. && x[i][1] < bounds[3] - const_part_*extent_y) {
              fadey = sin(M_PI*0.5*(x[i][1] - bounds[2])/(sinesq_part_*extent_y));
              fadey *= fadey;
            } else if (axis_[1] > 0. && x[i][1] > bounds[2] + const_part_*extent_y) {
              fadey = sin(M_PI*0.5*(x[i][1] - bounds[3])/(sinesq_part_*extent_y));
              fadey *= fadey;
            }
          }
          if (axis_[2] != 0.) {
            double extent_z = bounds[5] - bounds[4];

            if (axis_[2] < 0. && x[i][2] < bounds[5] - const_part_*extent_z) {
              fadez = sin(M_PI*0.5*(x[i][2] - bounds[4])/(sinesq_part_*extent_z));
              fadez *= fadez;
            } else if (axis_[2] > 0. && x[i][2] > bounds[4] + const_part_*extent_z) {
              fadez = sin(M_PI*0.5*(x[i][2] - bounds[5])/(sinesq_part_*extent_z));
              fadez *= fadez;
            }
          }

          double dtfm = dtf_ / rmass[i];
          vx = v[i][0] + dtfm * (f[i][0] + fadex * xvalue);
          vy = v[i][1] + dtfm * (f[i][1] + fadey * yvalue);
          vz = v[i][2] + dtfm * (f[i][2] + fadez * zvalue);

          double limit[3];
          limit[0]=limit[1]=limit[2]=1.;

          if (vx > vel_max_[0])
            limit[0] = (fadex*xvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_max_[0] - v[i][0]) / dtfm - f[i][0])/(fadex*xvalue))) );
          else if (vx < vel_min_[0])
            limit[0] = (fadex*xvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_min_[0] - v[i][0]) / dtfm - f[i][0])/(fadex*xvalue))) );

          if (vy > vel_max_[1])
            limit[1] = (fadey*yvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_max_[1] - v[i][1]) / dtfm - f[i][1])/(fadey*yvalue))) );
          else if (vy < vel_min_[1])
            limit[1] = (fadey*yvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_min_[1] - v[i][1]) / dtfm - f[i][1])/(fadey*yvalue))) );

          if (vz > vel_max_[2])
            limit[2] = (fadez*zvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_max_[2] - v[i][2]) / dtfm - f[i][2])/(fadez*zvalue))) );
          else if (vz < vel_min_[2])
            limit[2] = (fadez*zvalue == 0.) ? 1. : ( std::max(0., std::min(1., ((vel_min_[2] - v[i][2]) / dtfm - f[i][2])/(fadez*zvalue))) );

          double a = std::min(limit[0], std::min(limit[1], limit[2]));
          f[i][0] += a * fadex * xvalue;
          f[i][1] += a * fadey * yvalue;
          f[i][2] += a * fadez * zvalue;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixForceControlRegion::modify_param(int narg, char **arg)
{
    if (narg < 2) error->fix_error(FLERR,this,"Illegal fix_modify command");
    int start_id = atoi(arg[1]);
    int end_id = start_id;
    if (narg > 2) end_id = atoi(arg[2]);

    if (strcmp(arg[0],"activate") == 0) {
      while (start_id <= end_id) {
        active_.insert(actual_->cell(start_id));
        ++start_id;
      }
    } else if (strcmp(arg[0],"deactivate") == 0) {
      while (start_id <= end_id) {
        active_.erase(actual_->cell(start_id));
        ++start_id;
      }
    } else if (strcmp(arg[0],"massflow_correction_on") == 0) {
      while (start_id <= end_id) {
        int tcell = target_->cell(start_id);
        modifier_[tcell] = true;
        ++start_id;
      }
    } else if (strcmp(arg[0],"massflow_correction_off") == 0) {
      while (start_id <= end_id) {
        int tcell = target_->cell(start_id);
        modifier_[tcell] = false;
        ++start_id;
      }
    }
    else error->all(FLERR,"Illegal fix_modify command");
    return 2;
}

/* ---------------------------------------------------------------------- */

void FixForceControlRegion::reset_dt()
{
  dtv_ = update->dt;
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
