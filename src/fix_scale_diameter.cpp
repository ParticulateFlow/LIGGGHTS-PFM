/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_scale_diameter.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "error.h"
#include "fix_property_atom.h"
#include "region.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM,EXTERNAL};
enum{XCYLINDER, YCYLINDER, ZCYLINDER, SPHERE};

/* ---------------------------------------------------------------------- */

FixScaleDiameter::FixScaleDiameter(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  fix_property_(NULL),
  scale_region(NULL),
  idregion(NULL),
  radius_(1.0),
  radius_inner_(0.0),
  scale_to_(1.0),
  scale_to_var_(-1),
  scale_to_style_(NONE),
  scale_to_str_(NULL),
  scale_range_(1.0),
  scale_width_(0.0),
  scale_mass_(true)
{
  if (narg < 5) error->fix_error(FLERR,this,"Not enough arguments for fix scale/diameter command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->fix_error(FLERR,this,"Illegal fix scale/diameter command");

  int iarg = 4;
  bool hasargs = true;
  while(iarg < narg && hasargs) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"Not enough arguments for 'region' option");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->fix_error(FLERR,this,"Region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      scale_region = domain->regions[iregion];
      if(strcmp(scale_region->style,"sphere") == 0)
        region_style = SPHERE;
      else if(strcmp(scale_region->style,"cylinder") == 0)
        region_style = XCYLINDER;
      else
        error->fix_error(FLERR,this,"Region must be of style 'sphere' or 'cylinder'");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->fix_error(FLERR,this,"Illegal fix scale/diameter command");
      ++iarg;
      if (strstr(arg[iarg],"v_") == arg[iarg]) {
        const int n = strlen(&arg[iarg][2]) + 1;
        scale_to_str_ = new char[n];
        strcpy(scale_to_str_,&arg[iarg][2]);
        scale_to_style_ = EQUAL;
      } else {
        scale_to_ = force->numeric(FLERR,arg[iarg]);
        if(scale_to_ <= 0.)
          error->fix_error(FLERR,this,"scale > 0 required");
        scale_to_style_ = CONSTANT;
        scale_range_ = 1. - scale_to_;
      }
      ++iarg;
      scale_width_ = force->numeric(FLERR,arg[iarg]);
      if(scale_width_ > 1.0 || scale_width_ <= 0.0)
        error->fix_error(FLERR,this,"scale_width must be in range ]0;1]");
      ++iarg;
      hasargs = true;
    } else break;
  }


  if(!scale_region)
    error->fix_error(FLERR,this,"No region defined");

  if(region_style == SPHERE) {
    radius_ = 0.5*(scale_region->extent_xhi - scale_region->extent_xlo);
    center_[0] = scale_region->extent_xhi - radius_;
    center_[1] = scale_region->extent_yhi - radius_;
    center_[2] = scale_region->extent_zhi - radius_;
  } else {
    center_[0] = 0.5*(scale_region->extent_xhi + scale_region->extent_xlo);
    center_[1] = 0.5*(scale_region->extent_yhi + scale_region->extent_ylo);
    center_[2] = 0.5*(scale_region->extent_zhi + scale_region->extent_zlo);

    double dim[3] = {scale_region->extent_xhi - scale_region->extent_xlo,
                     scale_region->extent_yhi - scale_region->extent_ylo,
                     scale_region->extent_zhi - scale_region->extent_zlo};
    if(fabs(dim[0] - dim[1]) < 1e-6) {
      radius_ = 0.5*dim[0];
      region_style = ZCYLINDER;
    } else if(fabs(dim[0] - dim[2]) < 1e-6) {
      radius_ = 0.5*dim[0];
      region_style = YCYLINDER;
    } else {
      radius_ = 0.5*dim[1];
    }
  }

  rad_mass_vary_flag = 1;

  scale_width_ *= radius_;
  radius_inner_ = radius_ - scale_width_;

  fix_property_ = NULL;
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

FixScaleDiameter::~FixScaleDiameter()
{
  delete []idregion;
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::post_create()
{
  if (fix_property_ == NULL) {
    char *fixid = new char[14+strlen(id)];
    sprintf(fixid,"scale_d_orig_%s",id);

    const char *fixarg[9];
    // register property/atom
    fixarg[0]=fixid;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=fixid;
    fixarg[4]="scalar"; // 1 scalar per particle to be registered
    fixarg[5]="yes";    // restart yes
    fixarg[6]="no";     // communicate ghost no
    fixarg[7]="no";     // communicate rev no
    fixarg[8]="0.";     // take 0 as default radius
    modify->add_fix(9,const_cast<char**>(fixarg));

    fix_property_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid,"property/atom","scalar",0,0,style));
    delete [] fixid;
  }

  const double *radius = atom->radius;
  const int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; ++i) {
    fix_property_->vector_atom[i] = radius[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::pre_delete(bool unfixflag)
{
  if (unfixflag && fix_property_) modify->delete_fix(fix_property_->id);
}

/* ---------------------------------------------------------------------- */

int FixScaleDiameter::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::init()
{
  if (scale_to_str_) {
    scale_to_var_ = input->variable->find(scale_to_str_);
    if (scale_to_var_ < 0)
      error->fix_error(FLERR,this,"Variable name does not exist");
    if (!input->variable->equalstyle(scale_to_var_))
      error->fix_error(FLERR,this,"Variable is invalid style");
  }
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::setup_pre_force(int vflag)
{
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::pre_force(int vflag)
{
  if (nevery == 0) return;
  if (update->ntimestep % nevery) return;
  change_settings();
}

/* ---------------------------------------------------------------------- */

void FixScaleDiameter::set_scale(double scale)
{
  scale_to_= scale;
  scale_range_ = 1. - scale_to_;
  scale_to_style_ = EXTERNAL;
}

/* ----------------------------------------------------------------------
   change atom parameters based on distance to region center
------------------------------------------------------------------------- */

void FixScaleDiameter::change_settings()
{
  double **x = atom->x;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  const int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  if (scale_to_style_ == EQUAL) {
    scale_to_ = input->variable->compute_equal(scale_to_var_);
    if(scale_to_ <= 0.)
      error->fix_error(FLERR,this,"scale_to_ > 0 required");
    scale_range_ = 1. - scale_to_;
  }

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      if (!scale_region->match(x[i][0],x[i][1],x[i][2]))
        continue;

      const double delx = center_[0] - x[i][0];
      const double dely = center_[1] - x[i][1];
      const double delz = center_[2] - x[i][2];

      double rsq = 0.;

      switch(region_style) {
      case SPHERE:
        rsq = delx * delx + dely * dely + delz * delz;
        break;
      case XCYLINDER:
        rsq = dely * dely + delz * delz;
        break;
      case YCYLINDER:
        rsq = delx * delx + delz * delz;
        break;
      case ZCYLINDER:
        rsq = delx * delx + dely * dely;
        break;
      }

      double scale = 1.0;

      if (rsq < radius_inner_*radius_inner_) {
        if (scale_to_style_ == CONSTANT)
          continue;
        scale = scale_to_;
      } else {
        scale = scale_to_ + scale_range_*((sqrt(rsq)-radius_inner_)/scale_width_);
      }

      const double old_radius = radius[i];

      if (fix_property_->vector_atom[i] <= 0.0) {
        fix_property_->vector_atom[i] = old_radius;
      }

      radius[i] = scale * fix_property_->vector_atom[i];

      if (scale_mass_) {
        const double relative_scale = radius[i]/old_radius;
        rmass[i] *= relative_scale*relative_scale*relative_scale;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixScaleDiameter::compute_scalar()
{
  return scale_to_;
}
