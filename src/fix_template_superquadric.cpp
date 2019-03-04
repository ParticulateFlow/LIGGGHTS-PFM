/* ----------------------------------------------------------------------
    This file is part of the LIGGGHTS DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#include <math.h>
#include "math_extra.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_template_superquadric.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsertSuperquadric.h"
#include "region.h"
#include "domain.h"
#include "force.h"
#include "comm.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_region_variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;
using namespace FixConst;

#define LMP_DEBUGMODE_SUPERQUADRIC false

/* ---------------------------------------------------------------------- */

FixTemplateSuperquadric::FixTemplateSuperquadric(LAMMPS *lmp, int narg, char **arg) :
    FixTemplateSphere(lmp, narg, arg)
{
  if (!atom->superquadric_flag)
      error->all(FLERR,"Fix particletemplate/superquadric requires atom style superquadric");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->fix_error(FLERR,this,"not enough arguments");

  iarg = 4;

  // set default values
  atom_type = 1;
  vol_limit = 1e-12;

  pdf_density = NULL;
  pdf_shapex = NULL;
  pdf_shapey = NULL;
  pdf_shapez = NULL;
  pdf_size = NULL;
  pdf_blockiness1 = new PDF(error);
  pdf_blockiness2 = new PDF(error);

  pdf_blockiness1->set_params<RANDOM_CONSTANT>(2.0);
  pdf_blockiness2->set_params<RANDOM_CONSTANT>(2.0);

  n_pti_max = 0;
  pti_list = NULL;

  reg = NULL;
  reg_var = NULL;

  relative = false;

  //parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"atom_type") == 0)
    {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments");
      atom_type=atoi(arg[iarg+1]);
      if (atom_type < 1)
        error->fix_error(FLERR,this,"invalid atom type (must be >=1)");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region") == 0)
    {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments");
      int ireg = domain->find_region(arg[iarg+1]);
      if (ireg < 0) error->fix_error(FLERR,this,"illegal region");
      reg = domain->regions[ireg];
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"relative") == 0)
    {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'relative'");
      if(0 == strcmp(arg[iarg+1],"yes"))
        relative = true;
      else if(0 == strcmp(arg[iarg+1],"no"))
        relative = false;
      else
        error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'relative'");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region_variable") == 0)
    {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0)
        error->fix_error(FLERR,this,"illegal region/variable fix");
      reg_var = static_cast<FixRegionVariable*>(modify->fix[ifix]);
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"volume_limit") == 0)
    {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'volume_limit'");
      vol_limit = atof(arg[iarg+1]);
      if(vol_limit <= 0)
        error->fix_error(FLERR,this,"volume_limit > 0 required");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"shape") == 0) {
      pdf_shapex = new PDF(error);
      pdf_shapey = new PDF(error);
      pdf_shapez = new PDF(error);
      hasargs = true;
      if(strcmp(this->style,"particletemplate/superquadric"))
        error->fix_error(FLERR,this,"keyword shape only valid for particletemplate/superquadric");
      if (iarg+5 > narg)
        error->all(FLERR,"Illegal fix particletemplate/superquadric command, not enough arguments");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
        double values[] = { atof(arg[iarg+2])*force->cg(),
                            atof(arg[iarg+3])*force->cg(),
                            atof(arg[iarg+4])*force->cg()};
        if( values[0] <= 0. || values[1] <= 0. || values[2] <= 0.)
          error->all(FLERR,"Illegal fix particletemplate/superquadric command, shape parameters must be >= 0");
        pdf_shapex->set_params<RANDOM_CONSTANT>(values[0]);
        pdf_shapey->set_params<RANDOM_CONSTANT>(values[1]);
        pdf_shapez->set_params<RANDOM_CONSTANT>(values[2]);
        iarg += 5;
      }
      else if (strcmp(arg[iarg+1],"gaussian") == 0)
      {
        if (iarg+6 > narg)
           error->fix_error(FLERR,this,"not enough arguments");
        double mux = atof(arg[iarg+2]);
        double muy = atof(arg[iarg+3]);
        double muz = atof(arg[iarg+4]);
        double sigma = atof(arg[iarg+5]);
        if(mux <= 0.)
          error->fix_error(FLERR,this,"illegal mux value for density");
        if(muy <= 0.)
          error->fix_error(FLERR,this,"illegal muy value for density");
        if(muz <= 0.)
          error->fix_error(FLERR,this,"illegal muz value for density");
        if( sigma <= 0. ) error->all(FLERR,"illegal sigmax value for density");
        pdf_shapex->set_params<RANDOM_GAUSSIAN>(mux,sigma);
        pdf_shapey->set_params<RANDOM_GAUSSIAN>(muy,sigma);
        pdf_shapez->set_params<RANDOM_GAUSSIAN>(muz,sigma);
        iarg += 6;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
        if (iarg+8 > narg)
           error->fix_error(FLERR,this,"not enough arguments");
        double shxmin = atof(arg[iarg+2]);
        double shxmax = atof(arg[iarg+3]);
        double shymin = atof(arg[iarg+4]);
        double shymax = atof(arg[iarg+5]);
        double shzmin = atof(arg[iarg+6]);
        double shzmax = atof(arg[iarg+7]);
        if(shxmax < shxmin)
          error->fix_error(FLERR,this,"max shapex less than min shapex");
        if(shymax < shymin)
          error->fix_error(FLERR,this,"max shapey less than min shapey");
        if(shzmax < shzmin)
          error->fix_error(FLERR,this,"max shapez less than min shapez");
        if(shxmin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapex");
        if(shymin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapey");
        if(shzmin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapez");

        pdf_shapex->set_params<RANDOM_UNIFORM>(shxmin,shxmax);
        pdf_shapey->set_params<RANDOM_UNIFORM>(shymin,shymax);
        pdf_shapez->set_params<RANDOM_UNIFORM>(shzmin,shzmax);
        iarg += 8;
      }
      else
          error->fix_error(FLERR,this,"fix particletemplate/superquadric currently supports only constant, gaussian and uniform shape params");
      }
      else if (strcmp(arg[iarg],"shapex") == 0 || strcmp(arg[iarg],"shapey") == 0 || strcmp(arg[iarg],"shapez") == 0) {
        pdf_shapex = new PDF(error);
        pdf_shapey = new PDF(error);
        pdf_shapez = new PDF(error);
        hasargs = true;
        if(strcmp(this->style,"particletemplate/superquadric"))
          error->fix_error(FLERR,this,"keyword shape only valid for particletemplate/superquadric");
        if (iarg+3 > narg)
          error->all(FLERR,"Illegal fix particletemplate/superquadric command, not enough arguments");
        if (strcmp(arg[iarg+1],"constant") == 0)
        {
          double value = atof(arg[iarg+2])*force->cg();
          if( value <= 0. )
            error->all(FLERR,"Illegal fix particletemplate/superquadric command, shape parameters must be >= 0");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_CONSTANT>(value);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_CONSTANT>(value);
          if(strcmp(arg[iarg],"shapez") == 0)
              pdf_shapez->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
        }
        else if (strcmp(arg[iarg+1],"gaussian") == 0)
        {
          if (iarg+4 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if(mu <= 0.)
            error->fix_error(FLERR,this,"illegal mux value for density");
          if( sigma <= 0. ) error->all(FLERR,"illegal sigma value for density");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          if(strcmp(arg[iarg],"shapez") == 0)
            pdf_shapez->set_params<RANDOM_GAUSSIAN>(mu,sigma);
          iarg += 4;
        }
        else if (strcmp(arg[iarg+1],"uniform") == 0)
        {
          if (iarg+4 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          double shmin = atof(arg[iarg+2]);
          double shmax = atof(arg[iarg+3]);
          if(shmax < shmin)
            error->fix_error(FLERR,this,"max shape less than min shape");
          if(shmin <= 0.0)
            error->fix_error(FLERR,this,"Illegal value min shapex");
          if(strcmp(arg[iarg],"shapex") == 0)
            pdf_shapex->set_params<RANDOM_UNIFORM>(shmin,shmax);
          if(strcmp(arg[iarg],"shapey") == 0)
            pdf_shapey->set_params<RANDOM_UNIFORM>(shmin,shmax);
          if(strcmp(arg[iarg],"shapez") == 0)
            pdf_shapez->set_params<RANDOM_UNIFORM>(shmin,shmax);
          iarg += 4;
        }
        else
          error->fix_error(FLERR,this,"fix particletemplate/superquadric currently supports only constant, gaussian and uniform shape params");
    }
    else if (strcmp(arg[iarg],"size") == 0) {
      pdf_size = new PDF(error);
      hasargs = true;
      if(strcmp(this->style,"particletemplate/superquadric"))
        error->fix_error(FLERR,this,"keyword shape only valid for particletemplate/superquadric");
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal fix particletemplate/superquadric command, not enough arguments");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
        double value = atof(arg[iarg+2])*force->cg();
        if( value <= 0. )
          error->all(FLERR,"Illegal fix particletemplate/superquadric command, shape parameters must be >= 0");
        pdf_size->set_params<RANDOM_CONSTANT>(value);
        iarg += 3;
      }
      else if (strcmp(arg[iarg+1],"gaussian") == 0)
      {
        if (iarg+4 > narg)
          error->fix_error(FLERR,this,"not enough arguments");
        double mu = atof(arg[iarg+2]);
        double sigma = atof(arg[iarg+3]);
        if(mu <= 0.)
          error->fix_error(FLERR,this,"illegal mux value for density");
        if( sigma <= 0. ) error->all(FLERR,"illegal sigma value for density");
        pdf_size->set_params<RANDOM_GAUSSIAN>(mu,sigma);
        iarg += 4;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
        if (iarg+4 > narg)
          error->fix_error(FLERR,this,"not enough arguments");
        double shmin = atof(arg[iarg+2]);
        double shmax = atof(arg[iarg+3]);
        if(shmax < shmin)
          error->fix_error(FLERR,this,"max shape less than min shape");
        if(shmin <= 0.0)
          error->fix_error(FLERR,this,"Illegal value min shapex");
        pdf_size->set_params<RANDOM_UNIFORM>(shmin,shmax);
        iarg += 4;
      }
      else
        error->fix_error(FLERR,this,"fix particletemplate/superquadric currently supports only constant, gaussian and uniform shape params");
    }
    else if (strcmp(arg[iarg],"blockiness") == 0  || strcmp(arg[iarg],"roundness") == 0) {
          hasargs = true;
          if(strcmp(this->style,"particletemplate/superquadric"))
            error->fix_error(FLERR,this,"keyword blockiness only valid for particletemplate/superquadric");
          if (iarg+4 > narg)
            error->all(FLERR,"Illegal fix particletemplate/superquadric command, not enough arguments");
          if (strcmp(arg[iarg+1],"constant") == 0)
          {
              if(strcmp(arg[iarg],"roundness") == 0)
                  error->warning(FLERR,"Keyword 'roundness' will be deprecated in future, please use 'blockiness' istead");
              double values[] = { atof(arg[iarg+2]),
                                  atof(arg[iarg+3])};
              if( values[0] < 2. || values[1] < 2.)
                error->all(FLERR,"Illegal fix particletemplate/superquadric command, blockiness parameters must >= 2");
              pdf_blockiness1->set_params<RANDOM_CONSTANT>(values[0]);
              pdf_blockiness2->set_params<RANDOM_CONSTANT>(values[1]);
              iarg += 4;
          }
          else
              error->fix_error(FLERR,this,"fix particletemplate/superquadric currently supports only constant blockiness params");
    }
    else if (strcmp(arg[iarg],"density") == 0) {
      hasargs = true;
      if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments");
      pdf_density = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2]);
          if( value <= 0.) error->fix_error(FLERR,this,"density must be >= 0");
          pdf_density->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else
          error->fix_error(FLERR,this,"fix particletemplate/superquadric currently supports only constant density");
    }
    else if(strcmp(style,"particletemplate/superquadric") == 0)
        error->fix_error(FLERR,this,"unrecognized keyword");
  }

  if(pdf_density == NULL) error->fix_error(FLERR,this,"have to define 'density'");

  // end here for derived classes
  if(strcmp(this->style,"particletemplate/superquadric"))return;

  if(pdf_size == NULL) {
    if(pdf_shapex == NULL) error->fix_error(FLERR,this,"have to define 'shape' or 'size'");
    if(pdf_shapey == NULL) error->fix_error(FLERR,this,"have to define 'shape' or 'size'");
    if(pdf_shapez == NULL) error->fix_error(FLERR,this,"have to define 'shape' or 'size'");
  } else {
    if(pdf_shapex != NULL) error->fix_error(FLERR,this,"'shape' and 'size' cannot be used simultaneously in the same fix. Use either 'size', or 'shape'");
    if(pdf_shapey != NULL) error->fix_error(FLERR,this,"'shape' and 'size' cannot be used simultaneously in the same fix. Use either 'size', or 'shape'");
    if(pdf_shapez != NULL) error->fix_error(FLERR,this,"'shape' and 'size' cannot be used simultaneously in the same fix. Use either 'size', or 'shape'");
  }

  // set mass and volume expectancy
  double blockiness_expect[] = { expectancy(pdf_blockiness1),
                                expectancy(pdf_blockiness2) };
  if(pdf_size == NULL) {
    double shape_expect[] = { expectancy(pdf_shapex),
                              expectancy(pdf_shapey),
                              expectancy(pdf_shapez) };

    double rad;
    MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape_expect, blockiness_expect, &rad);
    MathExtraLiggghtsNonspherical::volume_superquadric(shape_expect, blockiness_expect, &volume_expect);
  } else {
    double shape_expect[] = { 0.5*expectancy(pdf_size),
                              0.5*expectancy(pdf_size),
                              0.5*expectancy(pdf_size) };
    double r3 = 0.125*cubic_expectancy(pdf_size);
    double rad;
    MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape_expect, blockiness_expect, &rad);
    MathExtraLiggghtsNonspherical::volume_superquadric(shape_expect, blockiness_expect, &volume_expect);
    volume_expect *= r3 / (shape_expect[0]*shape_expect[1]*shape_expect[2]);
  }
  mass_expect = expectancy(pdf_density) * volume_expect;

}

/* ---------------------------------------------------------------------- */

FixTemplateSuperquadric::~FixTemplateSuperquadric()
{
    delete pdf_shapex;
    delete pdf_shapey;
    delete pdf_shapez;
    delete pdf_blockiness1;
    delete pdf_blockiness2;

    if(strcmp(style,"particletemplate/superquadric") == 0)
    {
        if(pti_list) delete_ptilist();
    }
}

/* ----------------------------------------------------------------------*/

void FixTemplateSuperquadric::init_ptilist(int n_random_max)///, const bool enforce_single, FixPropertyAtom * const fix_release)
{
    if(pti_list) error->one(FLERR,"invalid FixTemplateSuperquadric::init_list()");
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");
    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsertSuperquadric(lmp);
}

/* ----------------------------------------------------------------------*/

void FixTemplateSuperquadric::randomize_ptilist(int n_random,int distribution_groupbit)///,int distorder)
{
    for(int i = 0; i < n_random; i++)
    {
        ParticleToInsertSuperquadric *ptisq_ptr = dynamic_cast<ParticleToInsertSuperquadric*>(pti_list[i]);
        // if(!ptisq_ptr) --> error

        pti_list[i]->atom_type = atom_type;

        // randomize shape
        double shape[3];
        if(pdf_size == NULL) {
          shape[0] = rand(pdf_shapex, random);
          shape[1] = rand(pdf_shapey, random);
          shape[2] = rand(pdf_shapez, random);
        } else {
          double particle_size = 0.5*rand(pdf_size, random);
          shape[0] = particle_size;
          shape[1] = particle_size;
          shape[2] = particle_size;
        }

        ptisq_ptr->shape_ins[0] = shape[0];
        ptisq_ptr->shape_ins[1] = shape[1];
        ptisq_ptr->shape_ins[2] = shape[2];
        double blockiness[] = { rand(pdf_blockiness1, random),
                               rand(pdf_blockiness2, random) };
        ptisq_ptr->blockiness_ins[0] = blockiness[0];
        ptisq_ptr->blockiness_ins[1] = blockiness[1];
        double radius_;
        MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape, blockiness, &radius_);
        pti_list[i]->radius_ins[0] = pti_list[i]->r_bound_ins = radius_;

        // randomize density
        pti_list[i]->density_ins = rand(pdf_density,random);

        // calculate volume, mass and main components of inertia tensor
        MathExtraLiggghtsNonspherical::volume_superquadric(ptisq_ptr->shape_ins, ptisq_ptr->blockiness_ins, &(ptisq_ptr->volume_ins));
        MathExtraLiggghtsNonspherical::area_superquadric(ptisq_ptr->shape_ins, ptisq_ptr->blockiness_ins, &(ptisq_ptr->area_ins));
        pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;
        MathExtraLiggghtsNonspherical::inertia_superquadric(ptisq_ptr->shape_ins, ptisq_ptr->blockiness_ins, ptisq_ptr->density_ins, ptisq_ptr->inertia_ins);

        // init insertion position
        vectorZeroize3D(pti_list[i]->x_ins[0]);
        vectorZeroize3D(pti_list[i]->v_ins);
        vectorZeroize3D(pti_list[i]->omega_ins);

        pti_list[i]->groupbit = groupbit | distribution_groupbit;
    }

}

/* ----------------------------------------------------------------------*/

double FixTemplateSuperquadric::min_rad() const
{
    double shape_min[3];
    if(pdf_size == NULL) {
      shape_min[0] = pdf_min(pdf_shapex);
      shape_min[1] = pdf_min(pdf_shapey);
      shape_min[2] = pdf_min(pdf_shapez);
    } else {
      shape_min[0] = 0.5*pdf_min(pdf_size);
      shape_min[1] = 0.5*pdf_min(pdf_size);
      shape_min[2] = 0.5*pdf_min(pdf_size);
    }
    double blockiness_min[] = { pdf_min(pdf_blockiness1),
                               pdf_min(pdf_blockiness2) };
    double rad;
    MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape_min, blockiness_min, &rad);
    return rad;
}

/* ----------------------------------------------------------------------*/

double FixTemplateSuperquadric::max_rad() const
{
    double shape_max[3];
    if(pdf_size == NULL) {
      shape_max[0] = pdf_max(pdf_shapex);
      shape_max[1] = pdf_max(pdf_shapey);
      shape_max[2] = pdf_max(pdf_shapez);
    } else {
      shape_max[0] = 0.5*pdf_max(pdf_size);
      shape_max[1] = 0.5*pdf_max(pdf_size);
      shape_max[2] = 0.5*pdf_max(pdf_size);
    }
    double blockiness_max[] = { pdf_max(pdf_blockiness1),
                               pdf_max(pdf_blockiness2) };
    double rad;
    MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(shape_max, blockiness_max, &rad);
    return rad;
}

/* ----------------------------------------------------------------------*/

double FixTemplateSuperquadric::max_r_bound() const
{
  return max_rad();
}

#endif
