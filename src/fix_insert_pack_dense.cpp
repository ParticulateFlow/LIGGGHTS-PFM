/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   contributing authors
   Philippe Seil (JKU)
   ---------------------------------------------------------------------- */


#include "fix_insert_pack_dense.h"

#include <cstdlib>

#include "region.h"
#include "modify.h"
#include "error.h"
#include "domain.h"
#include "fix_particledistribution_discrete.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixInsertPackDense::FixInsertPackDense(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  x_init(0),
  ins_region(0),
  idregion(0),
  fix_distribution(0),
  random(0),
  seed(-1)
{
  int iarg = 3;

  if(strcmp(arg[iarg++],"seed")) error->fix_error(FLERR,this,"expecting keyword 'seed'");
  seed = atoi(arg[iarg++]) + comm->me;
  if (seed <= 0) error->fix_error(FLERR,this,"illegal seed");

  // random number generator, seed depends on proc
  random = new RanPark(lmp,seed);

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"distributiontemplate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int ifix = modify->find_fix(arg[iarg+1]);
      if(ifix < 0 || strcmp(modify->fix[ifix]->style,"particledistribution/discrete"))
        error->fix_error(FLERR,this,"Fix insert requires you to define a valid ID for a fix of type particledistribution/discrete");
      fix_distribution = static_cast<FixParticledistributionDiscrete*>(modify->fix[ifix]);
      iarg += 2;
      hasargs = true;
    } else  if(strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      ins_region = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    }
  }

  if(!ins_region)
    error->fix_error(FLERR,this,"no insertion region provided");
  if(!fix_distribution)
    error->fix_error(FLERR,this,"no particle distribution provided");

  maxrad = fix_distribution->max_rad();

}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::post_create()
{
  x_init = new double[3];
  
  // TODO: account for extents of local subdomain.
  // this only works on single proc
  x_init[0] = 0.5*(ins_region->extent_xlo + ins_region->extent_xhi);
  x_init[1] = 0.5*(ins_region->extent_ylo + ins_region->extent_yhi);
  x_init[2] = 0.5*(ins_region->extent_zlo + ins_region->extent_zhi);

  if(!ins_region->inside(x_init[0],x_init[1],x_init[2])){
    // TODO: get a random point, change error to warning
    error->fix_error(FLERR,this,"starting point not in insertion region");
  }
}

FixInsertPackDense::~FixInsertPackDense()
{
  if(x_init) delete[] x_init;
}

/* ---------------------------------------------------------------------- */

