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

#include <stdlib.h>
#include <math.h>

#include "region.h"
#include "modify.h"
#include "error.h"
#include "domain.h"
#include "fix_particledistribution_discrete.h"
#include "random_park.h"
#include "update.h"
#include "particleToInsert.h"
#include "math_extra.h"
#include "region_neighbor_list.h"
#include "bounding_box.h"

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
      if(fix_distribution->has_multisphere())
        error->fix_error(FLERR,this,"no multisphere templates allowed for this insertion fix");
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

  // force reneighboring in next timestep
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep+1;

  maxrad = fix_distribution->max_rad();

}

/* ---------------------------------------------------------------------- */

int FixInsertPackDense::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  printf("setmask of fix_insert_pack_dense: mask = %d\n",mask);
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::post_create()
{
  printf("post_create of fix_insert_pack_dense\n");

  if(x_init) delete[] x_init;
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

  // TODO: change this to account for parallel stuff
  BoundingBox box(ins_region->extent_xlo,ins_region->extent_xhi,
                  ins_region->extent_ylo,ins_region->extent_yhi,
                  ins_region->extent_zlo,ins_region->extent_zhi);
  neighlist.setBoundingBox(box,fix_distribution->max_rad());
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::pre_exchange()
{
  // this fix should only run exactly once
  static bool insertion_done(false);
  if(insertion_done) return;
  insertion_done = true;

  // TODO: check if insertion region is empty
  
  ParticleToInsert *pti1 = fix_distribution->get_random_particle(groupbit);
  ParticleToInsert *pti2 = fix_distribution->get_random_particle(groupbit);
  ParticleToInsert *pti3 = fix_distribution->get_random_particle(groupbit);

  generate_initial_config(pti1,pti2,pti3);

  fix_distribution->pti_list.push_back(pti1);
  fix_distribution->pti_list.push_back(pti2);
  fix_distribution->pti_list.push_back(pti3);

  Particle p1 = particle_from_pti(pti1);
  Particle p2 = particle_from_pti(pti2);
  Particle p3 = particle_from_pti(pti3);

  neighlist.insert(p1);
  neighlist.insert(p2);
  neighlist.insert(p3);
  
  int nInserted = 3;

  printf("perform actual insertion\n");
  // actual insertion
  fix_distribution->pre_insert();
  printf("insert\n");
  fix_distribution->insert(nInserted);
  printf("finalize_insertion\n");

  if (atom->tag_enable)
  {
    atom->tag_extend();
    if (atom->map_style)
    {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  fix_distribution->finalize_insertion();
  printf("done\n");
  
}

void FixInsertPackDense::generate_initial_config(ParticleToInsert *&p1,
                                                 ParticleToInsert *&p2,
                                                 ParticleToInsert *&p3)
{

  printf("generate_initial_config\n");
  
  double x1[3],x2[3],x3[3];
  vectorZeroize3D(x1);
  vectorZeroize3D(x2);
  vectorZeroize3D(x3);

  // first, construct touching spheres
  double const r1=p1->radius_ins[0],r2=p2->radius_ins[0],r3=p3->radius_ins[0];
  x2[0] = r1+r2;

  double const a=r2+r3,b=r1+r3,c=r1+r2;
  double const alpha = acos((a*a-b*b-c*c)/(-2.*b*c));

  x3[0] = b*cos(alpha);
  x3[1] = b*sin(alpha);

  // then, compute COM & move COM to origin
  double com[3];
  for(int i=0;i<3;i++)
    com[i] = (x1[i]+x2[i]+x3[i])/3.;

  MathExtra::sub3(x1,com,x1);
  MathExtra::sub3(x2,com,x2);
  MathExtra::sub3(x3,com,x3);
  
  // maybe, at some point in the future, implement random orientation
  // of initial packing
  
  // then, move to starting point & write to PTI
  MathExtra::add3(x1,x_init,p1->x_ins[0]);
  MathExtra::add3(x2,x_init,p2->x_ins[0]);
  MathExtra::add3(x3,x_init,p3->x_ins[0]);

}

/* ---------------------------------------------------------------------- */

FixInsertPackDense::~FixInsertPackDense()
{
  if(x_init) delete[] x_init;
}

/* ---------------------------------------------------------------------- */

Particle FixInsertPackDense::particle_from_pti(ParticleToInsert* pti)
{
  Particle p(pti->x_ins[0],pti->radius_ins[0]);
  return p;
}
