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
#include "vector_liggghts.h"
#include "region_neighbor_list.h"
#include "bounding_box.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-10

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

  // insert first three particles to initialize the algorithm
  insert_first_particles();

  // insert_next_particle() returns false
  // if algorithm terminates
  while(insert_next_particle() && neighlist.count() < 1001)
    {
      printf("neighlist count %d\n",neighlist.count()); 
    }

  // actual insertion
  fix_distribution->pre_insert();
  fix_distribution->insert(neighlist.count());

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
  
}

void FixInsertPackDense::insert_first_particles()
{
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
  
  frontSpheres.push_back(p1);
  frontSpheres.push_back(p2);
  frontSpheres.push_back(p3);
}

bool FixInsertPackDense::insert_next_particle()
{
  Particle frontSphere = frontSpheres.front();

  bool inserted(false);
  while(!inserted){
    int nValid = try_front_sphere(frontSphere);
    inserted = nValid > 0;
    if(nValid < 2){
      frontSpheres.pop_front();
      if(frontSpheres.empty())
        return false;
    }
    if(!inserted)
      frontSphere = frontSpheres.front();
  }
  
  return true;
}

// returns number of valid candidate positions
// if positions > 0, also performs insertion
int FixInsertPackDense::try_front_sphere(Particle &p)
{
  ParticleToInsert *pti = get_next_pti();
  double r_insert = pti->radius_ins[0];

  candidatePoints.clear();
  
  RegionNeighborList::ParticleBin *particles = neighlist.getParticlesCloseTo(p.x);
  // get all possible combinations of spheres
  for(RegionNeighborList::ParticleBin::iterator i=particles->begin();i!=particles->end();++i){
    for(RegionNeighborList::ParticleBin::iterator j=i+1;j!=particles->end();++j){
      printf("p1 x %f %f %f | r %f\n",p.x[0],p.x[1],p.x[2],p.radius);
      printf("p2 x %f %f %f | r %f\n",(*i).x[0],(*i).x[1],(*i).x[2],(*i).radius);
      printf("p3 x %f %f %f | r %f\n",(*j).x[0],(*j).x[1],(*j).x[2],(*j).radius);
      // then, for each combination, compute candidate points
      // abuse Particle::r to hold distance to x_init in candidatePoints
      compute_and_store_candidate_points(p,*i,*j,r_insert);
    }
  }
  
  // then, check candidate points. At the same time, store which one is closest
  Particle *p_close = 0;
  double d_min = 1000;
  for(ParticleList::iterator it = candidatePoints.begin(); it != candidatePoints.end(); ++it){
    if(neighlist.hasOverlap((*it).x,r_insert-SMALL)){
      printf("removing particle\n");
      it = candidatePoints.erase(it);
      if(it != candidatePoints.begin()) --it;
      continue;
    }
    // abused radius to store distance to insertion center
    if((*it).radius < d_min){
      d_min = (*it).radius;
      if(!p_close)
        p_close = new Particle((*it).x,r_insert);
      else
        vectorCopy3D((*it).x,p_close->x);
    }
  }

  // clean up: delete particle bin
  delete particles;

  if(!p_close){ // no valid candidate point found
    rejectedSpheres.push_back(pti);
    return 0;
  }
  else{
    vectorCopy3D(p_close->x,pti->x_ins[0]);
    fix_distribution->pti_list.push_back(pti);
    frontSpheres.push_back(*p_close);
    neighlist.insert(*p_close);
    printf("inserting particle at %f %f %f | radius %f\n",
           p_close->x[0],p_close->x[1],p_close->x[2],p_close->radius);
  }
  
  return candidatePoints.size();
  
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

void FixInsertPackDense::compute_and_store_candidate_points(Particle const &p1,
                                                            Particle const &p2,
                                                            Particle const &p3,
                                                            double const r_insert)
{
  double const halo1 = p1.radius+r_insert;
  double const halo2 = p2.radius+r_insert;
  double const halo3 = p3.radius+r_insert;

  // exclude impossible combinations
  double const d_12_sqr = pointDistanceSqr(p1.x,p2.x);
  printf("d_12_sqr %f\n",d_12_sqr);
  if(d_12_sqr > (halo1+halo2)*(halo1+halo2) || d_12_sqr < SMALL*SMALL)
    return;
  double const d_13_sqr = pointDistanceSqr(p1.x,p3.x);
  printf("d_13_sqr %f\n",d_13_sqr);
  if(d_13_sqr > (halo1+halo3)*(halo1+halo3) || d_13_sqr < SMALL*SMALL)
    return;
  double const d_23_sqr = pointDistanceSqr(p2.x,p3.x);
  printf("d_23_sqr %f\n",d_23_sqr);
  if(d_23_sqr > (halo2+halo3)*(halo2+halo3) || d_23_sqr < SMALL*SMALL)
    return;

  // first, construct intersection circle between halos of particles 1,2
  double const alpha = 0.5*(1. - (halo2*halo2-halo1*halo1)/d_12_sqr );
  if(alpha < 0. || alpha > 1.) return;

  printf("alpha %f\n",alpha);
  
  // center and radius of intersection circle
  double c_c[3];
  for(int i=0;i<3;i++) c_c[i] = (1.-alpha)*p1.x[i] + alpha*p2.x[i];
  double const d_x1_cc_sqr = pointDistanceSqr(p1.x,c_c);
  double const r_c = sqrt(halo1*halo1 - d_x1_cc_sqr);

  // normal vector from p1 to p2
  double n[3];
  MathExtra::sub3(p1.x,p2.x,n);
  MathExtra::normalize3(n,n);
  
  // normal distance of p3 to intersection plane
  double tmp_vec[3];
  MathExtra::sub3(p3.x,c_c,tmp_vec);
  double const lambda = MathExtra::dot3(tmp_vec,n);
  printf("lambda %f\n",lambda);
  if(lambda > halo3 || -lambda > halo3) return;

  
  // intersection circle of circle plane and third sphere
  vectorCopy3D(n,tmp_vec);
  MathExtra::scale3(lambda,tmp_vec);
  double c_p[3];
  MathExtra::sub3(p3.x,tmp_vec,c_p);
  double const r_p = sqrt(halo3*halo3-lambda*lambda);

  double const dist_pc_sqr = pointDistanceSqr(c_c,c_p);
  if(dist_pc_sqr > (r_p+r_c)*(r_p+r_c)) return;

  printf("dist_pc_sqr %f",dist_pc_sqr);
  
  double const alpha2 = 0.5*(1. - (r_p*r_p-r_c*r_c)/dist_pc_sqr);
  double c_m[3];
  for(int i=0;i<3;i++) c_m[i] = (1.-alpha2)*c_c[i] + alpha2*c_p[i];

  double const d_sqr_cc_cm = pointDistanceSqr(c_c,c_m);
  if(d_sqr_cc_cm > r_c*r_c) return;

  printf("d_sqr_cc_cm %f\n",d_sqr_cc_cm);
  
  double const h = sqrt(r_c*r_c - d_sqr_cc_cm);

  
  if(h < SMALL){ // only one candidate point
    Particle candidate(c_m,0.);
    candidate.radius = pointDistance(c_m,x_init);
    candidatePoints.push_back(candidate);
    printf("found single candidate point at %f %f %f | distance to center: %f\n",
           candidate.x[0],candidate.x[1],candidate.x[2],candidate.radius);
  } else { // two candidate points
    Particle candidate1(c_m,0.),candidate2(c_m,0.);

    MathExtra::sub3(c_p,c_c,tmp_vec);
    MathExtra::normalize3(tmp_vec,tmp_vec);
    double vec_tmp2[3];
    MathExtra::cross3(n,tmp_vec,vec_tmp2);

    MathExtra::scale3(h,vec_tmp2);
    
    MathExtra::add3(candidate1.x,vec_tmp2,candidate1.x);
    MathExtra::sub3(candidate2.x,vec_tmp2,candidate2.x);

    candidate1.radius = pointDistance(candidate1.x,x_init);
    candidate2.radius = pointDistance(candidate2.x,x_init);

    candidatePoints.push_back(candidate1);
    candidatePoints.push_back(candidate2);
    
    printf("found candidate point at %f %f %f | distance to center: %f\n",
           candidate1.x[0],candidate1.x[1],candidate1.x[2],candidate1.radius);
    printf("found candidate point at %f %f %f | distance to center: %f\n",
           candidate2.x[0],candidate2.x[1],candidate2.x[2],candidate2.radius);
    printf("total candidate points: %d\n",candidatePoints.size());
  }
}

ParticleToInsert* FixInsertPackDense::get_next_pti()
{
  if(rejectedSpheres.empty())
    return fix_distribution->get_random_particle(groupbit);

  // this might be not all too efficient: if a particle gets rejected
  // multiple times, it also gets inserted and deleted from the list
  // each time. Will refactor if performance critical.
  ParticleToInsert *pti = rejectedSpheres.front();
  rejectedSpheres.pop_front();
  return pti;
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
