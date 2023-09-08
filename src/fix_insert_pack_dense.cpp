/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2017-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
   Daniel Queteschiner (JKU Linz)
   ---------------------------------------------------------------------- */


#include "fix_insert_pack_dense.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits>

#include "force.h"
#include "region.h"
#include "modify.h"
#include "error.h"
#include "domain.h"
#include "fix_particledistribution_discrete.h"
#include "fix_property_atom.h"
#include "input.h"
#include "random_park.h"
#include "update.h"
#include "variable.h"
#include "particleToInsert.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "vector_liggghts.h"
#include "region_neighbor_list.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace LIGGGHTS;
using namespace FixConst;

#define SMALL 1e-10

/* ---------------------------------------------------------------------- */

// this is the maximum volume fraction that the insertion algorithm
// can achieve - usually, it is lower.
const double FixInsertPackDense::max_volfrac = 0.57;

FixInsertPackDense::FixInsertPackDense(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  x_init(NULL),
  ins_region(NULL),
  idregion(NULL),
  var(NULL),
  fix_distribution(NULL),
  target_volfrac(max_volfrac),
  random(NULL),
  seed(-1),
  neighlist(lmp),
  insertion_done(false),
  is_inserter(true),
  insert_every(0),
  n_inserted(0),
  n_inserted_local(0),
  property_name(NULL),
  fix_property(NULL),
  fix_property_value(0.),
  fix_property_ivalue(0),
  property_index(-1),
  property_iindex(-1)
{
  int iarg = 3;

  if(strcmp(arg[iarg++],"seed")) error->fix_error(FLERR,this,"expecting keyword 'seed'");
  seed = atoi(arg[iarg++]) + comm->me;
  if (seed <= 0) error->fix_error(FLERR,this,"illegal seed");

  // random number generator, seed depends on proc
  random = new RanPark(lmp,seed);

  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"distributiontemplate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0 || strcmp(modify->fix[ifix]->style,"particledistribution/discrete"))
        error->fix_error(FLERR,this,"Fix insert requires you to define a valid ID for a fix of type particledistribution/discrete");
      fix_distribution = static_cast<FixParticledistributionDiscrete*>(modify->fix[ifix]);
      if (fix_distribution->has_multisphere())
        error->fix_error(FLERR,this,"no multisphere templates allowed for this insertion fix");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      ins_region = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg],"volumefraction_region") == 0) {
      if (iarg+2>narg) error->fix_error(FLERR,this,"");
      target_volfrac = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"insert_every") == 0 || strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      insert_every = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg],"target_variable") == 0) {
      if (iarg+4>narg) error->fix_error(FLERR,this,"");
      int n = strlen(arg[iarg+1]) + 1;
      var = new char[n];
      strcpy(var,arg[iarg+1]);
      var_insvalid = atof(arg[iarg+2]);
      var_threshold = atof(arg[iarg+3]);
      iarg += 4;
      hasargs = true;
    } else if (strcmp(arg[iarg],"set_property") == 0) {
      if (iarg+3 > narg) error->fix_error(FLERR,this,"");
      int n = strlen(arg[iarg+1]) + 1;
      property_name = new char[n];
      strcpy(property_name,arg[iarg+1]);
      fix_property_value = force->numeric(FLERR,arg[iarg+2]);
      fix_property_ivalue = static_cast<int>(fix_property_value);
      iarg += 3;
      hasargs = true;
    }
  }

  if (target_volfrac > max_volfrac) {
    char errmsg[500];
    sprintf(errmsg,"target_volfrac higher than %f - might not be achieved. Reseting to maximum value.",max_volfrac);
    error->warning(FLERR,errmsg);
    target_volfrac = max_volfrac;
  }

  radius_factor = cbrt(max_volfrac/target_volfrac);
  if (comm->me == 0 && screen)
    fprintf(screen, "radius scaling factor: %f\n",radius_factor);

  if (!ins_region)
    error->fix_error(FLERR,this,"no insertion region provided");

  if (!fix_distribution)
    error->fix_error(FLERR,this,"no particle distribution provided");

  // force reneighboring in next timestep
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep+1;
  most_recent_ins_step = -1;

  maxrad = fix_distribution->max_rad()*radius_factor;

  // set bounding box
  ins_bbox = BoundingBox(ins_region->extent_xlo,ins_region->extent_xhi,
                         ins_region->extent_ylo,ins_region->extent_yhi,
                         ins_region->extent_zlo,ins_region->extent_zhi);

}

/* ---------------------------------------------------------------------- */

int FixInsertPackDense::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::post_create()
{

  delete[] x_init;
  x_init = new double[3];

  double sublo[3],subhi[3];
  double const insreg_cutoff = maxrad;
  for (int i=0; i<3; ++i) {
    sublo[i] = domain->sublo[i]+insreg_cutoff;
    subhi[i] = domain->subhi[i]-insreg_cutoff;
  }

  ins_bbox.shrinkToSubbox(sublo,subhi);

  // for insertion of inital config we need an insertion box with larger cutoff:
  // circumscribed sphere (of 3 spheres with maxrad radius placed on the corners
  // of an equilateral triangle) has a radius of (1. + 2./sqrt(3.))*maxrad ~ 2.155*maxrad
  // where 2./sqrt(3.)*maxrad is the circumscribed circle of the triangle
  double const init_cutoff = 2.155*maxrad;
  for (int i=0; i<3; ++i) {
    sublo[i] = domain->sublo[i]+init_cutoff;
    subhi[i] = domain->subhi[i]-init_cutoff;
  }

  BoundingBox init_config_bbox(ins_region->extent_xlo+init_cutoff,ins_region->extent_xhi-init_cutoff,
                               ins_region->extent_ylo+init_cutoff,ins_region->extent_yhi-init_cutoff,
                               ins_region->extent_zlo+init_cutoff,ins_region->extent_zhi-init_cutoff);
  init_config_bbox.shrinkToSubbox(sublo,subhi);
  init_config_bbox.getCenter(x_init);

  if (!init_config_bbox.hasVolume()) {
    // no insertion volume on this subdomain
    is_inserter = false;
  } else if (!init_config_bbox.isInside(x_init) ||
             !ins_region->match_shrinkby_cut(x_init, init_cutoff)) {
    bool const subdomain_flag(true);
    ins_region->generate_random_shrinkby_cut(x_init,init_cutoff,subdomain_flag);
  }

  neighlist.reset();
  neighlist.setBoundingBox(ins_bbox,fix_distribution->max_rad()*radius_factor);
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::init()
{
    if(property_name)
    {
        fix_property = static_cast<FixPropertyAtom*>(modify->find_fix_property(property_name,"property/atom","scalar",1,1,this->style,false));
        if(!fix_property)
        {
            if(strstr(property_name,"i_") == property_name)
            {
                int flag;
                property_iindex = atom->find_custom(&property_name[2],flag);
                if(property_iindex < 0 || flag != 0)
                {
                    char errmsg[500];
                    sprintf(errmsg,"Could not locate a property storing value(s) for %s as requested by %s.",property_name,this->style);
                    error->all(FLERR,errmsg);
                }
            }
            else if(strstr(property_name,"d_") == property_name)
            {
                int flag;
                property_index = atom->find_custom(&property_name[2],flag);
                if(property_index < 0 || flag != 1)
                {
                    char errmsg[500];
                    sprintf(errmsg,"Could not locate a property storing value(s) for %s as requested by %s.",property_name,this->style);
                    error->all(FLERR,errmsg);
                }
            }
            else
            {
                char errmsg[500];
                sprintf(errmsg,"Could not locate a property storing value(s) for %s as requested by %s.",property_name,this->style);
                error->all(FLERR,errmsg);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::pre_exchange()
{
  // various checks if insertion should take place
  // 1. time to insert?
  // 2. periodic insertion?
  // 3. insertion according to external variable

  if (next_reneighbor != update->ntimestep || most_recent_ins_step == update->ntimestep) return;
  most_recent_ins_step = update->ntimestep;

  if (insert_every > 0) next_reneighbor += insert_every;
  else if (insertion_done) return;

  if (var)
  {
    int ivar = input->variable->find(var);
    if (ivar < 0) error->fix_error(FLERR,this,"target variable not found");
    double value = input->variable->compute_equal(ivar);
    if (fabs(value-var_insvalid) > var_threshold) return;
  }

  insertion_done = true;

  n_inserted = 0;
  n_inserted_local = 0;

  bool prepared = prepare_insertion();

  if (screen) {
    // only proc 0 writes general output
    if (comm->me == 0) {
      fprintf(screen,"fix insert/pack/dense %s will attempt to insert approximately %d particles in region %s\n",
              id,n_insert_estim,idregion);
      fprintf(screen,"maximum volume fraction: %f\n",max_volfrac);
      fprintf(screen,"target volume fraction:  %f\n",target_volfrac);
    }

    fprintf(screen,"process %d : attempting to insert approximately %d particles\n",
         comm->me,n_insert_estim_local);
  }

  // insert first three particles to initialize the algorithm
  if (prepared && is_inserter && insert_first_particles()) {

    const int n_write = n_insert_estim_local/10;
    int n_next_write = n_write;

    while(!frontSpheres.empty()) {
      if (screen && !(fix_distribution->pti_list.size() == static_cast<FixParticledistributionDiscrete::pti_list_type::size_type>(n_inserted_local))) {
        fprintf(screen, "pti_list.size() %lu | n_inserted_local %d\n",fix_distribution->pti_list.size(),n_inserted_local);
      }
      assert(fix_distribution->pti_list.size() == static_cast<FixParticledistributionDiscrete::pti_list_type::size_type>(n_inserted_local));
      handle_next_front_sphere();

      if (n_inserted_local > n_next_write) {
        double percent = static_cast<double>(n_inserted_local)/static_cast<double>(n_insert_estim_local)*100;
        if (screen) fprintf(screen,"process %d : %2.0f%% done, inserted %d/%d particles\n",
                            comm->me,percent,n_inserted_local,n_insert_estim_local);
        n_next_write += n_write;
      }
    }
  }

  // actual insertion
  fix_distribution->pre_insert(n_inserted_local,fix_property,fix_property_value,property_index,fix_property_ivalue,property_iindex);

  fix_distribution->insert(n_inserted_local);

  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  fix_distribution->finalize_insertion(true);

  n_inserted = n_inserted_local;
  MPI_Sum_Scalar(n_inserted,world);

  if (comm->me == 0 && screen) {
    double percent = static_cast<double>(n_inserted)/static_cast<double>(n_insert_estim)*100.;
    fprintf(screen,"inserted %d of %d particles (%4.2f%%)\n",n_inserted,n_insert_estim,percent);
    double const volfrac_inserted
        = target_volfrac*static_cast<double>(n_inserted)/static_cast<double>(n_insert_estim);
    fprintf(screen,"volume fracton: %f\n",volfrac_inserted);
  }
}

/* ---------------------------------------------------------------------- */

bool FixInsertPackDense::prepare_insertion()
{
  // need to iterate twice: first time to find largest radius, second
  // time to compute volume and actually insert already present
  // particles into region neighbor list with now appropriate bin size
  double volume_present_local = 0.;

  neighlist.reset();
  neighlist.setBoundingBox(ins_bbox,fix_distribution->max_rad()*radius_factor);
  distfield.reset();

  if (is_inserter) {
    double rad_max_present = 0.;
    for(int i=0;i<atom->nlocal;i++){
      if (ins_region->match_expandby_cut(atom->x[i],atom->radius[i])) {
        if(atom->radius[i] > rad_max_present) rad_max_present = atom->radius[i];
      }
    }

    if (rad_max_present > fix_distribution->max_rad()) {
      neighlist.reset();
      neighlist.setBoundingBox(ins_bbox,rad_max_present*radius_factor);
    }

    for (int i=0;i<atom->nlocal;i++) {
      if (ins_region->match_expandby_cut(atom->x[i],atom->radius[i])) {
        neighlist.insert(atom->x[i],atom->radius[i],atom->type[i]);
        volume_present_local += MathConst::MY_4PI3*atom->radius[i]*atom->radius[i]*atom->radius[i];
      }
    }
  }

  double volume_present = volume_present_local;
  MPI_Sum_Scalar(volume_present,world);

  // estimate # of particles to insert
  int const ntry_mc = 100000;
  bool const cutflag(false);
  ins_region->volume_mc(ntry_mc,cutflag,fix_distribution->max_r_bound(),
                          region_volume,region_volume_local);

  double const vol_part_ave = fix_distribution->vol_expect();
  n_insert_estim = floor((region_volume-volume_present)*target_volfrac/vol_part_ave);

  if (!is_inserter) {
    n_insert_estim_local = 0;
    return false;
  }

  n_insert_estim_local = floor((region_volume_local-volume_present_local)*target_volfrac/vol_part_ave);

  // calculate distance field
  distfield.build(ins_region,ins_bbox,fix_distribution->max_rad()*radius_factor);

  return true;
}

/* ---------------------------------------------------------------------- */

bool FixInsertPackDense::insert_first_particles()
{
  ParticleToInsert *pti1 = fix_distribution->get_random_particle(groupbit);
  ParticleToInsert *pti2 = fix_distribution->get_random_particle(groupbit);
  ParticleToInsert *pti3 = fix_distribution->get_random_particle(groupbit);

  // check if starting point does not conflict with any pre-existing particles
  double const maxrad_init = 2.155*maxrad;
  if(neighlist.hasOverlap(x_init,maxrad_init,pti1->get_atom_type()) ||
     ( pti1->get_atom_type() != pti2->get_atom_type() && neighlist.hasOverlap(x_init,maxrad_init,pti2->get_atom_type()) ) ||
     ( pti1->get_atom_type() != pti3->get_atom_type() && pti2->get_atom_type() != pti3->get_atom_type() && neighlist.hasOverlap(x_init,maxrad_init,pti3->get_atom_type()) )) {
    int const nAttemptMax = 100;
    int nAttempt = 0;
    for (; nAttempt<nAttemptMax; ++nAttempt) {
      ins_region->generate_random_shrinkby_cut(x_init,maxrad_init,true);
      if(!neighlist.hasOverlap(x_init,maxrad_init,pti1->get_atom_type()) &&
         ( pti1->get_atom_type() == pti2->get_atom_type() || !neighlist.hasOverlap(x_init,maxrad_init,pti2->get_atom_type()) ) &&
         ( pti1->get_atom_type() == pti3->get_atom_type() || pti2->get_atom_type() == pti3->get_atom_type() || !neighlist.hasOverlap(x_init,maxrad_init,pti3->get_atom_type()) ))
        break;
    }
    if (nAttempt == nAttemptMax) {
      char errmsg[500] = {0};
      snprintf(errmsg, 499,"fix '%s' failed to find suitable starting point for insertion on process %d",id,comm->me);
      error->warning(FLERR,errmsg);
      return false;
    }
  }

  generate_initial_config(pti1,pti2,pti3);

  fix_distribution->pti_list.push_back(pti1);
  fix_distribution->pti_list.push_back(pti2);
  fix_distribution->pti_list.push_back(pti3);

  Particle p1 = particle_from_pti(pti1);
  Particle p2 = particle_from_pti(pti2);
  Particle p3 = particle_from_pti(pti3);

  p1.radius *= radius_factor;
  p2.radius *= radius_factor;
  p3.radius *= radius_factor;

  neighlist.insert(p1.x,p1.radius,p1.type);
  neighlist.insert(p2.x,p2.radius,p2.type);
  neighlist.insert(p3.x,p3.radius,p3.type);

  frontSpheres.push(p1);
  frontSpheres.push(p2);
  frontSpheres.push(p3);

  n_inserted_local += 3;

  return true;
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::handle_next_front_sphere()
{
  Particle current = frontSpheres.front();
  RegionNeighborList::ParticleBin particles;

  do {
    ParticleToInsert *pti = get_next_pti();
    double const r_insert = pti->radius_ins[0]*radius_factor;
    double const type_insert = pti->get_atom_type();
    double const cutoff_dist = current.radius+2.*r_insert;

    particles.clear();
    neighlist.getParticlesCloseTo(current.x,cutoff_dist,particles);

    // identify candidate points
    candidatePoints.clear();
    for (unsigned int i=0; i<particles.size()-1; ++i) {
      for (unsigned int j=i+1; j<particles.size(); ++j) {
        compute_and_append_candidate_points(current,particles[i],particles[j],r_insert,type_insert);
      }
    }

    if(candidatePoints.empty()){
      rejectedSpheres.push(pti);
      break;
    }

    // then, search for candidate point closest to insertion center
    double dist_min_sq = std::numeric_limits<double>::max();
    ParticleVector::iterator closest_candidate;
    for(ParticleVector::iterator it = candidatePoints.begin(); it != candidatePoints.end(); ++it){
      double dist_sq = pointDistanceSquared((*it).x,x_init);
      if(dist_sq < dist_min_sq){
        dist_min_sq = dist_sq;
        closest_candidate = it;
      }
    }

    vectorCopy3D((*closest_candidate).x,pti->x_ins[0]);
    fix_distribution->pti_list.push_back(pti);
    frontSpheres.push(*closest_candidate);
    neighlist.insert((*closest_candidate).x,(*closest_candidate).radius,(*closest_candidate).type);
    n_inserted_local++;

  } while(candidatePoints.size() > 1);

  frontSpheres.pop();
}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::generate_initial_config(ParticleToInsert *&p1,
                                                 ParticleToInsert *&p2,
                                                 ParticleToInsert *&p3)
{
  double x1[3],x2[3],x3[3];
  vectorZeroize3D(x1);
  vectorZeroize3D(x2);
  vectorZeroize3D(x3);

  // first, construct touching spheres
  double const r1 = p1->radius_ins[0]*radius_factor;
  double const r2 = p2->radius_ins[0]*radius_factor;
  double const r3 = p3->radius_ins[0]*radius_factor;

  x2[0] = r1+r2;

  double const a = r2+r3;
  double const b = r1+r3;
  double const c = r1+r2;
  double const alpha = acos((a*a-b*b-c*c)/(-2.*b*c));

  x3[0] = b*cos(alpha);
  x3[1] = b*sin(alpha);

  // rotate by random quaternion
  double rotquat[4];
  MathExtraLiggghts::random_unit_quat(random,rotquat);
  MathExtraLiggghts::vec_quat_rotate(x2,rotquat);
  MathExtraLiggghts::vec_quat_rotate(x3,rotquat);

  // then, compute COM & move COM to origin
  double com[3];
  for(int i=0;i<3;i++)
    com[i] = (x1[i]+x2[i]+x3[i])/3.;

  MathExtra::sub3(x1,com,x1);
  MathExtra::sub3(x2,com,x2);
  MathExtra::sub3(x3,com,x3);

  // then, move to starting point & write to PTI
  MathExtra::add3(x1,x_init,p1->x_ins[0]);
  MathExtra::add3(x2,x_init,p2->x_ins[0]);
  MathExtra::add3(x3,x_init,p3->x_ins[0]);

}

/* ---------------------------------------------------------------------- */

void FixInsertPackDense::compute_and_append_candidate_points(Particle const &p1,
                                                             Particle const &p2,
                                                             Particle const &p3,
                                                             double const r_insert,
                                                             double const type_insert)
{
  double const halo1 = p1.radius+r_insert+SMALL;
  double const halo2 = p2.radius+r_insert+SMALL;
  double const halo3 = p3.radius+r_insert+SMALL;

  // exclude impossible combinations
  double const dist_12_sq = pointDistanceSquared(p1.x,p2.x);
  if(dist_12_sq > (halo1+halo2)*(halo1+halo2) || dist_12_sq < SMALL*SMALL)
    return;
  double const dist_13_sq = pointDistanceSquared(p1.x,p3.x);
  if(dist_13_sq > (halo1+halo3)*(halo1+halo3) || dist_13_sq < SMALL*SMALL)
    return;
  double const dist_23_sq = pointDistanceSquared(p2.x,p3.x);
  if(dist_23_sq > (halo2+halo3)*(halo2+halo3) || dist_23_sq < SMALL*SMALL)
    return;

  // first, construct intersection circle between halos of particles 1,2
  double const alpha = 0.5*(1. - (halo2*halo2-halo1*halo1)/dist_12_sq );
  if(alpha < 0. || alpha > 1.) return;

  // center and radius of intersection circle
  double c_c[3];
  for(int i=0;i<3;i++) c_c[i] = (1.-alpha)*p1.x[i] + alpha*p2.x[i];
  double const dist_x1_cc_sq = pointDistanceSquared(p1.x,c_c);
  double const r_c_sq = halo1*halo1 - dist_x1_cc_sq;
  double const r_c = sqrt(r_c_sq);

  // normal vector from p1 to p2
  double n[3];
  MathExtra::sub3(p1.x,p2.x,n);
  MathExtra::normalize3(n,n);

  // normal distance of p3 to intersection plane
  double tmp_vec[3];
  MathExtra::sub3(p3.x,c_c,tmp_vec);
  double const lambda = MathExtra::dot3(tmp_vec,n);
  if(lambda > halo3 || -lambda > halo3) return;


  // intersection circle of circle plane and third sphere
  vectorCopy3D(n,tmp_vec);
  MathExtra::scale3(lambda,tmp_vec);
  double c_p[3];
  MathExtra::sub3(p3.x,tmp_vec,c_p);
  double const r_p_sq = halo3*halo3-lambda*lambda;
  double const r_p = sqrt(r_p_sq);

  double const dist_pc_sq = pointDistanceSquared(c_c,c_p);
  if(dist_pc_sq > (r_p+r_c)*(r_p+r_c)) return;

  double const alpha2 = 0.5*(1. - (r_p_sq-r_c_sq)/dist_pc_sq);
  double c_m[3];
  for(int i=0;i<3;i++) c_m[i] = (1.-alpha2)*c_c[i] + alpha2*c_p[i];

  double const dist_cc_cm_sq = pointDistanceSquared(c_c,c_m);
  if(dist_cc_cm_sq > r_c_sq) return;

  double const h = sqrt(r_c_sq - dist_cc_cm_sq);


  if(h < SMALL){ // only one candidate point
    Particle candidate(c_m,r_insert,type_insert);

    if(candidate_point_is_valid(candidate)){
      candidatePoints.push_back(candidate);
    }
  } else { // two candidate points
    Particle candidate1(c_m,r_insert,type_insert);
    Particle candidate2(c_m,r_insert,type_insert);

    MathExtra::sub3(c_p,c_c,tmp_vec);
    MathExtra::normalize3(tmp_vec,tmp_vec);
    double vec_tmp2[3];
    MathExtra::cross3(n,tmp_vec,vec_tmp2);

    MathExtra::scale3(h,vec_tmp2);

    MathExtra::add3(candidate1.x,vec_tmp2,candidate1.x);
    MathExtra::sub3(candidate2.x,vec_tmp2,candidate2.x);

    if(candidate_point_is_valid(candidate1)){
      candidatePoints.push_back(candidate1);
    }
    if(candidate_point_is_valid(candidate2)){
      candidatePoints.push_back(candidate2);
    }
  }
}

/* ---------------------------------------------------------------------- */

ParticleToInsert* FixInsertPackDense::get_next_pti()
{
  if(rejectedSpheres.empty())
    return fix_distribution->get_random_particle(groupbit);

  // this might be not all too efficient: if a particle gets rejected
  // multiple times, it also gets inserted and deleted from the list
  // each time. Will refactor if performance critical.
  ParticleToInsert *pti = rejectedSpheres.front();
  rejectedSpheres.pop();
  return pti;
}

/* ---------------------------------------------------------------------- */

FixInsertPackDense::~FixInsertPackDense()
{
  delete[] x_init;
  delete[] idregion;
  delete[] property_name;
  delete random;
}

/* ---------------------------------------------------------------------- */

Particle FixInsertPackDense::particle_from_pti(ParticleToInsert* pti)
{
  Particle p(pti->x_ins[0],pti->radius_ins[0],pti->get_atom_type());
  return p;
}

/* ---------------------------------------------------------------------- */

bool FixInsertPackDense::is_completely_in_subregion(Particle &p)
{
  if(distfield.isInside(p.x))
    return true;
  if(distfield.isOutside(p.x))
    return false;

  // position was either not found or is in boundary cell of distfield
  return ins_bbox.isInside(p.x) && ins_region->match_shrinkby_cut(p.x,p.radius);
}

/* ---------------------------------------------------------------------- */

bool FixInsertPackDense::candidate_point_is_valid(Particle &p)
{
  return ( !neighlist.hasOverlap(p.x,p.radius,p.type) && is_completely_in_subregion(p) );
}
