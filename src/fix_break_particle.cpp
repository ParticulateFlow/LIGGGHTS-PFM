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
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_break_particle.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_fragments.h"
#include "particleToInsert_fragments.h"
#include "pair_gran.h"
#include "force.h"
#include "neigh_list.h"
#include <string>
#include <algorithm>
#include "math_const.h"

using namespace MathConst;
using namespace LAMMPS_NS;
using namespace FixConst;

enum{BC_ENERGY, BC_FORCE, BC_VON_MISES};

#define LMP_DEBUGMODE_FIX_BREAK_PARTICLE true
#define LMP_DEBUG_OUT_FIX_BREAK_PARTICLE screen

/* ---------------------------------------------------------------------- */

FixBreakParticle::FixBreakParticle(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while (iarg < narg && hasargs) {
    hasargs = false;
    if (strcmp(arg[iarg],"fMat") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing fMat value");
      fMat = atof(arg[iarg+1]);
      if(fMat <= 0. ) error->fix_error(FLERR,this,"'fMat' must be > 0");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"energy_threshold") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing energy_threshold value");
      threshold = atof(arg[iarg+1]);
      if (threshold < 0.) error->fix_error(FLERR,this,"'energy_threshold' must be >= 0");
      breakage_criterion = BC_ENERGY;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"force_threshold") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing force_threshold value");
      threshold = atof(arg[iarg+1]);
      if (threshold <= 0.) error->fix_error(FLERR,this,"'force_threshold' must be > 0");
      breakage_criterion = BC_FORCE;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"von_mises_stress") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing von_mises_stress value");
      threshold = atof(arg[iarg+1]);
      if (threshold <= 0.) error->fix_error(FLERR,this,"'von_mises_stress' must be > 0");
      breakage_criterion = BC_VON_MISES;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"min_radius") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing min_rad value");
      min_break_rad = atof(arg[iarg+1]);
      if (min_break_rad < 0.) error->fix_error(FLERR,this,"'min_break_rad' must be >= 0");
      iarg += 2;
      hasargs = true;
    } else {
      std::string unknown("unknown keyword ");
      unknown += arg[iarg];
      error->fix_error(FLERR,this,unknown.c_str());
    }
  }

  // no fixed insertion target (such as # particles) since insertion triggered by break-ups
  ninsert_exists = 0;

  // turn off overlap check and turn off start stats since both not required
  check_ol_flag = 0;
  print_stats_start_flag = 0;

  breakdata = NULL;
  fix_break = NULL;
  fix_breaker = NULL;
  fix_collision_factor = NULL;
  tag_max = 0;

  // execute end of step
  nevery = 1;

  n_break = 0;
  mass_break = 0.;

  /////if(maxrad >= 1.) error->fix_error(FLERR,this,"Particle distribution must be relative, max radius mus be < 1");
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::post_create()
{
  if (!fix_break) {
    char * breakvar_name = new char[7+strlen(id)];
    sprintf(breakvar_name,"break_%s",id);

    const char * fixarg[9];
    fixarg[0] = breakvar_name;
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = breakvar_name;
    fixarg[4] = "scalar"; //NP 1 scalar per particle to be registered
    fixarg[5] = "yes";    //NP restart yes
    fixarg[6] = "yes";    //NP communicate ghost yes
    fixarg[7] = "no";     //NP communicate rev no
    fixarg[8] = "0.";
    fix_break = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    delete [] breakvar_name;
  }
  if (!fix_breaker) {
    char * breakvar_name = new char[9+strlen(id)];
    sprintf(breakvar_name,"breaker_%s",id);

    const char * fixarg[9];
    fixarg[0] = breakvar_name;
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = breakvar_name;
    fixarg[4] = "scalar"; //NP 1 scalar per particle to be registered
    fixarg[5] = "yes";    //NP restart yes
    fixarg[6] = "yes";    //NP communicate ghost yes
    fixarg[7] = "no";     //NP communicate rev no
    fixarg[8] = "0.";
    fix_breaker = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    delete [] breakvar_name;
  }
  if (!fix_collision_factor) {
    char * breakvar_name = new char[10+strlen(id)];
    sprintf(breakvar_name,"break_cf_%s",id);

    const char * fixarg[9];
    fixarg[0] = breakvar_name;
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = breakvar_name;
    fixarg[4] = "scalar"; //NP 1 scalar per particle to be registered
    fixarg[5] = "yes";    //NP restart yes
    fixarg[6] = "yes";    //NP communicate ghost yes
    fixarg[7] = "no";     //NP communicate rev no
    fixarg[8] = "1.";
    fix_collision_factor = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    delete [] breakvar_name;
  }
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::pre_delete(bool unfixflag)
{
  modify->delete_fix(fix_break->id);
}

/* ---------------------------------------------------------------------- */

FixBreakParticle::~FixBreakParticle()
{
  if(breakdata) memory->destroy(breakdata);
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::init_defaults()
{
  threshold = -1.0;
  fMat = 1.0;
  breakage_criterion = BC_ENERGY;
  fix_fragments = NULL;
  min_break_rad = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::init()
{
  FixInsert::init();
  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

  dnum = pair_gran->dnum();
  deltaMaxOffset = pair_gran->get_history_value_offset("deltaMax");
  siblingOffset = pair_gran->get_history_value_offset("sibling");
  collisionFactorOffset = pair_gran->get_history_value_offset("collisionFactor");
  impactEnergyOffset = pair_gran->get_history_value_offset("impactEnergy");
  forceMaxOffset = pair_gran->get_history_value_offset("forceMax");
  if (deltaMaxOffset < 0 || impactEnergyOffset < 0)
    error->fix_error(FLERR,this,"failed to find deltaMaxOffset/impactEnergy value in contact history");
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixBreakParticle::calc_insertion_properties()
{
  // error checks
  if (threshold < 0.)
    error->fix_error(FLERR,this,"you have to specify a threshold value");
  if (nflowrate > 0. || massflowrate > 0.)
    error->fix_error(FLERR,this,"specifying 'nflowrate' or 'massflowrate' is not allowed");
  if (ninsert > 0 || massinsert > 0.)
    error->fix_error(FLERR,this,"specifying 'nparticles' or 'mass' is not allowed");
  if (insert_every <= 0)
    error->fix_error(FLERR,this,"specifying 'every' must be > 0");

  // fix holding particle fragments
  if(fix_distribution->n_particletemplates() != 1)
    error->fix_error(FLERR,this,"fix of type particledistribution/discrete must hold exactly one template");
  if(strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/fragments"))
    error->fix_error(FLERR,this,"fix of type particledistribution/discrete must hold exactly one template of type fix particletemplate/fragments");

  fix_fragments = static_cast<FixTemplateFragments*>(fix_distribution->particletemplates()[0]);

  // do not need ninsert_per
}

/* ---------------------------------------------------------------------- */

int FixBreakParticle::setmask()
{
  int mask = FixInsert::setmask();
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

double FixBreakParticle::insertion_fraction()
{
  if(n_break_this == 0) return 0.0;
  return static_cast<double>(n_break_this_local)/static_cast<double>(n_break_this);
}

/* ---------------------------------------------------------------------- */

inline int FixBreakParticle::is_nearby(int i)
{
  // need not check overlap with existing particles since we
  // know space originally taken by deleted particles is free
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::pre_force(int)
{
  // set sibling contact flags and collision factor in contact_history
  if(n_break_this_local > 0) {
    double **x = atom->x;
    double *radius = atom->radius;
    int *mask = atom->mask;

    int inum = pair_gran->list->inum;
    int * ilist = pair_gran->list->ilist;
    int * numneigh = pair_gran->list->numneigh;
    int ** firstneigh = pair_gran->list->firstneigh;
    double ** firsthist = pair_gran->listgranhistory->firstdouble;

    for (int ii = 0; ii < inum; ++ii) {

      const int i = ilist[ii];

      if (mask[i] & groupbit) {

        const double xtmp = x[i][0];
        const double ytmp = x[i][1];
        const double ztmp = x[i][2];
        const double radi = radius[i];

        double * const allhist = firsthist[i];
        int * const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; ++jj) {

          const int j = jlist[jj];
          const double delx = xtmp - x[j][0];
          const double dely = ytmp - x[j][1];
          const double delz = ztmp - x[j][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          const double radj = radius[j];
          const double radsum = radi + radj;

          if (rsq < radsum * radsum) {

            if (mask[j] & groupbit) {
              double * contact_history = &allhist[dnum*jj];

              if (contact_history[deltaMaxOffset] == 0.0) {
                contact_history[siblingOffset] = 1.0;
                contact_history[collisionFactorOffset] = fix_collision_factor->vector_atom[i];
              }
            } else {
              error->warning(FLERR, "Inserted fragments overlap with alien particle");
            }
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::end_of_step()
{
  double *flag = fix_break->vector_atom;
  double *breaker = fix_breaker->vector_atom;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag;
  double *radius = atom->radius;

  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;

  int ** firstneigh = pair_gran->list->firstneigh;
  double ** firsthist = pair_gran->listgranhistory->firstdouble;

  switch (breakage_criterion) {
  case BC_ENERGY:
    {
      std::vector<double> breaker_energy(nlocal, 0.0);
      std::vector<int> breaker_tag(nlocal, 0);
      for (int ii = 0; ii < inum; ++ii) {
        const int i = ilist[ii];
        double * const allhist = firsthist[i];
        int * const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; ++jj) {
          const int j = jlist[jj];
          double * contact_history = &allhist[dnum*jj];
          const double impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[i]);
          const double impact_energy_limited_j = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[j]);
          if(breaker_energy[i] < impact_energy_limited_i) {
            breaker_energy[i] = impact_energy_limited_i;
            breaker_tag[i] = tag[j];
          }
          if(breaker_energy[j] < impact_energy_limited_j) {
            breaker_energy[j] = impact_energy_limited_j;
            breaker_tag[j] = tag[i];
          }
          flag[i] += impact_energy_limited_i;
          flag[j] += impact_energy_limited_j;
        }
        if(mask[i] & groupbit && radius[i] > min_break_rad && breaker_energy[i] > 0.0) {
          const double probability = 1.0 - exp(-fMat * 2.0*radius[i] * flag[i]);
          if (probability > random->uniform()) {
            breaker[i] = static_cast<double>(breaker_tag[i]);
          }
        }
      }
    }
    break;
  case BC_FORCE:
    {
      std::vector<double> forceMax(nlocal, 0.0);
      for (int ii = 0; ii < inum; ++ii) {
        const int i = ilist[ii];
        double * const allhist = firsthist[i];
        int * const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; ++jj) {
          const int j = jlist[jj];
          double * contact_history = &allhist[dnum*jj];

          forceMax[i] = std::max(forceMax[i], contact_history[forceMaxOffset]);
          forceMax[j] = std::max(forceMax[j], contact_history[forceMaxOffset]);
        }
        if(mask[i] & groupbit && radius[i] > min_break_rad && forceMax[i] > 0.0) {
          double probability = 1.0 - exp(-fMat * forceMax[i] / threshold);
          if (flag[i] == 0.0) {
            flag[i] = random->uniform();
          }
          if (probability > flag[i]) {
            flag[i] = -probability;  // sign indicates breakage
          }
        }
      }
    }
    break;
  case BC_VON_MISES:
      /// TODO
      break;
  default:
    break;
  }
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::pre_insert()
{
  int i,ibreak;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *flag = fix_break->vector_atom;
  double *breaker = fix_breaker->vector_atom;
  AtomVec *avec = atom->avec;

  // count # of particles to break

  n_break_this_local = 0;
  mass_break_this_local = 0.;

  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;
  int ** firstneigh = pair_gran->list->firstneigh;
  double ** firsthist = pair_gran->listgranhistory->firstdouble;

  switch (breakage_criterion) {
  case BC_ENERGY:
    {
      std::vector<bool> found_breaker(nlocal, false);
      for (int ii = 0; ii < inum; ++ii) {
        const int i = ilist[ii];
        double * const allhist = firsthist[i];
        int * const jlist = firstneigh[i];
        const int jnum = numneigh[i];

        for (int jj = 0; jj < jnum; ++jj) {
          const int j = jlist[jj];

          if (breaker[i] > 0 && atom->tag[j] == breaker[i]) {
            found_breaker[i] = true;
            double * contact_history = &allhist[dnum*jj];
            if (contact_history[deltaMaxOffset] < 0.0) {
              ++n_break_this_local;
              mass_break_this_local += rmass[i];
              flag[i] = -(1.0 - exp(-fMat * 2.0*radius[i] * flag[i])); // sign indicates breakage
            }
          }
          if (breaker[j] > 0 && atom->tag[i] == breaker[j]) {
            found_breaker[j] = true;
            double * contact_history = &allhist[dnum*jj];
            if (contact_history[deltaMaxOffset] < 0.0) {
              ++n_break_this_local;
              mass_break_this_local += rmass[j];
              flag[j] = -(1.0 - exp(-fMat * 2.0*radius[j] * flag[j])); // sign indicates breakage
            }
          }
        }
        if (mask[i] & groupbit) {
          if(breaker[i] > 0 && !found_breaker[i]) {
            fprintf(screen,"FixBreakParticle::pre_insert: breaker not found!\n");
          }
        }
      }
    }
    break;
  case BC_FORCE:
    for(int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit && flag[i] < 0.0) {
        ++n_break_this_local;
        mass_break_this_local += rmass[i];
      }
    }
    break;
  case BC_VON_MISES:
    /// TODO
    break;
  default:
    break;
  }

  // tally stats
  MPI_Sum_Scalar(n_break_this_local,n_break_this,world);
  n_break += n_break_this;
  MPI_Sum_Scalar(mass_break_this_local,mass_break_this,world);
  mass_break += mass_break_this;

  // virtual contact force calculation to restore overlap energy of contacting particles
  delta_v.clear();

  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];

    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj];

      if (flag[i] < 0.0 || flag[j] < 0.0) {

        virtual_x_i.clear();
        virtual_x_i.push_back(x[i][0]);
        virtual_x_i.push_back(x[i][1]);
        virtual_x_i.push_back(x[i][2]);
        virtual_v_i.clear();
        virtual_v_i.push_back(v[i][0]);
        virtual_v_i.push_back(v[i][1]);
        virtual_v_i.push_back(v[i][2]);

        virtual_x_j.clear();
        virtual_x_j.push_back(x[j][0]);
        virtual_x_j.push_back(x[j][1]);
        virtual_x_j.push_back(x[j][2]);
        virtual_v_j.clear();
        virtual_v_j.push_back(v[j][0]);
        virtual_v_j.push_back(v[j][1]);
        virtual_v_j.push_back(v[j][2]);

        virtual_f_ij.clear();
        virtual_f_ij.push_back(0.0);
        virtual_f_ij.push_back(0.0);
        virtual_f_ij.push_back(0.0);

        while(true) {
            if(virtual_force(i, j) <= 0.0) break;
            virtual_initial_integrate(i, j);
            if(virtual_force(i, j) <= 0.0) break;
            virtual_final_integrate(i, j);
        }

        {
          std::map<int, std::vector<double> >::iterator it;
          it = delta_v.find(i);
          if (it == delta_v.end()) {
            std::vector<double> delta_v0(3,0.0);
            delta_v[i] = delta_v0;
          }
          it = delta_v.find(j);
          if (it == delta_v.end()) {
            std::vector<double> delta_v0(3,0.0);
            delta_v[j] = delta_v0;
          }
        }

        delta_v[i][0] += (virtual_v_i[0] - v[i][0]);
        delta_v[i][1] += (virtual_v_i[1] - v[i][1]);
        delta_v[i][2] += (virtual_v_i[2] - v[i][2]);

        delta_v[j][0] += (virtual_v_j[0] - v[j][0]);
        delta_v[j][1] += (virtual_v_j[1] - v[j][1]);
        delta_v[j][2] += (virtual_v_j[2] - v[j][2]);
      }
    }
  }

  for (std::map<int, std::vector<double> >::iterator it=delta_v.begin(); it!=delta_v.end(); ++it) {
    v[it->first][0] += it->second[0];
    v[it->first][1] += it->second[1];
    v[it->first][2] += it->second[2];
  }

  // allocate breakage data
  if (breakdata) memory->destroy(breakdata);
  memory->create(breakdata,n_break_this_local,10,"FixBreakParticle::breakdata");

  // fill breakage data and remove particles
  const double eMF = 0.485;
  i = ibreak = 0;
  while (i < nlocal) {
    if (mask[i] & groupbit && flag[i] < 0.0) {
      // copy data needed for insertion
      vectorCopy3D(x[i],&breakdata[ibreak][0]);
      vectorCopy3D(v[i],&breakdata[ibreak][3]);
      vectorScalarMult3D(&breakdata[ibreak][3], eMF);
      breakdata[ibreak][6] = radius[i];
      breakdata[ibreak][7] = -flag[i];
      breakdata[ibreak][8] = atom->tag[i];
      breakdata[ibreak][9] = (1.0 - eMF*eMF) * 0.5 * rmass[i] * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      ++ibreak;

      // delete particle
      avec->copy(nlocal-1,i,1);
      --nlocal;
    } else {
      ++i;
    }
  }

  // update local and global # particles
  //NP wait with tags, global map etc
  //NP since done by fix insert anyway

  atom->nlocal = nlocal;
  double rlocal = static_cast<double>(atom->nlocal);
  MPI_Allreduce(&rlocal,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);

  fix_fragments->pre_insert(n_break_this_local, breakdata);

  // print stats
  print_stats_breakage_during();
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::print_stats_breakage_during()
{
  int step = update->ntimestep;

  if (me == 0 && n_break_this > 0) {
    if (screen) {
      fprintf(screen ,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
              n_break_this,mass_break_this,step,n_break,mass_break);
    }

    if (logfile) {
      fprintf(logfile,"Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
              n_break_this,mass_break_this,step,n_break,mass_break);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBreakParticle::calc_ninsert_this()
{
  // number of ptis to insert this timestep
  // will effectively insert n_break_this * n_fragments spheres
  return n_break_this;
}

/* ----------------------------------------------------------------------
   generate new particles at positions where old particles were deleted
   function is executed locally on each process as opposed to
   FixInsertPack::x_v_omega()

   overlap check is not needed since space around broken particles is empty

   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixBreakParticle::x_v_omega(int ninsert_this, int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this)
{
  double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4];
  int nins;
  ParticleToInsertFragments *pti;

  vectorZeroize3D(omega_ins);
  vectorZeroize4D(quat_ins);

  // local insertion
  double mass_inserted_this_local = 0.;
  int ninserted_this_local = 0;
  int ninserted_spheres_this_local = 0;

  // global insertion
  ninserted_this = ninserted_spheres_this = 0;
  mass_inserted_this = 0.;

  // n_break_this ptis with n_fragments spheres each
  // n_break_this * n_fragments spheres to be generated

  for (int iparticle = 0; iparticle < n_break_this_local; ++iparticle) {

    // get position, velocity and radius of broken particle
    vectorCopy3D(&breakdata[iparticle][0],pos_ins);
    vectorCopy3D(&breakdata[iparticle][3],v_ins);

    // get pti and scale it down with radius of broken particle
    pti = static_cast<ParticleToInsertFragments*>(fix_distribution->pti_list[iparticle]);
    pti->fix_property_atom_id = fix_collision_factor->id;
    nins = pti->set_x_v_omega(pos_ins,v_ins,omega_ins,quat_insert);

    // tally stats
    ninserted_spheres_this_local += nins;
    mass_inserted_this_local += pti->mass_ins;
    ++ninserted_this_local;
  }

  // tally stats, have to do this since operation is locally on each process
  // as opposed to e.g. FixInsertPack::x_v_omega()

  MPI_Sum_Scalar(ninserted_spheres_this_local,ninserted_spheres_this,world);
  MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
  MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);
  tag_max = atom->tag_max();
}


double FixBreakParticle::virtual_force(int i, int j)
{
  virtual_f_ij.clear();
  virtual_f_ij.push_back(0.0);
  virtual_f_ij.push_back(0.0);
  virtual_f_ij.push_back(0.0);
  double deltan = 0.0;

  if(true) { // hertz/break
    const int itype = atom->type[i];
    const int jtype = atom->type[j];
    double ri = atom->radius[i];
    double rj = atom->radius[j];
    double reff = ri*rj/(ri+rj);//cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
    double mi = atom->rmass[i];
    double mj = atom->rmass[j];
    double meff = mi * mj / (mi + mj);

    const double delx = virtual_x_i[0] - virtual_x_j[0];
    const double dely = virtual_x_i[1] - virtual_x_j[1];
    const double delz = virtual_x_i[2] - virtual_x_j[2];
    const double rsq = delx * delx + dely * dely + delz * delz;
    const double radsum = ri + rj;

    if (rsq < radsum * radsum) {
      const double r = sqrt(rsq);
      const double rinv = 1.0 / r;
      const double enx = delx * rinv;
      const double eny = dely * rinv;
      const double enz = delz * rinv;

      deltan = radsum - r;

      double sqrtval = sqrt(reff*deltan);

      const int max_type = pair_gran->mpg->max_type();
      const double *Y, *nu;
      const double * const * e;
      Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
      nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
      e = static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,style))->get_array();

      const double Yeff = 1./((1.-nu[itype-1]*nu[itype-1])/Y[itype-1]+(1.-nu[jtype-1]*nu[jtype-1])/Y[jtype-1]);
      const double loge = log(e[itype-1][jtype-1]);
      const double betaeff = loge/sqrt(loge * loge + MY_PI * MY_PI);
      const double Sn = 2.*Yeff*sqrtval;

      double kn = 4./3. * Yeff * sqrtval;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      double gamman = -2. * sqrtFiveOverSix * betaeff * sqrt(Sn*meff);

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;

      // relative translational velocity
      const double vr1 = virtual_v_i[0] - virtual_v_j[0];
      const double vr2 = virtual_v_i[1] - virtual_v_j[1];
      const double vr3 = virtual_v_i[2] - virtual_v_j[2];

      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double Fn_damping = -gamman * vn;
      const double Fn_contact = kn*deltan;
      double Fn = Fn_damping + Fn_contact;

      // limit force to avoid the artefact of negative repulsion force
      bool limitForce = false;
      if (limitForce && (Fn < 0.0) ) {
        Fn = 0.0;
      }

      // apply normal force
      if (false) { /*(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];*/
      } else {
        virtual_f_ij[0] = Fn * enx;
        virtual_f_ij[1] = Fn * eny;
        virtual_f_ij[2] = Fn * enz;
      }
    }
  } else { // hooke/break
    /// TODO
  }
  return deltan;
}


void FixBreakParticle::virtual_initial_integrate(int i, int j)
{
    const double dtv = update->dt;
    const double dtf = 0.5 * update->dt * force->ftm2v;
    double *rmass = atom->rmass;
    double dtfm;

    // update v and x of atoms

    dtfm = dtf / rmass[i];
    virtual_v_i[0] += dtfm * virtual_f_ij[0];
    virtual_v_i[1] += dtfm * virtual_f_ij[1];
    virtual_v_i[2] += dtfm * virtual_f_ij[2];
    virtual_x_i[0] += dtv * virtual_v_i[0];
    virtual_x_i[1] += dtv * virtual_v_i[1];
    virtual_x_i[2] += dtv * virtual_v_i[2];

    dtfm = dtf / rmass[j];
    virtual_v_j[0] -= dtfm * virtual_f_ij[0];
    virtual_v_j[1] -= dtfm * virtual_f_ij[1];
    virtual_v_j[2] -= dtfm * virtual_f_ij[2];
    virtual_x_j[0] += dtv * virtual_v_j[0];
    virtual_x_j[1] += dtv * virtual_v_j[1];
    virtual_x_j[2] += dtv * virtual_v_j[2];
}


void FixBreakParticle::virtual_final_integrate(int i, int j)
{
    const double dtf = 0.5 * update->dt * force->ftm2v;
    double *rmass = atom->rmass;
    double dtfm;

    // update v of atoms

    dtfm = dtf / rmass[i];
    virtual_v_i[0] += dtfm * virtual_f_ij[0];
    virtual_v_i[1] += dtfm * virtual_f_ij[1];
    virtual_v_i[2] += dtfm * virtual_f_ij[2];

    dtfm = dtf / rmass[j];
    virtual_v_j[0] -= dtfm * virtual_f_ij[0];
    virtual_v_j[1] -= dtfm * virtual_f_ij[1];
    virtual_v_j[2] -= dtfm * virtual_f_ij[2];
}


