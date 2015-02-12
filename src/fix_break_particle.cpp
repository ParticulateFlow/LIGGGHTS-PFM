/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
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
#include "input.h"
#include "variable.h"
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
#include "math_extra.h"
#include "fix_wall_gran.h"
#include "primitive_wall.h"
#include "tri_mesh_contacts.h"

using namespace MathConst;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LIGGGHTS::ContactModels;

enum{BC_ENERGY, BC_FORCE, BC_VON_MISES};
enum{BD_CONSTANT, BD_UNIFORM, BD_GAUSSIAN};
enum{NONE, CONSTANT, EQUAL, ATOM};

#define LMP_DEBUGMODE_FIX_BREAK_PARTICLE true
#define LMP_DEBUG_OUT_FIX_BREAK_PARTICLE screen

/* ---------------------------------------------------------------------- */

FixBreakParticle::FixBreakParticle(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();
  fMatstr = thresholdstr = NULL;
  fMatAtom = thresholdAtom = NULL;
  fMatstyle = thresholdstyle = NONE;
  maxatom1 = maxatom2 = 0;

  bool hasargs = true;
  while (iarg < narg && hasargs) {
    hasargs = false;
    if (strcmp(arg[iarg],"fMat") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing fMat value");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        fMatstr = new char[n];
        strcpy(fMatstr,&arg[iarg+1][2]);
      } else {
        fMat = atof(arg[iarg+1]);
        if(fMat <= 0. ) error->fix_error(FLERR,this,"'fMat' must be > 0");
        fMatstyle = CONSTANT;
      }
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"energy_threshold") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing energy_threshold value");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        thresholdstr = new char[n];
        strcpy(thresholdstr,&arg[iarg+1][2]);
      } else {
        threshold = atof(arg[iarg+1]);
        if (threshold < 0.) error->fix_error(FLERR,this,"'energy_threshold' must be >= 0");
        thresholdstyle = CONSTANT;
      }
      breakage_criterion = BC_ENERGY;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"force_threshold") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing force_threshold value");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        thresholdstr = new char[n];
        strcpy(thresholdstr,&arg[iarg+1][2]);
      } else {
        threshold = atof(arg[iarg+1]);
        if (threshold <= 0.) error->fix_error(FLERR,this,"'force_threshold' must be > 0");
        thresholdstyle = CONSTANT;
      }
      breakage_criterion = BC_FORCE;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"von_mises_stress") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: missing von_mises_stress value");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        thresholdstr = new char[n];
        strcpy(thresholdstr,&arg[iarg+1][2]);
      } else {
        threshold = atof(arg[iarg+1]);
        if (threshold <= 0.) error->fix_error(FLERR,this,"'von_mises_stress' must be > 0");
        thresholdstyle = CONSTANT;
      }
      breakage_criterion = BC_VON_MISES;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"breakability") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments: breakability");
      if (strcmp(arg[iarg+1],"constant") == 0) {
        if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments: breakability");
        const_breakability = atof(arg[iarg+2]);
        breakability_distribution = BD_CONSTANT;
        if (const_breakability < 0. || const_breakability >= 1.) error->fix_error(FLERR,this,"constant breakability must be >= 0 and < 1");
        ++iarg;
      } else if (strcmp(arg[iarg+1],"uniform") == 0) {
        breakability_distribution = BD_UNIFORM;
      } else if (strcmp(arg[iarg+1],"gaussian") == 0) {
        if (iarg+4 > narg) error->fix_error(FLERR,this,"not enough arguments: breakability");
        rand_expected_value = atof(arg[iarg+2]);
        rand_full_width = atof(arg[iarg+3]);
        if (rand_expected_value < 0. || rand_expected_value > 1.) error->fix_error(FLERR,this,"breakability must be >= 0 and < 1");
        breakability_distribution = BD_GAUSSIAN;
        iarg +=2;
      } else {
        error->fix_error(FLERR,this,"unknown option for: breakability");
      }
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
  fix_breakability = NULL;
  fix_breaker = NULL;
  fix_breaker_wall = NULL;
  fix_collision_factor = NULL;
  fix_stress = NULL;

  // execute end of step
  nevery = 1;

  n_break = n_break_this = n_break_this_local = 0;
  mass_break = mass_break_this = mass_break_this_local = 0.;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::post_create()
{
  if (!fix_break) { // breaking flag
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
  if (!fix_breakability) { // breakability
    char * breakvar_name = new char[14+strlen(id)];
    sprintf(breakvar_name,"breakability_%s",id);

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
    fix_breakability = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    delete [] breakvar_name;
  }
  if (!fix_breaker) { // holds the id of the atom that caused breakage
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
  if (!fix_breaker_wall) { // holds the hash of the wall id that causes breakage
    char * breakvar_name = new char[14+strlen(id)];
    sprintf(breakvar_name,"breaker_wall_%s",id);

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
    fix_breaker_wall = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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

  if (!fix_stress && breakage_criterion == BC_VON_MISES) {
    // stress tensor in Voigt notation:
    char * breakvar_name = new char[14+strlen(id)];
    sprintf(breakvar_name,"break_stress_%s",id);

    const char* fixarg[14];
    fixarg[0] = breakvar_name;
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = breakvar_name;
    fixarg[4] = "vector"; // 1 vector per particle to be registered
    fixarg[5] = "no";     // restart
    fixarg[6] = "no";     // communicate ghost
    fixarg[7] = "no";     // communicate rev
    fixarg[8] = "0.";     // sigma_x = sigma11
    fixarg[9] = "0.";     // sigma_y = sigma22
    fixarg[10]= "0.";     // sigma_z = sigma33
    fixarg[11]= "0.";     // tau_yz  = sigma23
    fixarg[12]= "0.";     // tau_xz  = sigma13
    fixarg[13]= "0.";     // tau_xy  = sigma12
    fix_stress = modify->add_fix_property_atom(14,const_cast<char**>(fixarg),style);
    delete [] breakvar_name;
  }
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::pre_delete(bool unfixflag)
{
  modify->delete_fix(fix_break->id);
  modify->delete_fix(fix_breakability->id);
  modify->delete_fix(fix_breaker->id);
  modify->delete_fix(fix_breaker_wall->id);
  modify->delete_fix(fix_collision_factor->id);
  if (fix_stress) modify->delete_fix(fix_stress->id);
}

/* ---------------------------------------------------------------------- */

FixBreakParticle::~FixBreakParticle()
{
  memory->destroy(breakdata);
  memory->destroy(fMatAtom);
  memory->destroy(thresholdAtom);
  delete [] fMatstr;
  delete [] thresholdstr;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::init_defaults()
{
  threshold = -1.0;
  fMat = 1.0;
  breakability_distribution = BD_UNIFORM;
  rand_expected_value = 0.5;
  rand_full_width = 1.0;
  const_breakability = 0.5;
  breakage_criterion = BC_ENERGY;
  fix_fragments = NULL;
  min_break_rad = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::init()
{
  FixInsert::init();

  if (fMatstr) {
    fMatvar = input->variable->find(fMatstr);
    if (fMatvar < 0)
      error->all(FLERR,"Variable name for fMat of fix break/particle does not exist");
    if (input->variable->equalstyle(fMatvar)) fMatstyle = EQUAL;
    else if (input->variable->atomstyle(fMatvar)) fMatstyle = ATOM;
    else error->all(FLERR,"Variable for fMat of fix break/particle is invalid style");
  }
  if (thresholdstr) {
    thresholdvar = input->variable->find(thresholdstr);
    if (thresholdvar < 0)
      error->all(FLERR,"Variable name for threshold of fix break/particle does not exist");
    if (input->variable->equalstyle(thresholdvar)) thresholdstyle = EQUAL;
    else if (input->variable->atomstyle(thresholdvar)) thresholdstyle = ATOM;
    else error->all(FLERR,"Variable for threshold of fix break/particle is invalid style");
  }

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  max_type = pair_gran->mpg->max_type();

  dnum = pair_gran->dnum();
  deltaMaxOffset = pair_gran->get_history_value_offset("deltaMax");
  siblingOffset = pair_gran->get_history_value_offset("sibling");
  collisionFactorOffset = pair_gran->get_history_value_offset("collisionFactor");
  impactEnergyOffset = pair_gran->get_history_value_offset("impactEnergy");
  forceMaxOffset = pair_gran->get_history_value_offset("forceMax");
  normalOffset = pair_gran->get_history_value_offset("enx");

  if (deltaMaxOffset < 0 || siblingOffset < 0 || collisionFactorOffset < 0 ||
      impactEnergyOffset < 0 || forceMaxOffset < 0 || normalOffset < 0)
    error->fix_error(FLERR,this,"failed to find value offset in contact history");
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixBreakParticle::calc_insertion_properties()
{
  // error checks
  if (!thresholdstr && threshold < 0.)
    error->fix_error(FLERR,this,"you have to specify a threshold value");
  if (nflowrate > 0. || massflowrate > 0.)
    error->fix_error(FLERR,this,"specifying 'nflowrate' or 'massflowrate' is not allowed");
  if (ninsert > 0 || massinsert > 0.)
    error->fix_error(FLERR,this,"specifying 'nparticles' or 'mass' is not allowed");
  if (insert_every <= 0)
    error->fix_error(FLERR,this,"specifying 'every' must be > 0");

  // fix holding particle fragments
  if (fix_distribution->n_particletemplates() != 1)
    error->fix_error(FLERR,this,"fix of type particledistribution/discrete must hold exactly one template");
  if (strcmp(fix_distribution->particletemplates()[0]->style,"particletemplate/fragments"))
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
  if (n_break_this == 0) return 0.0;
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
  if (n_break_this > 0) {
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

      if (mask[i] & groupbit && fix_collision_factor->vector_atom[i] != 1.0) {

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
              contact_history[siblingOffset] = 1.0;
              contact_history[siblingOffset+1] = 0.0;
              contact_history[deltaMaxOffset] = 0.0;
              contact_history[impactEnergyOffset] = 0.0;
              contact_history[collisionFactorOffset] = fix_collision_factor->vector_atom[i];
            } else {
              error->warning(FLERR, "Inserted fragments overlap with alien particle");
            }
          }
        }
      }
    }

    fix_collision_factor->set_all(1.0);
    n_break_this = 0;
  }
}

/* ---------------------------------------------------------------------- */

unsigned int JSHash(const std::string& str)
{
  unsigned int hash = 1315423911;

  for (std::size_t i = 0; i < str.length(); ++i) {
    hash ^= ((hash << 5) + str[i] + (hash >> 2));
  }

  return hash;
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::end_of_step()
{
  if (next_reneighbor-1 != update->ntimestep) return;

  modify->clearstep_compute();
  if (fMatstyle == EQUAL) {
    fMat = input->variable->compute_equal(fMatvar);
  } else if (fMatstyle == ATOM) {
    if (atom->nlocal > maxatom1) {
      maxatom1 = atom->nmax;
      memory->destroy(fMatAtom);
      memory->create(fMatAtom,maxatom1,"break/particle:fMatAtom");
    }
    input->variable->compute_atom(fMatvar,igroup,fMatAtom,1,0);
  }

  if (thresholdstyle == EQUAL) {
    threshold = input->variable->compute_equal(thresholdvar);
  } else if (thresholdstyle == ATOM) {
    if (atom->nlocal > maxatom2) {
      maxatom2 = atom->nmax;
      memory->destroy(thresholdAtom);
      memory->create(thresholdAtom,maxatom2,"break/particle:thresholdAtom");
    }
    input->variable->compute_atom(thresholdvar,igroup,thresholdAtom,1,0);
  }
  modify->addstep_compute(update->ntimestep + nevery);

  switch (breakage_criterion) {
  case BC_ENERGY:
    check_energy_criterion();
    break;
  case BC_FORCE:
    check_force_criterion();
    break;
  case BC_VON_MISES:
    check_von_mises_criterion();
    break;
  default:
    break;
  }
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::check_energy_criterion()
{
  double *flag = fix_break->vector_atom;
  double *breakability = fix_breakability->vector_atom;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag;
  double *radius = atom->radius;
  int nall = nlocal + atom->nghost;

  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;

  int ** firstneigh = pair_gran->list->firstneigh;
  double ** firsthist = pair_gran->listgranhistory->firstdouble;

  int n_wall_fixes = modify->n_fixes_style("wall/gran");

  std::vector<double> breaker_energy(nall, 0.0);
  std::vector<int> breaker_tag(nall, 0);
  // in case of impact energy criterion, we may already have a breaker but are still waiting for largest overlap
  std::vector<bool> already_doomed(nall, false);
  for (int i = 0; i < nall; ++i) {
    if (mask[i] & groupbit &&
        radius[i] > min_break_rad &&
        (fix_breaker->get_vector_atom_int(i) > 0 || fix_breaker_wall->get_vector_atom_int(i) > 0)) {
      already_doomed[i] = true;
    }
  }

  // particle - particle
  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    double * const allhist = firsthist[i];
    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj];
      double * contact_history = &allhist[dnum*jj];
      double impact_energy_limited_i;
      double impact_energy_limited_j;

      if (thresholdstyle == ATOM) {
        impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*thresholdAtom[i]/radius[i]);
        impact_energy_limited_j = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*thresholdAtom[j]/radius[j]);
      } else {
        impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[i]);
        impact_energy_limited_j = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[j]);
      }

      if (!already_doomed[i]) { // no breaker yet
        if (breaker_energy[i] < impact_energy_limited_i) {
          breaker_energy[i] = impact_energy_limited_i;
          breaker_tag[i] = tag[j];
        }
      }
      if (!already_doomed[j]) { // no breaker yet
        if (breaker_energy[j] < impact_energy_limited_j) {
          breaker_energy[j] = impact_energy_limited_j;
          breaker_tag[j] = tag[i];
        }
      }

      flag[i] += impact_energy_limited_i;
      flag[j] += impact_energy_limited_j;
    }
  }

  for (int i = 0; i < nlocal; ++i) {
    if (!already_doomed[i]) {
      if (mask[i] & groupbit && radius[i] > min_break_rad && breaker_energy[i] > 0.0) {

        double probability;
        if (fMatstyle == ATOM) {
          probability = 1.0 - exp(-fMatAtom[i] * 2.0*radius[i] * flag[i]);
        } else {
          probability = 1.0 - exp(-fMat        * 2.0*radius[i] * flag[i]);
        }

        if (breakability[i] == 0.0) {
          switch (breakability_distribution) {
          case BD_GAUSSIAN:
            breakability[i] = random->gaussian()*0.15*rand_full_width + rand_expected_value;
            break;
          case BD_CONSTANT:
            breakability[i] = const_breakability;
            break;
          default:
            breakability[i] = random->uniform();
            break;
          }
        }
        if (probability > breakability[i]) {
          fix_breaker->set_vector_atom_int(i, breaker_tag[i]);
        }
      }
    }
  }

  // particle - wall
  for (int ifix = 0; ifix < n_wall_fixes; ++ifix) {
    FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));

    if (fwg->is_mesh_wall()) {

      int n_FixMesh = fwg->n_meshes();

      for (int iMesh = 0; iMesh < n_FixMesh; iMesh++) {

        TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
        int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
        FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
        if (!fix_contact) continue;

        // get neighborList and numNeigh
        FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
        if (!meshNeighlist) continue;

        // loop owned and ghost triangles
        for (int iTri = 0; iTri < nTriAll; iTri++) {

          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for (int iCont = 0; iCont < numneigh; iCont++) {

            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (already_doomed[iPart]) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
            if (radius[iPart] < min_break_rad) continue;

            double *contact_history = NULL;
            {
              // get contact history of particle iPart and triangle idTri
              //NP depends on naming in fix_wall_gran!
              std::string fix_nneighs_name("n_neighs_mesh_");
              fix_nneighs_name += mesh->mesh_id();
              FixPropertyAtom* fix_nneighs = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name.c_str(),"property/atom","scalar",0,0,this->style));
              if (fix_nneighs) {
                int idTri = mesh->id(iTri);
                const int nneighs = fix_nneighs->get_vector_atom_int(iPart);
                for (int j = 0; j < nneighs; ++j) {
                  if (fix_contact->partner(iPart, j) == idTri) {
                    contact_history = fix_contact->contacthistory(iPart, j);
                    break;
                  }
                }
              }
            }

            if (contact_history) {
              double impact_energy_limited_i;

              if (thresholdstyle == ATOM) {
                impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*thresholdAtom[iPart]/radius[iPart]);
              } else {
                impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[iPart]);
              }
              if (impact_energy_limited_i > 0.0) {
                flag[iPart] += impact_energy_limited_i;

                if (breaker_energy[iPart] < impact_energy_limited_i) {
                  breaker_energy[iPart] = impact_energy_limited_i;
                  breaker_tag[iPart] = static_cast<int>(JSHash(fwg->id));

                  double probability;
                  if (fMatstyle == ATOM) {
                    probability = 1.0 - exp(-fMatAtom[iPart] * 2.0*radius[iPart] * flag[iPart]);
                  } else {
                    probability = 1.0 - exp(-fMat            * 2.0*radius[iPart] * flag[iPart]);
                  }

                  if (breakability[iPart] == 0.0) {
                    switch (breakability_distribution) {
                    case BD_GAUSSIAN:
                      breakability[iPart] = random->gaussian()*0.15*rand_full_width + rand_expected_value;
                      break;
                    case BD_CONSTANT:
                      breakability[iPart] = const_breakability;
                      break;
                    default:
                      breakability[iPart] = random->uniform();
                      break;
                    }
                  }
                  if (probability > breakability[iPart]) {
                    fix_breaker->set_vector_atom_int(iPart, 0); // remove any particle breaker
                    fix_breaker_wall->set_vector_atom_int(iPart, breaker_tag[iPart]);
                  }
                }
              }
            }
          }
        }
      }
    } else { // primitive wall

      double **c_history = NULL;

      if (fwg->dnum() > 0) {
        //NP depends on naming in fix_wall_gran!
        std::string hist_name("history_");
        hist_name += fwg->id;
        FixPropertyAtom *fix_history_primitive_ =
        static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",fwg->dnum(),0,fwg->style));
        if (fix_history_primitive_) {
          c_history = fix_history_primitive_->array_atom;
        }
      }

      if (c_history) {
        // loop neighbor list
        int *neighborList;
        int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

        for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++) {
          int iPart = *neighborList;
          // do not need to handle ghost particles
          if (iPart >= nlocal) continue;
          if (already_doomed[iPart]) continue;
          if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
          if (radius[iPart] < min_break_rad) continue;

          double *contact_history = c_history[iPart];
          double impact_energy_limited_i;

          if (thresholdstyle == ATOM) {
            impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*thresholdAtom[iPart]/radius[iPart]);
          } else {
            impact_energy_limited_i = std::max(0.0, contact_history[impactEnergyOffset] - 0.5*threshold/radius[iPart]);
          }
          flag[iPart] += impact_energy_limited_i;

          if (breaker_energy[iPart] < impact_energy_limited_i) {
            breaker_energy[iPart] = impact_energy_limited_i;
            breaker_tag[iPart] = static_cast<int>(JSHash(fwg->id));

            double probability;
            if (fMatstyle == ATOM) {
              probability = 1.0 - exp(-fMatAtom[iPart] * 2.0*radius[iPart] * flag[iPart]);
            } else {
              probability = 1.0 - exp(-fMat            * 2.0*radius[iPart] * flag[iPart]);
            }

            if (breakability[iPart] == 0.0) {
              switch (breakability_distribution) {
              case BD_GAUSSIAN:
                breakability[iPart] = random->gaussian()*0.15*rand_full_width + rand_expected_value;
                break;
              case BD_CONSTANT:
                breakability[iPart] = const_breakability;
                break;
              default:
                breakability[iPart] = random->uniform();
                break;
              }
            }
            if (probability > breakability[iPart]) {
              fix_breaker->set_vector_atom_int(iPart, 0); // remove any particle breaker
              fix_breaker_wall->set_vector_atom_int(iPart, breaker_tag[iPart]);
            }
          }
        }
      }
    }
  }

  // identify breaking particles
  {
    std::vector<bool> found_breaker(nlocal, false);

    // check particles
    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];
      double * const allhist = firsthist[i];
      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; ++jj) {
        const int j = jlist[jj];
        double * contact_history = &allhist[dnum*jj];

        if (i < nlocal && tag[j] == fix_breaker->get_vector_atom_int(i)) {
          found_breaker[i] = true;
          if (contact_history[deltaMaxOffset] < 0.0) {
            // sign indicates breakage
            if (fMatstyle == ATOM) {
              flag[i] = -(1.0 - exp(-fMatAtom[i] * 2.0*radius[i] * flag[i]));
            } else {
              flag[i] = -(1.0 - exp(-fMat        * 2.0*radius[i] * flag[i]));
            }
          }
        }
        if (j < nlocal && tag[i] == fix_breaker->get_vector_atom_int(j)) {
          found_breaker[j] = true;
          if (contact_history[deltaMaxOffset] < 0.0) {
            // sign indicates breakage
            if (fMatstyle == ATOM) {
              flag[j] = -(1.0 - exp(-fMatAtom[j] * 2.0*radius[j] * flag[j]));
            } else {
              flag[j] = -(1.0 - exp(-fMat        * 2.0*radius[j] * flag[j]));
            }
          }
        }
      }
    }

    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        if (fix_breaker->get_vector_atom_int(i) > 0 && !found_breaker[i]) {
          // breaker vanished from neighlist, break now
          // sign indicates breakage
          if (fMatstyle == ATOM) {
            flag[i] = -(1.0 - exp(-fMatAtom[i] * 2.0*radius[i] * flag[i]));
          } else {
            flag[i] = -(1.0 - exp(-fMat        * 2.0*radius[i] * flag[i]));
          }
        }
      }
    }

    // check walls
    for (int ifix = 0; ifix < n_wall_fixes; ++ifix) {

      FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));

      if (fwg->is_mesh_wall()) {
        int n_FixMesh = fwg->n_meshes();

        for (int iMesh = 0; iMesh < n_FixMesh; iMesh++) {

          TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
          int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
          FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
          if (!fix_contact) continue;

          // get neighborList and numNeigh
          FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
          if (!meshNeighlist) continue;

          // loop owned and ghost triangles
          for (int iTri = 0; iTri < nTriAll; iTri++) {

            const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
            const int numneigh = neighborList.size();

            for (int iCont = 0; iCont < numneigh; iCont++) {

              const int iPart = neighborList[iCont];

              // do not need to handle ghost particles
              if (iPart >= nlocal) continue;
              if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
              if (radius[iPart] < min_break_rad) continue;
              if (fix_breaker_wall->get_vector_atom_int(iPart) != static_cast<int>(JSHash(fwg->id))) continue;

              double *contact_history = NULL;
              {
                // get contact history of particle iPart and triangle idTri
                int idTri = mesh->id(iTri);
                std::string fix_nneighs_name("n_neighs_mesh_");
                fix_nneighs_name += mesh->mesh_id();
                FixPropertyAtom* fix_nneighs = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name.c_str(),"property/atom","scalar",0,0,this->style));
                const int nneighs = fix_nneighs->get_vector_atom_int(iPart);
                for (int j = 0; j < nneighs; ++j) {
                  if (fix_contact->partner(iPart, j) == idTri) {
                    contact_history = fix_contact->contacthistory(iPart, j);
                    break;
                  }
                }
              }

              if (contact_history) {
                if (contact_history[deltaMaxOffset] < 0.0) {
                  // sign indicates breakage
                  if (fMatstyle == ATOM) {
                    flag[iPart] = -(1.0 - exp(-fMatAtom[iPart] * 2.0*radius[iPart] * flag[iPart]));
                  } else {
                    flag[iPart] = -(1.0 - exp(-fMat            * 2.0*radius[iPart] * flag[iPart]));
                  }
                }
              }
            }
          }
        }
      } else { // primitive wall

        double **c_history = NULL;

        if (fwg->dnum() > 0) {
          std::string hist_name("history_");
          hist_name += fwg->id;
          FixPropertyAtom *fix_history_primitive_ =
          static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",fwg->dnum(),0,fwg->style));
          if (fix_history_primitive_) {
            c_history = fix_history_primitive_->array_atom;
          }
        }

        if (c_history) {
          // loop neighbor list
          int *neighborList;
          int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

          for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++) {
            int iPart = *neighborList;
            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
            if (radius[iPart] < min_break_rad) continue;
            if (fix_breaker_wall->get_vector_atom_int(iPart) != static_cast<int>(JSHash(fwg->id))) continue;

            double *contact_history = c_history[iPart];
            if (contact_history[deltaMaxOffset] < 0.0) {
              // sign indicates breakage
              if (fMatstyle == ATOM) {
                flag[iPart] = -(1.0 - exp(-fMatAtom[iPart] * 2.0*radius[iPart] * flag[iPart]));
              } else {
                flag[iPart] = -(1.0 - exp(-fMat            * 2.0*radius[iPart] * flag[iPart]));
              }
            }
          }
        }
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::check_force_criterion()
{
  double *flag = fix_break->vector_atom;
  double *breakability = fix_breakability->vector_atom;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *radius = atom->radius;
  int nall = nlocal + atom->nghost;

  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;

  int ** firstneigh = pair_gran->list->firstneigh;
  double ** firsthist = pair_gran->listgranhistory->firstdouble;

  int n_wall_fixes = modify->n_fixes_style("wall/gran");

  std::vector<double> forceMax(nall, 0.0);

  // particle - particle
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
  }

  // particle - wall
  for (int ifix = 0; ifix < n_wall_fixes; ++ifix) {
    FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
    if (fwg->is_mesh_wall()) {

      int n_FixMesh = fwg->n_meshes();

      for (int iMesh = 0; iMesh < n_FixMesh; iMesh++) {

        TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
        int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
        FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
        if (!fix_contact) continue;

        // get neighborList and numNeigh
        FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
        if (!meshNeighlist) continue;

        // loop owned and ghost triangles
        for (int iTri = 0; iTri < nTriAll; iTri++) {

          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for (int iCont = 0; iCont < numneigh; iCont++) {

            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
            if (radius[iPart] < min_break_rad) continue;

            double *contact_history = NULL;
            {
              // get contact history of particle iPart and triangle idTri
              //NP depends on naming in fix_wall_gran!
              std::string fix_nneighs_name("n_neighs_mesh_");
              fix_nneighs_name += mesh->mesh_id();
              FixPropertyAtom* fix_nneighs = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name.c_str(),"property/atom","scalar",0,0,this->style));
              if (fix_nneighs) {
                int idTri = mesh->id(iTri);
                const int nneighs = fix_nneighs->get_vector_atom_int(iPart);
                for (int j = 0; j < nneighs; ++j) {
                  if (fix_contact->partner(iPart, j) == idTri) {
                    contact_history = fix_contact->contacthistory(iPart, j);
                    break;
                  }
                }
              }
            }

            if (contact_history) {
              forceMax[iPart] = std::max(forceMax[iPart], contact_history[forceMaxOffset]);
            }
          }
        }
      }
    } else { // primitive wall

      double **c_history = NULL;

      if (fwg->dnum() > 0) {
        //NP depends on naming in fix_wall_gran!
        std::string hist_name("history_");
        hist_name += fwg->id;
        FixPropertyAtom *fix_history_primitive_ =
          static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",fwg->dnum(),0,fwg->style));
        if (fix_history_primitive_) {
          c_history = fix_history_primitive_->array_atom;
        }
      }

      if (c_history) {
        // loop neighbor list
        int *neighborList;
        int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

        for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++) {
          int iPart = *neighborList;
          // do not need to handle ghost particles
          if (iPart >= nlocal) continue;
          if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
          if (radius[iPart] < min_break_rad) continue;

          double *contact_history = c_history[iPart];

          forceMax[iPart] = std::max(forceMax[iPart], contact_history[forceMaxOffset]);
        }
      }
    }
  }

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit && radius[i] > min_break_rad && forceMax[i] > 0.0) {
      double probability;
      // P = 1 - exp(-fMat * (r/r0)^(3-2m) * (f/f0)^m)
      // m ... Weibull modulus (shape parameter)
      // here: m = 1; fMat is supposed to be scaled by 1/r0
      if (fMatstyle == ATOM) {
        if (thresholdstyle == ATOM) {
          probability = 1.0 - exp(-fMatAtom[i] * radius[i] * forceMax[i] / thresholdAtom[i]);
        } else {
          probability = 1.0 - exp(-fMatAtom[i] * radius[i] * forceMax[i] / threshold);
        }
      } else {
        if (thresholdstyle == ATOM) {
          probability = 1.0 - exp(-fMat        * radius[i] * forceMax[i] / thresholdAtom[i]);
        } else {
          probability = 1.0 - exp(-fMat        * radius[i] * forceMax[i] / threshold);
        }
      }

      if (breakability[i] == 0.0) {
        switch (breakability_distribution) {
        case BD_GAUSSIAN:
          breakability[i] = random->gaussian()*0.15*rand_full_width + rand_expected_value;
          break;
        case BD_CONSTANT:
          breakability[i] = const_breakability;
          break;
        default:
          breakability[i] = random->uniform();
          break;
        }
      }
      if (probability > breakability[i]) {
        flag[i] = -probability;  // sign indicates breakage
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::check_von_mises_criterion()
{
  double *flag = fix_break->vector_atom;
  double *breakability = fix_breakability->vector_atom;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double *radius = atom->radius;
  int nall = nlocal + atom->nghost;

  int inum = pair_gran->list->inum;
  int * ilist = pair_gran->list->ilist;
  int * numneigh = pair_gran->list->numneigh;

  int ** firstneigh = pair_gran->list->firstneigh;
  double ** firsthist = pair_gran->listgranhistory->firstdouble;

  int n_wall_fixes = modify->n_fixes_style("wall/gran");

  double **stress = fix_stress->array_atom;
  std::vector<double> r_over_vol(nall, 0.0);
  for (int i = 0; i < nall; ++i) {
    for (int j = 0; j < 6; ++j) {
      stress[i][j] = 0.0;
    }
    r_over_vol[i] = 1./(MY_4PI3 * radius[i] * radius[i]);
  }

  // particle - particle
  for (int ii = 0; ii < inum; ++ii) {
    const int i = ilist[ii];
    double * const allhist = firsthist[i];
    int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      const int j = jlist[jj];
      double * contact_history = &allhist[dnum*jj];

      double Fn = contact_history[forceMaxOffset];
      //NP en points from j to i
      double enx = contact_history[normalOffset];
      double eny = contact_history[normalOffset+1];
      double enz = contact_history[normalOffset+2];

      if (mask[i] & groupbit && radius[i] > min_break_rad) {
        stress[i][0] += r_over_vol[i] * enx * Fn * enx;
        stress[i][1] += r_over_vol[i] * eny * Fn * eny;
        stress[i][2] += r_over_vol[i] * enz * Fn * enz;
        stress[i][3] += r_over_vol[i] * eny * Fn * enz;
        stress[i][4] += r_over_vol[i] * enx * Fn * enz;
        stress[i][5] += r_over_vol[i] * enx * Fn * eny;
      }

      if (mask[j] & groupbit && radius[j] > min_break_rad) {
        stress[j][0] += r_over_vol[j] * enx * Fn * enx;
        stress[j][1] += r_over_vol[j] * eny * Fn * eny;
        stress[j][2] += r_over_vol[j] * enz * Fn * enz;
        stress[j][3] += r_over_vol[j] * eny * Fn * enz;
        stress[j][4] += r_over_vol[j] * enx * Fn * enz;
        stress[j][5] += r_over_vol[j] * enx * Fn * eny;
      }
    }
  }

  // particle - wall
  for (int ifix = 0; ifix < n_wall_fixes; ++ifix) {
    FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
    if (fwg->is_mesh_wall()) {

      int n_FixMesh = fwg->n_meshes();

      for (int iMesh = 0; iMesh < n_FixMesh; iMesh++) {

        TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
        int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
        FixContactHistoryMesh *fix_contact = fwg->mesh_list()[iMesh]->contactHistory();
        if (!fix_contact) continue;

        // get neighborList and numNeigh
        FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();
        if (!meshNeighlist) continue;

        // loop owned and ghost triangles
        for (int iTri = 0; iTri < nTriAll; iTri++) {

          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for (int iCont = 0; iCont < numneigh; iCont++) {

            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
            if (radius[iPart] < min_break_rad) continue;

            double *contact_history = NULL;
            {
              // get contact history of particle iPart and triangle idTri
              //NP depends on naming in fix_wall_gran!
              std::string fix_nneighs_name("n_neighs_mesh_");
              fix_nneighs_name += mesh->mesh_id();
              FixPropertyAtom* fix_nneighs = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name.c_str(),"property/atom","scalar",0,0,this->style));
              if (fix_nneighs) {
                int idTri = mesh->id(iTri);
                const int nneighs = fix_nneighs->get_vector_atom_int(iPart);
                for (int j = 0; j < nneighs; ++j) {
                  if (fix_contact->partner(iPart, j) == idTri) {
                    contact_history = fix_contact->contacthistory(iPart, j);
                    break;
                  }
                }
              }
            }

            if (contact_history) {
              double Fn = contact_history[forceMaxOffset];
              //NP en points from j to i
              double enx = contact_history[normalOffset];
              double eny = contact_history[normalOffset+1];
              double enz = contact_history[normalOffset+2];

              stress[iPart][0] += r_over_vol[iPart] * enx * Fn * enx;
              stress[iPart][1] += r_over_vol[iPart] * eny * Fn * eny;
              stress[iPart][2] += r_over_vol[iPart] * enz * Fn * enz;
              stress[iPart][3] += r_over_vol[iPart] * eny * Fn * enz;
              stress[iPart][4] += r_over_vol[iPart] * enx * Fn * enz;
              stress[iPart][5] += r_over_vol[iPart] * enx * Fn * eny;
            }
          }
        }
      }
    } else { // primitive wall

      double **c_history = NULL;

      if (fwg->dnum() > 0) {
        //NP depends on naming in fix_wall_gran!
        std::string hist_name("history_");
        hist_name += fwg->id;
        FixPropertyAtom *fix_history_primitive_ =
          static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",fwg->dnum(),0,fwg->style));
        if (fix_history_primitive_) {
          c_history = fix_history_primitive_->array_atom;
        }
      }

      if (c_history) {
        // loop neighbor list
        int *neighborList;
        int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

        for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++) {
          int iPart = *neighborList;
          // do not need to handle ghost particles
          if (iPart >= nlocal) continue;
          if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
          if (radius[iPart] < min_break_rad) continue;

          double *contact_history = c_history[iPart];

          double Fn = contact_history[forceMaxOffset];
          //NP en points from j to i
          double enx = contact_history[normalOffset];
          double eny = contact_history[normalOffset+1];
          double enz = contact_history[normalOffset+2];

          stress[iPart][0] += r_over_vol[iPart] * enx * Fn * enx;
          stress[iPart][1] += r_over_vol[iPart] * eny * Fn * eny;
          stress[iPart][2] += r_over_vol[iPart] * enz * Fn * enz;
          stress[iPart][3] += r_over_vol[iPart] * eny * Fn * enz;
          stress[iPart][4] += r_over_vol[iPart] * enx * Fn * enz;
          stress[iPart][5] += r_over_vol[iPart] * enx * Fn * eny;
        }
      }
    }
  }

  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit && radius[i] > min_break_rad) {
      // deviator stress
      double trace3 = (stress[i][0] + stress[i][1] + stress[i][2]) / 3.0;
      stress[i][0] -= trace3;
      stress[i][1] -= trace3;
      stress[i][2] -= trace3;
      double von_Mises_stress = 0.5 * (stress[i][0]*stress[i][0] + stress[i][1]*stress[i][1] + stress[i][2]*stress[i][2])
                                     + stress[i][3]*stress[i][3] + stress[i][4]*stress[i][4] + stress[i][5]*stress[i][5];
      if (von_Mises_stress > 0.0) {
        von_Mises_stress = sqrt(3.0 * von_Mises_stress);
        if (breakability[i] == 0.0) {
          switch (breakability_distribution) {
          case BD_GAUSSIAN:
           breakability[i] = random->gaussian()*0.15*rand_full_width + rand_expected_value;
            break;
          case BD_CONSTANT:
            breakability[i] = const_breakability;
            break;
          default:
            breakability[i] = random->uniform();
            break;
          }
        }
        double probability;
        // P = 1 - exp(-fMat * (r/r0)^3 * (sigma/sigma0)^m)
        // m ... Weibull modulus (shape parameter)
        // here: m = 1; fMat is supposed to be scaled by 1/r0^3
        if (fMatstyle == ATOM) {
          if (thresholdstyle == ATOM) {
            probability = 1.0 - exp(-fMatAtom[i] * radius[i]*radius[i]*radius[i] * (von_Mises_stress / thresholdAtom[i]));
          } else {
            probability = 1.0 - exp(-fMatAtom[i] * radius[i]*radius[i]*radius[i] * (von_Mises_stress / threshold));
          }
        } else {
          if (thresholdstyle == ATOM) {
            probability = 1.0 - exp(-fMat * radius[i]*radius[i]*radius[i] * (von_Mises_stress / thresholdAtom[i]));
          } else {
            probability = 1.0 - exp(-fMat * radius[i]*radius[i]*radius[i] * (von_Mises_stress / threshold));
          }
        }
        if (probability > breakability[i]) {
          flag[i] = -probability;  // sign indicates breakage
        }
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::pre_insert()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *flag = fix_break->vector_atom;
  int *tag = atom->tag;

  // count # of particles to break

  n_break_this_local = 0;
  mass_break_this_local = 0.;

  std::set<int> breaking_atom_tags_this;

  {
    // find out which particles are breaking this timestep
    std::vector<int> breaking_atom_tags_this_local;

    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit && flag[i] < 0.0) {
        ++n_break_this_local;
        mass_break_this_local += rmass[i];
        breaking_atom_tags_this_local.push_back(tag[i]);
      }
    }

    // tally stats
    MPI_Sum_Scalar(n_break_this_local,n_break_this,world);
    n_break += n_break_this;
    MPI_Sum_Scalar(mass_break_this_local,mass_break_this,world);
    mass_break += mass_break_this;

    if(n_break_this > 0) {
      int *breaking_atom_tags_this_recv = NULL;
      MPI_Allgather_Vector(&breaking_atom_tags_this_local[0], n_break_this_local, breaking_atom_tags_this_recv, world);
      for (int i = 0; i < n_break_this; ++i) {
        breaking_atom_tags_this.insert(breaking_atom_tags_this_recv[i]);
      }
      delete [] breaking_atom_tags_this_recv;
    }
  }

  std::multimap<int,std::vector<double> > contacting_atoms; // atom_tag, x & radius
  std::multimap<int,TriMeshContacts*> contacting_meshes; // atom_tag, mesh & triangles
  std::multimap<int,PrimitiveWall*> contacting_prim_walls;

  if (n_break_this > 0) {

    int inum = pair_gran->list->inum;
    int * ilist = pair_gran->list->ilist;
    int * numneigh = pair_gran->list->numneigh;
    int ** firstneigh = pair_gran->list->firstneigh;
    double ** firsthist = pair_gran->listgranhistory->firstdouble;

    int n_wall_fixes = modify->n_fixes_style("wall/gran");

    // estimate max fragment radius from max overlap
    std::map<int, double> deltaMax;
    for (int i = 0; i < nlocal; ++i)
      deltaMax[tag[i]] = 0.0;

    // virtual contact force calculation to restore overlap energy of contacting particles
    delta_v.clear();

    // particle - particle energy
    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];

      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      double * const allhist = firsthist[i];

      for (int jj = 0; jj < jnum; ++jj) {
        const int j = jlist[jj];

        if (breaking_atom_tags_this.find(tag[i]) != breaking_atom_tags_this.end() ||
            breaking_atom_tags_this.find(tag[j]) != breaking_atom_tags_this.end()) {

          const double delx = x[j][0] - x[i][0];
          const double dely = x[j][1] - x[i][1];
          const double delz = x[j][2] - x[i][2];
          const double rsq = delx * delx + dely * dely + delz * delz;
          const double radsum = radius[j] + radius[i];

          if (rsq < radsum * radsum) {
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

            if (breaking_atom_tags_this.find(tag[i]) != breaking_atom_tags_this.end()) {
              std::vector<double> xr_j(virtual_x_j);
              xr_j.push_back(radius[j]);
              contacting_atoms.insert(std::pair<int, std::vector<double> >(tag[i],xr_j));
            }
            if (breaking_atom_tags_this.find(tag[j]) != breaking_atom_tags_this.end()) {
              std::vector<double> xr_i(virtual_x_i);
              xr_i.push_back(radius[i]);
              contacting_atoms.insert(std::pair<int, std::vector<double> >(tag[j],xr_i));
            }

            double deltan = radsum - sqrt(rsq);
            deltaMax[tag[i]] = std::max(deltaMax[tag[i]], deltan);
            deltaMax[tag[j]] = std::max(deltaMax[tag[j]], deltan);

            double * contact_history = &allhist[dnum*jj];
            double siblingDeltaMax = contact_history[siblingOffset+1];
            double collision_factor = contact_history[collisionFactorOffset];

            while (true) {
              if (virtual_force(i, j, collision_factor, &siblingDeltaMax) <= 0.0) break;
              virtual_initial_integrate(i, virtual_v_i, virtual_x_i);
              MathExtra::negate3(&virtual_f_ij[0]);
              virtual_initial_integrate(j, virtual_v_j, virtual_x_j);

              if (virtual_force(i, j, collision_factor, &siblingDeltaMax) <= 0.0) break;
              virtual_final_integrate(i, virtual_v_i);
              MathExtra::negate3(&virtual_f_ij[0]);
              virtual_final_integrate(j, virtual_v_j);
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
    }

    // particle - wall energy
    if (n_break_this_local > 0) {
      for (int ifix = 0; ifix < n_wall_fixes; ++ifix) {
        FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
        const int64_t normalmodel = pair_gran->hashcode() & 0x0f;
        if (fwg->is_mesh_wall()) {

          double bary[3];
          double delta[3],deltan;

          for (int iMesh = 0; iMesh < fwg->n_meshes(); iMesh++) {
            TriMesh *mesh = fwg->mesh_list()[iMesh]->triMesh();
            int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

            // get neighborList and numNeigh
            FixNeighlistMesh * meshNeighlist = fwg->mesh_list()[iMesh]->meshNeighlist();

            int atom_type_wall = fwg->mesh_list()[iMesh]->atomTypeWall();

            std::map<int,std::set<int> > contacting_triangles; // atom_tag, triangle
            // loop owned and ghost triangles
            for (int iTri = 0; iTri < nTriAll; iTri++) {
              const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
              const int numneigh = neighborList.size();
              for (int iCont = 0; iCont < numneigh; iCont++) {
                const int iPart = neighborList[iCont];

                // do not need to handle ghost particles
                if (iPart >= nlocal) continue;
                if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
                if (breaking_atom_tags_this.find(tag[iPart]) == breaking_atom_tags_this.end()) continue;

                deltan = mesh->resolveTriSphereContactBary(iPart, iTri, radius[iPart], x[iPart], delta, bary);

                // to ensure that the fragments end up on the same side of the mesh as the original particle
                // all contacting triangles are required in ParticleSpatialDistribution
                // triangles that have an edge/corner contact but edge/corner are not active return (LARGE_TRIMESH - radius)
                if (deltan + radius[iPart] > 999999) {
                  contacting_triangles[tag[iPart]].insert(iTri);
                } else if (deltan < 0.0) {
                  contacting_triangles[tag[iPart]].insert(iTri);

                  deltan = -deltan;
                  deltaMax[tag[iPart]] = std::max(deltaMax[tag[iPart]], deltan);

                  double Yeff;
                  {
                    int itype = atom->type[iPart];
                    int jtype = atom_type_wall;
                    const double *Y, *nu;
                    Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
                    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
                    Yeff = 1./((1.-nu[itype-1]*nu[itype-1])/Y[itype-1]+(1.-nu[jtype-1]*nu[jtype-1])/Y[jtype-1]);
                  }

                  double E_el = 0.0;
                  if (normalmodel == 6) { // hertz/break
                    const double sqrtval = sqrt(radius[iPart]*deltan);
                    const double kn = 4./3. * Yeff * sqrtval;
                    E_el = 0.4*kn*deltan*deltan;
                  } else if (normalmodel == 7) { // hooke/break
                    const double sqrtval = sqrt(radius[iPart]);
                    const double meff = rmass[iPart];
                    const double charVel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,style))->compute_scalar();
                    const double kn = (16./15.) * sqrtval * Yeff * pow(15. * meff * charVel * charVel /(16. * sqrtval * Yeff), 0.2);
                    E_el = 0.5 * kn * deltan * deltan;
                  }

                  double v_el = sqrt(2.0 * E_el / rmass[iPart]);
                  const double rsq = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
                  const double r = sqrt(rsq);
                  const double rinv = 1.0 / r;
                  std::map<int, std::vector<double> >::iterator it;
                  it = delta_v.find(iPart);
                  if (it == delta_v.end()) {
                    std::vector<double> delta_v0(3,0.0);
                    delta_v[iPart] = delta_v0;
                  }
                  delta_v[iPart][0] += v_el * delta[0]*rinv;
                  delta_v[iPart][1] += v_el * delta[1]*rinv;
                  delta_v[iPart][2] += v_el * delta[2]*rinv;
                }
              }
            } // end triangle loop

            for (std::map<int,std::set<int> >::const_iterator it=contacting_triangles.begin(); it!=contacting_triangles.end(); ++it) {
              TriMeshContacts *tmc = new TriMeshContacts(mesh);
              tmc->contacts = it->second;
              contacting_meshes.insert(std::pair<int,TriMeshContacts*>(it->first, tmc));
            }
          }
        } else { // primitive wall

          // loop neighbor list
          int *neighborList;
          int nNeigh = fwg->primitiveWall()->getNeighbors(neighborList);

          for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++) {
            int iPart = *neighborList;

            // do not need to handle ghost particles
            if (iPart >= nlocal) continue;
            if (!(mask[iPart] & groupbit) || !(mask[iPart] & fwg->groupbit)) continue;
            if (breaking_atom_tags_this.find(tag[iPart]) == breaking_atom_tags_this.end()) continue;

            double delta[3];
            double deltan = fwg->primitiveWall()->resolveContact(x[iPart], radius[iPart], delta);

            if (deltan < 0.0) {
              contacting_prim_walls.insert(std::pair<int,PrimitiveWall*>(tag[iPart],fwg->primitiveWall()));

              deltan = -deltan;
              deltaMax[tag[iPart]] = std::max(deltaMax[tag[iPart]], deltan);

              double Yeff;
              {
                  int itype = atom->type[iPart];
                  int jtype = fwg->atom_type_wall();
                  const double *Y, *nu;
                  Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
                  nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
                  Yeff = 1./((1.-nu[itype-1]*nu[itype-1])/Y[itype-1]+(1.-nu[jtype-1]*nu[jtype-1])/Y[jtype-1]);
              }

              double E_el = 0.0;
              if (normalmodel == 6) { // hertz/break
                const double sqrtval = sqrt(radius[iPart]*deltan);
                const double kn = 4./3. * Yeff * sqrtval;
                E_el = 0.4*kn*deltan*deltan;
              } else if (normalmodel == 7) { // hooke/break
                const double sqrtval = sqrt(radius[iPart]);
                const double meff = rmass[iPart];
                const double charVel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,style))->compute_scalar();
                const double kn = (16./15.) * sqrtval * Yeff * pow(15. * meff * charVel * charVel /(16. * sqrtval * Yeff), 0.2);
                E_el = 0.5 * kn * deltan * deltan;
              }

              double v_el = sqrt(2.0 * E_el / rmass[iPart]);
              const double rsq = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
              const double r = sqrt(rsq);
              const double rinv = 1.0 / r;
              std::map<int, std::vector<double> >::iterator it;
              it = delta_v.find(iPart);
              if (it == delta_v.end()) {
                std::vector<double> delta_v0(3,0.0);
                delta_v[iPart] = delta_v0;
              }
              delta_v[iPart][0] += v_el * delta[0]*rinv;
              delta_v[iPart][1] += v_el * delta[1]*rinv;
              delta_v[iPart][2] += v_el * delta[2]*rinv;
            }
          }
        }
      }
    }

    // correct velocities
    for (std::map<int, std::vector<double> >::iterator it=delta_v.begin(); it!=delta_v.end(); ++it) {
      v[it->first][0] += it->second[0];
      v[it->first][1] += it->second[1];
      v[it->first][2] += it->second[2];
    }

    // allocate breakage data if required
    if (breakdata) {
      memory->destroy(breakdata);
      breakdata = NULL;
    }

    if (n_break_this_local > 0) {

      memory->create(breakdata,n_break_this_local,11,"FixBreakParticle::breakdata");

      // fill breakage data and remove particles
      //NP leave ghosts alone
      AtomVec *avec = atom->avec;
      const double eMF = 0.485;
      int i=0, ibreak=0;
      while (i < nlocal) {
        if ((mask[i] & groupbit) && (flag[i] < 0.0)) {
          // copy data needed for insertion
          vectorCopy3D(x[i],&breakdata[ibreak][0]);
          vectorCopy3D(v[i],&breakdata[ibreak][3]);
          vectorScalarMult3D(&breakdata[ibreak][3], eMF);
          breakdata[ibreak][6] = radius[i];
          breakdata[ibreak][7] = -flag[i];
          breakdata[ibreak][8] = tag[i];
          breakdata[ibreak][9] = (1.0 - eMF*eMF) * 0.5 * rmass[i] * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
          breakdata[ibreak][10] = deltaMax[tag[i]];

          ++ibreak;

          // delete particle
          avec->copy(nlocal-1,i,1);
          --nlocal;
        } else {
          ++i;
        }
      }

      // update local # particles
      atom->nlocal = nlocal;
    }
  }

  // update global # particles
  //NP wait with tags, global map etc
  //NP since done by fix insert anyway
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  fix_fragments->pre_insert(n_break_this_local, breakdata, contacting_atoms, contacting_prim_walls, contacting_meshes);
  for (std::multimap<int,TriMeshContacts*>::iterator it=contacting_meshes.begin(); it != contacting_meshes.end(); ++it)
    delete (it->second);

  // print stats
  print_stats_breakage_during();
}

/* ---------------------------------------------------------------------- */

void FixBreakParticle::print_stats_breakage_during()
{
  if (me == 0 && n_break_this > 0) {
    int step = update->ntimestep;

    if (screen) {
      fprintf(screen ,"INFO: Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
              n_break_this,mass_break_this,step,n_break,mass_break);
    }

    if (logfile) {
      fprintf(logfile,"INFO: Particle breakage: broke %d particles (mass %f) at step %d\n - a total of %d particles (mass %f) broken so far.\n",
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

void FixBreakParticle::x_v_omega(int ninsert_this_local, int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
  double pos_ins[3],v_ins[3],omega_ins[3],quat_ins[4];
  int nins;
  ParticleToInsertFragments *pti;

  vectorZeroize3D(omega_ins);
  vectorZeroize4D(quat_ins);

  ninserted_this_local = ninserted_spheres_this_local = 0;
  mass_inserted_this_local = 0.;

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
}


double FixBreakParticle::virtual_force(int i, int j, double collision_factor, double *siblingDeltaMax)
{
  virtual_f_ij.clear();
  virtual_f_ij.push_back(0.0);
  virtual_f_ij.push_back(0.0);
  virtual_f_ij.push_back(0.0);

  CollisionData cdata;

  const double delx = virtual_x_i[0] - virtual_x_j[0];
  const double dely = virtual_x_i[1] - virtual_x_j[1];
  const double delz = virtual_x_i[2] - virtual_x_j[2];

  cdata.rsq = delx * delx + dely * dely + delz * delz;
  cdata.radi = atom->radius[i];
  cdata.radj = atom->radius[j];
  cdata.radsum = cdata.radi + cdata.radj;

  if (cdata.rsq < cdata.radsum * cdata.radsum) {
    cdata.i = i;
    cdata.j = j;
    cdata.is_wall = false;
    cdata.itype = atom->type[i];
    cdata.jtype = atom->type[j];
    double mi = atom->rmass[i];
    double mj = atom->rmass[j];
    cdata.meff = mi * mj / (mi + mj);
    cdata.r = sqrt(cdata.rsq);
    const double rinv = 1.0 / cdata.r;
    cdata.rinv = rinv;
    cdata.deltan = cdata.radsum - cdata.r;
    cdata.en[0]   = delx * rinv;
    cdata.en[1]   = dely * rinv;
    cdata.en[2]   = delz * rinv;

    return virtual_collision(cdata, collision_factor, siblingDeltaMax);
  }

  return 0.0;
}


double FixBreakParticle::virtual_collision(CollisionData & cdata, double collision_factor, double *siblingDeltaMax)
{
  double deltan = 0.0;
  const int64_t normalmodel = pair_gran->hashcode() & 0x0f;

  if (normalmodel == 6) { // hertz/break
    const int itype = cdata.itype;
    const int jtype = cdata.jtype;
    const double ri = cdata.radi;
    const double rj = cdata.radj;
    const double reff = cdata.is_wall ? ri : (ri*rj/(ri+rj));
    const double meff = cdata.meff;

    {
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      deltan = cdata.deltan;

      if (!cdata.is_wall && *siblingDeltaMax > 0.0) {
        if (deltan > *siblingDeltaMax) {
          deltan -= *siblingDeltaMax;
          deltan += *siblingDeltaMax*collision_factor;
        } else {
          *siblingDeltaMax = deltan;
          deltan *= collision_factor;
        }
      }

      double sqrtval = sqrt(reff*deltan);

      double Yeff;
      double betaeff;
      {
        const double *Y, *nu;
        const double * const * e;
        Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
        nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
        e = static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,style))->get_array();

        Yeff = 1./((1.-nu[itype-1]*nu[itype-1])/Y[itype-1]+(1.-nu[jtype-1]*nu[jtype-1])/Y[jtype-1]);
        const double loge = log(e[itype-1][jtype-1]);
        betaeff = loge/sqrt(loge * loge + MY_PI * MY_PI);
      }

      const double Sn = 2. * Yeff*sqrtval;

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
      const double Fn_contact = kn * deltan;
      double Fn = Fn_damping + Fn_contact;

      // limit force to avoid the artefact of negative repulsion force
      bool limitForce = false;
      if (limitForce && (Fn < 0.0) ) {
        Fn = 0.0;
      }

      // apply normal force
      if (cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        virtual_f_ij[0] = Fn_ * enx;
        virtual_f_ij[1] = Fn_ * eny;
        virtual_f_ij[2] = Fn_ * enz;
      } else {
        virtual_f_ij[0] = Fn * enx;
        virtual_f_ij[1] = Fn * eny;
        virtual_f_ij[2] = Fn * enz;
      }
    }
  } else if (normalmodel == 7) { // hooke/break
    const int itype = cdata.itype;
    const int jtype = cdata.jtype;
    const double ri = cdata.radi;
    const double rj = cdata.radj;
    const double reff=cdata.is_wall ? ri : (ri*rj/(ri+rj));
    const double meff=cdata.meff;
    double coeffRestLogChosen;

    {
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      deltan = cdata.deltan;

      if (!cdata.is_wall && *siblingDeltaMax > 0.0) {
        if (deltan > *siblingDeltaMax) {
          deltan -= *siblingDeltaMax;
          deltan += *siblingDeltaMax*collision_factor;
        } else {
          *siblingDeltaMax = deltan;
          deltan *= collision_factor;
        }
      }

      const double sqrtval = sqrt(reff);

      double Yeff;
      const double charVel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,style))->compute_scalar();
      double coeffRestLog;
      {
        const double *Y, *nu;
        const double * const * e;
        Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
        nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
        e = static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,style))->get_array();

        Yeff = 1./((1.-nu[itype-1]*nu[itype-1])/Y[itype-1]+(1.-nu[jtype-1]*nu[jtype-1])/Y[jtype-1]);
        coeffRestLog = log(e[itype-1][jtype-1]);
      }

      /*if (viscous)  {
        // Stokes Number from MW Schmeeckle (2001)
        const double stokes=cdata.meff*cdata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
        // Empirical from Legendre (2006)
        coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
      } else*/ {
        coeffRestLogChosen=coeffRestLog;/// coeffRestLog[itype][jtype];
      }

      double kn = 16./15.*sqrtval*Yeff*pow(15.*meff*charVel*charVel/(16.*sqrtval*Yeff),0.2);
      const double gamman=sqrt(4.*meff*kn/(1.+(M_PI/coeffRestLogChosen)*(M_PI/coeffRestLogChosen)));

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;

      // relative translational velocity
      const double vr1 = virtual_v_i[0] - virtual_v_j[0];
      const double vr2 = virtual_v_i[1] - virtual_v_j[1];
      const double vr3 = virtual_v_i[2] - virtual_v_j[2];

      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double Fn_damping = -gamman * vn;
      const double Fn_contact = kn * deltan;
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      bool limitForce = false;
      if (limitForce && (Fn<0.0) ) {
        Fn = 0.0;
      }

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        virtual_f_ij[0] = Fn_ * enx;
        virtual_f_ij[1] = Fn_ * eny;
        virtual_f_ij[2] = Fn_ * enz;
      } else {
        virtual_f_ij[0] = Fn * enx;
        virtual_f_ij[1] = Fn * eny;
        virtual_f_ij[2] = Fn * enz;
      }
    }
  }
  return deltan;
}


void FixBreakParticle::virtual_initial_integrate(int i, std::vector<double> &virtual_v, std::vector<double> &virtual_x)
{
    const double dtv = update->dt;
    const double dtf = 0.5 * update->dt * force->ftm2v;
    double *rmass = atom->rmass;
    double dtfm;

    // update v and x of atom

    dtfm = dtf / rmass[i];
    virtual_v[0] += dtfm * virtual_f_ij[0];
    virtual_v[1] += dtfm * virtual_f_ij[1];
    virtual_v[2] += dtfm * virtual_f_ij[2];
    virtual_x[0] += dtv * virtual_v[0];
    virtual_x[1] += dtv * virtual_v[1];
    virtual_x[2] += dtv * virtual_v[2];
}


void FixBreakParticle::virtual_final_integrate(int i, std::vector<double> &virtual_v)
{
    const double dtf = 0.5 * update->dt * force->ftm2v;
    double *rmass = atom->rmass;
    double dtfm;

    // update v of atom

    dtfm = dtf / rmass[i];
    virtual_v[0] += dtfm * virtual_f_ij[0];
    virtual_v[1] += dtfm * virtual_f_ij[1];
    virtual_v[2] += dtfm * virtual_f_ij[2];
}

