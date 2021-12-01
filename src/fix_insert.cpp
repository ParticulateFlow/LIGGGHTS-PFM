/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2015-     JKU Linz

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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_multisphere.h"
#include "fix_particledistribution.h"
#include "fix_template_sphere.h"
#include "fix_property_atom.h"
#include "fix_insert.h"
#include "math_extra_liggghts.h"
#include "mpi_liggghts.h"
#include "vector_liggghts.h"
#include "volume_mesh.h"
#include "probability_distribution.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

#define LMP_DEBUGMODE_FIXINSERT false //(667 == update->ntimestep)//  true
#define LMP_DEBUG_OUT_FIXINSERT screen

/* ---------------------------------------------------------------------- */

FixInsert::FixInsert(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  neighList(lmp)
{
  if (narg < 7) error->fix_error(FLERR,this,"not enough arguments");

  //NP prohibit changing timestep size
  //NP not sure if this is necessary at all?
  time_depend = 1;

  restart_global = 1;

  setup_flag = false;

  /*NL*/ //VolumeMesh<4,4,3> *vm;
  /*NL*/ //TetMesh *tm;

  fix_distribution = NULL;
  fix_multisphere = NULL;
  multisphere = NULL;

  // required args
  iarg = 3;

  if(strcmp(arg[iarg++],"seed")) error->fix_error(FLERR,this,"expecting keyword 'seed'");
  seed = atoi(arg[iarg++]);
  if (seed <= 0) error->fix_error(FLERR,this,"illegal seed");

  // random number generator, seed independent of proc
  randomAll = new RanPark(lmp,seed);

  seed += comm->me;
  // random number generator, seed depends on proc
  random = new RanPark(lmp,seed);

  // set defaults
  init_defaults();

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  check_obb_flag = 1;
#endif

  // parse args
  //NP args processed by this class parsed here
  //NP let derived classes parse args rest of args
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"distributiontemplate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int ifix = modify->find_fix(arg[iarg+1]);
      if(ifix < 0 || strncmp(modify->fix[ifix]->style,"particledistribution", 20))
        error->fix_error(FLERR,this,"Fix insert requires you to define a valid ID for a fix of type particledistribution/discrete");
      fix_distribution = static_cast<FixParticledistribution*>(modify->fix[ifix]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"maxattempt") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      maxattempt = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"nparticles") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"INF") == 0)
        ninsert_exists = 0;
      else ninsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"INF") == 0)
        ninsert_exists = 0;
      else massinsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"massrate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      massflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"particlerate") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      nflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"insert_every") == 0 || strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"once") == 0) insert_every = 0;
      else insert_every = atoi(arg[iarg+1]);
      if(insert_every < 0) error->fix_error(FLERR,this,"insert_every must be >= 0");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      first_ins_step = atoi(arg[iarg+1]);
      if(first_ins_step < update->ntimestep + 1 && !modify->fix_restart_in_progress())
        error->fix_error(FLERR,this,"'start' step can not be before current step");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"overlapcheck") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes")==0) check_ol_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) check_ol_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"all_in") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes")==0) all_in_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) all_in_flag = 0;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
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
    } else if (strcmp(arg[iarg],"random_distribute") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"uncorrelated")==0) exact_number = 0;
      else if(strcmp(arg[iarg+1],"exact")==0) exact_number = 1;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"verbose") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"no")==0) print_stats_during_flag = 0;
      else if(strcmp(arg[iarg+1],"yes")==0) print_stats_during_flag = 1;
      else error->fix_error(FLERR,this,"");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+5 > narg) error->fix_error(FLERR,this,"not enough keyword for 'vel'");
      if (strcmp(arg[iarg+1],"constant") == 0)  {
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          iarg += 5;
      } else if (strcmp(arg[iarg+1],"uniform") == 0) {
          if (iarg+8 > narg) error->fix_error(FLERR,this,"not enough keyword for 'uniform'");
          v_randomSetting = RANDOM_UNIFORM;
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          v_insertFluct[0] = atof(arg[iarg+5]);
          v_insertFluct[1] = atof(arg[iarg+6]);
          v_insertFluct[2] = atof(arg[iarg+7]);
          iarg += 8;
      } else if (strcmp(arg[iarg+1],"gaussian") == 0) {
          if (iarg+8 > narg) error->fix_error(FLERR,this,"not enough keyword for 'gaussian'");
          v_randomSetting = RANDOM_GAUSSIAN;
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
          v_insertFluct[0] = atof(arg[iarg+5]);
          v_insertFluct[1] = atof(arg[iarg+6]);
          v_insertFluct[2] = atof(arg[iarg+7]);
          iarg += 8;
      } else
          error->fix_error(FLERR,this,"expecting keyword 'constant' or 'uniform' or 'gaussian' after keyword 'vel'");
      hasargs = true;
    } else if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+5 > narg) error->fix_error(FLERR,this,"");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          omega_insert[0] = atof(arg[iarg+2]);
          omega_insert[1] = atof(arg[iarg+3]);
          omega_insert[2] = atof(arg[iarg+4]);
      } else error->fix_error(FLERR,this,"expecting keyword 'constant' after keyword 'omega'");
      iarg += 5;
      hasargs = true;
    } else if (strcmp(arg[iarg],"orientation") == 0) {
      if (iarg+2 > narg)
        error->fix_error(FLERR,this,"not enough arguments for 'orientation'");
      iarg++;
      if(strcmp(arg[iarg],"random") == 0)
      {
          quat_random_ = true;
          iarg++;
      }
      else if(strcmp(arg[iarg],"template") == 0)
      {
          quat_random_ = false;
          iarg++;
      }
      else if (strcmp(arg[iarg],"constant") == 0)
      {
          iarg++;
          if (iarg+4 > narg) error->fix_error(FLERR,this,"");
          quat_insert[0] = atof(arg[iarg++]);
          quat_insert[1] = atof(arg[iarg++]);
          quat_insert[2] = atof(arg[iarg++]);
          quat_insert[3] = atof(arg[iarg++]);
      } else error->fix_error(FLERR,this,"expecting 'random', template' or 'constant' after keyword 'quat'");
      hasargs = true;
    }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    else if (strcmp(arg[iarg],"check_obb") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if(strcmp(arg[iarg+1],"yes")==0) check_obb_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) check_obb_flag = 0;
      else error->fix_error(FLERR,this,"");
      if(check_ol_flag==0) check_obb_flag = 0;
      iarg += 2;
      hasargs = true;
    }
#endif
    //NP throw error only if not derived class
    else if(strcmp(style,"insert") == 0) error->fix_error(FLERR,this,"unknown keyword");
  }

  // memory not allocated initially
  ninsert_this_max_local = 0;

  // allgather arrays
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // set next reneighbor
  force_reneighbor = 1;
  next_reneighbor = first_ins_step;
  most_recent_ins_step = -1;

  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;

  print_stats_start_flag = 1;

}

/* ---------------------------------------------------------------------- */

FixInsert::~FixInsert()
{
  delete random;
  delete randomAll;
  delete [] recvcounts;
  delete [] displs;
  delete [] property_name;
}

/* ---------------------------------------------------------------------- */

void FixInsert::post_create()
{
  // check for missing or contradictory settings
  sanity_check();

  //min/max type to be inserted, need that to check if material properties defined for all materials
  type_max = fix_distribution->max_type();
  type_min = fix_distribution->min_type();

  // calc max insertion radius
  int ntypes = atom->ntypes;
  maxrad = 0.;
  minrad = 1000.;
  for(int i = 1; i <= ntypes; i++)
  {
     maxrad = MathExtraLiggghts::max(maxrad,max_rad(i));
     minrad = MathExtraLiggghts::min(minrad,min_rad(i));
  }
}

/* ---------------------------------------------------------------------- */

void FixInsert::setup(int vflag)
{
  //NP do these things here to let derived classes parse args first
  //NP do not do this in post_create() since mesh data is not available in
  //NP parallel at this point

  // do this only once
  if(setup_flag) return;
  else setup_flag = true;

  // calculate ninsert, insert_every, ninsert_per
  calc_insertion_properties();

  // calc last step of insertion
  if(ninsert_exists)
  {
      if(ninsert <= ninsert_per)
        final_ins_step = first_ins_step;
      else
        final_ins_step = first_ins_step +
                static_cast<int>(static_cast<double>(ninsert)/ninsert_per) *  static_cast<double>(insert_every);

      /*NL*/ //if (screen) fprintf(screen,"ninsert %d, ninsert_per %f insert_every %d first_ins_step %d final_ins_step %d\n",ninsert,ninsert_per,insert_every,first_ins_step,final_ins_step);

      if(final_ins_step < 0)
        error->fix_error(FLERR,this,"Particle insertion: Overflow - need too long for particle insertion. "
                                    "Please decrease # particles to insert or increase insertion rate");
      if(ninsert < 0)
        error->fix_error(FLERR,this,"Particle insertion: Overflow - too many particles for particle insertion. "
                                    "Please decrease # particles to insert.");
  }
  else
    final_ins_step = -1;

  // print statistics
  print_stats_start();

}

/* ---------------------------------------------------------------------- */

void FixInsert::init_defaults()
{
  // default is that total # of particles to insert by this command is known
  ninsert_exists = 1;

  ninsert = ninserted = 0;
  massinsert = massinserted = 0.;
  nflowrate = massflowrate = 0.;

  insert_every = -1;
  ninsert_per = 0.;

  // 1st insertion on next timestep is default
  first_ins_step = update->ntimestep + 1;

  maxattempt = 50;

  check_ol_flag = 1;
  all_in_flag = 0;

  exact_number = 1;

  v_randomSetting = RANDOM_CONSTANT;
  vectorZeroize3D(v_insert);
  vectorZeroize3D(v_insertFluct);
  vectorZeroize3D(omega_insert);

  //NP initialize as unit maxtrix
  quatUnitize4D(quat_insert);
  quat_random_ = false;

  print_stats_during_flag = 1;
  warn_boxentent = true;

  property_name = 0;
  fix_property = 0;
  fix_property_value = 0.;
  fix_property_ivalue = 0;
  property_index = -1;
  property_iindex = -1;
}

/* ---------------------------------------------------------------------- */

void FixInsert::sanity_check()
{
    if(fix_distribution == NULL)
      error->fix_error(FLERR,this,"have to define a 'distributiontemplate'");

    if(MathExtraLiggghts::abs(vectorMag4DSquared(quat_insert)-1.) > 1e-10)
      error->fix_error(FLERR,this,"quaternion not valid");

    if(ninsert > 0 && massinsert > 0.)
      error->fix_error(FLERR,this,"must not define both 'nparticles' and 'mass'");
    if(nflowrate > 0. && massflowrate > 0.)
      error->fix_error(FLERR,this,"must not define both 'particlerate' and 'massrate'");

    if(insert_every == 0 && (massflowrate > 0. || nflowrate > 0.))
      error->fix_error(FLERR,this,"must not define 'particlerate' or 'massrate' for 'insert_every' = 0");
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_start()
{
  if (me == 0 && print_stats_start_flag) {

    if(ninsert_exists)
    {
        if (screen)
            fprintf(screen ,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f  (mass rate %f)\n"
                            "      %d particles (mass %f) within %d steps\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);

        if (logfile)
            fprintf(logfile,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f, (mass rate %f)\n"
                            "      %d particles (mass %f) within %d steps\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);
    }
    else if(massflowrate > 0.)
    {
        if (screen)
            fprintf(screen ,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f  (mass rate %f)\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate);

        if (logfile)
            fprintf(logfile,"INFO: Particle insertion %s: %f particles every %d steps - particle rate %f, (mass rate %f)\n",
                id,ninsert_per,insert_every,nflowrate,massflowrate);
    }
    else
    {
        if (screen)
            fprintf(screen ,"INFO: Particle insertion %s: inserting every %d steps\n",id,insert_every);

        if (logfile)
            fprintf(logfile ,"INFO: Particle insertion %s: inserting every %d steps\n",id,insert_every);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_during(int ninsert_this, double mass_inserted_this)
{
  if (me == 0 && print_stats_during_flag)
  {
    bigint step = update->ntimestep;

    if (screen)
      fprintf(screen ,"INFO: Particle insertion %s: inserted %d particle templates (mass %f) at step " BIGINT_FORMAT "\n - a total of %d particle templates (mass %f) inserted so far.\n",
              id,ninsert_this,mass_inserted_this,step,ninserted,massinserted);

    if (logfile)
      fprintf(logfile,"INFO: Particle insertion %s: inserted %d particle templates (mass %f) at step " BIGINT_FORMAT "\n - a total of %d particle templates (mass %f) inserted so far.\n",
              id,ninsert_this,mass_inserted_this,step,ninserted,massinserted);
  }
}

/* ---------------------------------------------------------------------- */

int FixInsert::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsert::init()
{
    int ntimestep = update->ntimestep;

    if (!atom->radius_flag || !atom->rmass_flag)
        error->fix_error(FLERR,this,"Fix insert requires atom attributes radius, rmass");
    if (domain->triclinic)
        error->fix_error(FLERR,this,"Cannot use with triclinic box");
    if (domain->dimension != 3)
        error->fix_error(FLERR,this,"Can use fix insert for 3d simulations only");
    //NP gravity fix is checked in derived class if needed

    fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere", 0));
    if(!fix_multisphere) multisphere = NULL;
    else multisphere = &fix_multisphere->data();

    // in case of new fix insert in a restarted simulation, have to add current time-step
    if(next_reneighbor > 0 && next_reneighbor < ntimestep)
        error->fix_error(FLERR,this,"'start' step can not be before current step");

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

int FixInsert::min_type() const
{
    return type_min;
}

/* ---------------------------------------------------------------------- */

int FixInsert::max_type() const
{
    return type_max;
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_rad(int type) const
{
    return fix_distribution->max_rad(type);
}

/* ---------------------------------------------------------------------- */

double FixInsert::min_rad(int type) const
{
    return fix_distribution->min_rad(type);
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_r_bound() const
{
    return fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

double FixInsert::extend_cut_ghost() const
{
    if(!fix_multisphere)
        return 0.;
    //NP this is to extend ghost region
    return 2.*fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

int FixInsert::calc_ninsert_this()
{
  if(ninsert_per == 0.) error->fix_error(FLERR,this,"ninsert_per == 0.");

  // number of bodies to insert this timestep
  int ninsert_this = static_cast<int>(ninsert_per + randomAll->uniform());
  if (ninsert_exists && ninserted + ninsert_this > ninsert) ninsert_this = ninsert - ninserted;

  /*NL*/ // if (screen) fprintf(screen,"ninsert_per %f, ninsert_this %d\n",ninsert_per,ninsert_this);

  return ninsert_this;
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixInsert::pre_exchange()
{
  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 1\n");}

  int ninsert_this, ninsert_this_local; // global and local # bodies to insert this time-step

  // just return if should not be called on this timestep
  //NP second check insures that no double execution takes place in case of
  //NP multiple execution of pre_exchange() in case of load-balancing
  if (next_reneighbor != update->ntimestep || most_recent_ins_step == update->ntimestep) return;
  most_recent_ins_step = update->ntimestep;

  // things to be done before inserting new particles
  pre_insert();

  // number of particles to insert this timestep
  ninsert_this = calc_ninsert_this();
  /*NL*///if (screen) fprintf(screen,"ninsert_this %d\n",ninsert_this);

  // limit to max number of particles that shall be inserted
  // to avoid that max # may be slightly exceeded by random processes
  // in fix_distribution->randomize_list, set exact_number to 1
  if (ninsert_exists && ninserted + ninsert_this >= ninsert)
  {
      ninsert_this = ninsert - ninserted;
      if(ninsert_this < 0)
        ninsert_this = 0;
      exact_number = 1;
  }

  // distribute ninsert_this across processors
  ninsert_this_local = distribute_ninsert_this(ninsert_this);
  /*NL*///if (screen) fprintf(screen,"ninsert_this %d ninsert_this_local %d\n",ninsert_this,ninsert_this_local);

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 2\n");}

  // re-allocate list if necessary
  //NP list is local on processes
  if(ninsert_this_local > ninsert_this_max_local)
  {
      fix_distribution->random_init_list(ninsert_this_local);
      ninsert_this_max_local = ninsert_this_local;
  }

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 3\n");}

  // generate list of insertions
  // number of inserted particles can change if exact_number = 0
  //NP will generate the global pti list in fix_distribution
  //NP radius necessary for overlapcheck
  //NP generates a pti list that contain ninsert_this bodies
  //NP can be more spheres since in multi-sphere case

  //NP list is the same on all processors with global
  //NP need to make the list global as

  ninsert_this_local = fix_distribution->randomize_list(ninsert_this_local,groupbit,exact_number);
  /*NL*///if (screen) fprintf(screen,"fix %s inits distribution %s with ninsert_this_local %d\n",id,fix_distribution->id,ninsert_this_local);

  MPI_Sum_Scalar(ninsert_this_local,ninsert_this,world);

  if(ninsert_this == 0)
  {
      // warn if flowrate should be fulfilled
      if((nflowrate > 0. || massflowrate > 0.) && comm->me == 0)
        error->warning(FLERR,"Particle insertion: Inserting no particle - check particle insertion settings");

      // schedule next insertion
      if (insert_every && (!ninsert_exists || ninserted < ninsert))
        next_reneighbor += insert_every;
      else next_reneighbor = 0;

      return;
  }
  else if(ninsert_this < 0)
  {
      /*NL*/ //if (screen) fprintf(screen,"ninsert_this %d\n",ninsert_this);
      error->fix_error(FLERR,this,"Particle insertion: Internal error");
  }

  double min_subbox_extent;
  int min_dim;
  domain->min_subbox_extent(min_subbox_extent,min_dim);

  if(warn_boxentent && min_subbox_extent < 2.2 *max_r_bound())
  {
      char msg[200];
      sprintf(msg,"Particle insertion on proc %d: sub-domain too small to insert particles: \nMax. bounding "
                  "sphere diameter is %f sub-domain extent in %s direction is only %f ",
                  comm->me,2.*max_r_bound(),0==min_dim?"x":(1==min_dim?"y":"z"),min_subbox_extent);
      error->warning(FLERR,msg);
  }

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 4a\n");}

  // warn if max # insertions exceeded by random processes
  if (ninsert_exists && ninserted + ninsert_this > ninsert)
  {
      error->warning(FLERR,"INFO: Particle insertion: Number of particles to insert was slightly exceeded by random process");
  }

  // fill xnear array with particles to check overlap against
  //NP in case overlap needs to checked, count # particles in insertion volume
  // add particles in insertion volume to xnear list
  neighList.reset();

  if(check_ol_flag)
    load_xnear(ninsert_this_local);

  // insertion counters in this step
  int ninserted_this = 0, ninserted_spheres_this = 0;
  int ninserted_this_local = 0, ninserted_spheres_this_local = 0;
  double mass_inserted_this = 0.;
  double mass_inserted_this_local = 0.;

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 4b\n");}

  // randomize insertion positions and set v, omega
  // also performs overlap check via xnear if requested
  // returns # bodies and # spheres that could actually be inserted
  x_v_omega(ninsert_this_local,ninserted_this_local,ninserted_spheres_this_local,mass_inserted_this_local);

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 5\n");}

  // actual particle insertion

  fix_distribution->pre_insert(ninserted_this_local,fix_property,fix_property_value,property_index,fix_property_ivalue,property_iindex);

  //NP pti list is body list, so use ninserted_this as arg
  ninserted_spheres_this_local = fix_distribution->insert(ninserted_this_local);

  // warn if max # insertions exceeded by random processes
  if (ninsert_exists && ninserted + ninsert_this > ninsert)
  {
      error->warning(FLERR,"INFO: Particle insertion: Number of particles to insert was slightly exceeded by random process");
  }

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 6\n");}

  // set tag # of new particles beyond all previous atoms, reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  if (atom->tag_enable)
  {
    atom->tag_extend();
    atom->natoms += static_cast<double>(ninserted_spheres_this);
    if (atom->map_style)
    {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 7\n");}

  // give particle distributions the chance to do some wrap-up
  //NP multisphere things here if needed
  //NP setup inserted particles, overwrites particle velocity, which needs to be set to fulfill rigid body constraint
  //NP also sets molecule id
  fix_distribution->finalize_insertion();

  if (atom->molecular && atom->molecule_flag)
  {
    atom->mol_extend();
  }

  // give derived classes the chance to do some wrap-up
  //NP do this after distribution wrap up
  //NP since per-atom "body" is set for multisphere via particledistribution
  //NP finalize_insertion() in fix insert/stream needs "body" to be set
  finalize_insertion(ninserted_spheres_this_local);

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 8\n");}

  // tally stats
  MPI_Sum_Scalar(ninserted_this_local,ninserted_this,world);
  ninserted += ninserted_this;
  MPI_Sum_Scalar(mass_inserted_this_local,mass_inserted_this,world);
  massinserted += mass_inserted_this;
  print_stats_during(ninserted_this,mass_inserted_this);

  if(ninserted_this < ninsert_this && comm->me == 0)
      error->warning(FLERR,"Particle insertion: Less insertions than requested");

  // next timestep to insert
  if (insert_every && (!ninsert_exists || ninserted < ninsert)) next_reneighbor += insert_every;
  else next_reneighbor = 0;

  /*NL*/ if(LMP_DEBUGMODE_FIXINSERT) {MPI_Barrier(world); fprintf(LMP_DEBUG_OUT_FIXINSERT,"FixInsert::pre_exchange 9\n");}
}

/* ----------------------------------------------------------------------
   distribute insertions across processors
------------------------------------------------------------------------- */

int FixInsert::distribute_ninsert_this(int ninsert_this)
{
    if(ninsert_this == 0) return 0;

    int me, nprocs, ngap, ninsert_this_local, *ninsert_this_local_all;
    double fraction_local, fraction_local_all_sum, *fraction_local_all, *remainder, r, rsum;

    me = comm->me;
    nprocs = comm->nprocs;

    fraction_local = insertion_fraction();
    /*NL*/ //if (screen) fprintf(screen,"called distribute_ninsert_this, exact_number %d, fraction_local %f\n",exact_number,fraction_local);

    if(!exact_number)
        return static_cast<int>(fraction_local*static_cast<double>(ninsert_this) + random->uniform());

    // for exact_number==1, have to allgather to exactly match ninsert_this

    fraction_local_all = new double[nprocs];
    remainder = new double[nprocs];
    ninsert_this_local_all = new int[nprocs];

    // allgather local fractions
    MPI_Allgather(&fraction_local,1,MPI_DOUBLE,fraction_local_all,1,MPI_DOUBLE,world);

    // proc0 calculates ninsert_this_local for all processes
    //NP important to do this on one proc since random generation is involved
    //NP sum of fraction_local will not be 1 since random generators may have different states
    //NP so have to normalize here
    if(me == 0)
    {
        // remove fractions < 2% / nprocs
        // have to normalize so not all portions get cancelled away for higher proc counts
        // normalize fraction_local_all so sum across processors is 1

        double lower_thresh = 0.02 / static_cast<double>(nprocs);

        fraction_local_all_sum = 0.;
        for(int iproc = 0; iproc < nprocs; iproc++)
        {
            if(fraction_local_all[iproc] < lower_thresh)
                fraction_local_all[iproc] = 0.;
            fraction_local_all_sum += fraction_local_all[iproc];
        }

        if(fraction_local_all_sum == 0.)
            error->one(FLERR,"Internal error distributing particles to processes");

        for(int iproc = 0; iproc < nprocs; iproc++)
            fraction_local_all[iproc] /= fraction_local_all_sum;

        rsum = 0.;
        for(int iproc = 0; iproc < nprocs; iproc++)
        {
            ninsert_this_local_all[iproc] = static_cast<int>(fraction_local_all[iproc]*static_cast<double>(ninsert_this));
            remainder[iproc] = fraction_local_all[iproc]*static_cast<double>(ninsert_this) - ninsert_this_local_all[iproc];
            rsum += remainder[iproc];
            /*NL*/ //if (screen) fprintf(screen,"proc %d (fraction_local %f): ninsert_this_local %d, remainder %f \n",
            /*NL*/ //                iproc,fraction_local_all[iproc],ninsert_this_local_all[iproc],remainder[iproc]);
        }

        ngap = round(rsum);
        /*NL*/// if (screen) fprintf(screen,"ngap %d rsum %f\n",ngap,rsum);
        for(int i = 0; i < ngap; i++)
        {
            r = random->uniform() * static_cast<double>(ngap);
            int iproc = 0;
            rsum = remainder[iproc];

            while(iproc < (nprocs-1) && rsum < r)
            {
                iproc++;
                rsum += remainder[iproc];
            }
            ninsert_this_local_all[iproc]++;
        }
    }

    // Bcast the result
    MPI_Bcast(ninsert_this_local_all,nprocs, MPI_INT,0,world);
    ninsert_this_local = ninsert_this_local_all[me];

    /*NL*/ //if (screen) fprintf(screen,"proc %d: fraction_local %f ninsert_this_local %d ninsert_this %d\n",me,fraction_local,ninsert_this_local,ninsert_this);

    delete []fraction_local_all;
    delete []remainder;
    delete []ninsert_this_local_all;

    return ninsert_this_local;
}

/* ----------------------------------------------------------------------
   count # of particles that could overlap
   must loop local + ghost particles
------------------------------------------------------------------------- */

int FixInsert::count_nnear()
{
    int nall = atom->nlocal + atom->nghost;
    int ncount = 0;

    for(int i = 0; i < nall; i++)
        ncount += is_nearby(i);

    return ncount;
}

/* ----------------------------------------------------------------------
   fill neighbor list with nearby particles
------------------------------------------------------------------------- */

int FixInsert::load_xnear(int)
{
  // load up neighbor list with local and ghosts

  neighList.reset();
  if(maxrad <= 0.)
    return 0;

  BoundingBox bb = getBoundingBox();

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  neighList.set_obb_flag(check_obb_flag);
#endif

#ifdef LIGGGHTS_DEBUG
  printf("subdomain bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n", domain->sublo[0], domain->subhi[0], domain->sublo[1], domain->subhi[1], domain->sublo[2], domain->subhi[2]);
#endif

  if(neighList.setBoundingBox(bb, maxrad)) {
    double **x = atom->x;
    double *radius = atom->radius;
    int *type = atom->type;
    const int nall = atom->nlocal + atom->nghost;

    for (int i = 0; i < nall; ++i)
    {
      if (is_nearby(i))
      {
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(atom->superquadric_flag && check_obb_flag)
          neighList.insert_superquadric(x[i], radius[i], type[i], atom->quaternion[i], atom->shape[i], atom->blockiness[i]);
        else
          neighList.insert(x[i], radius[i], type[i]);
#else
        neighList.insert(x[i], radius[i], type[i]);
#endif
      }
    }
  }

  return neighList.count();
}

/* ----------------------------------------------------------------------
   generate random velocity based on random setting
------------------------------------------------------------------------- */

void FixInsert::generate_random_velocity(double * velocity) {
  switch(v_randomSetting) {
    case RANDOM_UNIFORM:
      velocity[0] = v_insert[0] + v_insertFluct[0] * 2.0 * (random->uniform()-0.50);
      velocity[1] = v_insert[1] + v_insertFluct[1] * 2.0 * (random->uniform()-0.50);
      velocity[2] = v_insert[2] + v_insertFluct[2] * 2.0 * (random->uniform()-0.50);
      break;

    case RANDOM_GAUSSIAN:
      velocity[0] = v_insert[0] + v_insertFluct[0] * random->gaussian();
      velocity[1] = v_insert[1] + v_insertFluct[1] * random->gaussian();
      velocity[2] = v_insert[2] + v_insertFluct[2] * random->gaussian();
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixInsert::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = static_cast<double>(random->state());
  list[n++] = static_cast<double>(ninserted);
  list[n++] = static_cast<double>(first_ins_step);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = massinserted;

  /*NL*/ //if (screen) fprintf(screen,"next_reneighbor %d ninserted %d\n",next_reneighbor,ninserted);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixInsert::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  bigint next_reneighbor_re;

  seed = static_cast<int> (list[n++]);
  ninserted = static_cast<int> (list[n++]);
  first_ins_step = static_cast<int> (list[n++]);
  next_reneighbor_re = static_cast<bigint> (list[n++]);
  massinserted = list[n++];

  randomAll->reset(seed);
  seed += comm->me;
  random->reset(seed);

  // in order to be able to continue pouring with increased number of particles
  // if insert was already finished in run to be restarted
  if(next_reneighbor_re != 0) next_reneighbor = next_reneighbor_re;

  /*NL*/// if (screen) fprintf(screen,"next_reneighbor " BIGINT_FORMAT ", next_reneighbor_re " BIGINT_FORMAT "  "
  /*NL*///       "ninserted %d ninsert %d\n",next_reneighbor,next_reneighbor_re,ninserted,ninsert);
}

/* ----------------------------------------------------------------------
   output
------------------------------------------------------------------------- */

double FixInsert::compute_vector(int index)
{
    if(index == 0) return static_cast<double>(ninserted);
    if(index == 1) return massinserted;
    return 0.0;
}
