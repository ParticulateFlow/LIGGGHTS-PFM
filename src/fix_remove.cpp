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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_remove.h"
#include "atom.h"
#include "error.h"
#include "update.h"
#include "region.h"
#include "domain.h"
#include "atom_vec.h"
#include "comm.h"
#include "random_park.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{REMOVE_SHRINK,REMOVE_DELETE};

/*NL*/ #define DEBUG_COREX_MATERIALREMOVE false
/*NL*/ #define DEBUG__OUT_COREX_MATERIALREMOVE screen

/* ---------------------------------------------------------------------- */

FixRemove::FixRemove(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 9) error->fix_error(FLERR,this,"");

  int iarg = 3;
  if (strcmp(arg[iarg++],"nevery") != 0)
      error->fix_error(FLERR,this,"expecting 'nevery' keyword");
  nevery = atoi(arg[iarg++]);
  if (nevery <= 0) error->fix_error(FLERR,this,"");

  if (strcmp(arg[iarg++],"massrate") != 0)
      error->fix_error(FLERR,this,"expecting 'massrate' keyword");
  rate_remove = atof(arg[iarg++]);
  if (rate_remove <= 0.) error->fix_error(FLERR,this,"");

  if (strcmp(arg[iarg++],"style") != 0)
      error->fix_error(FLERR,this,"expecting 'style' keyword");
  if (strcmp(arg[iarg],"shrink") == 0)
  {
      style = REMOVE_SHRINK;
      rad_mass_vary_flag = 1;
      iarg++;
      if (strcmp(arg[iarg++],"delete_below") != 0)
          error->fix_error(FLERR,this,"expecting 'delete_below' keyword");
      delete_below = atof(arg[iarg++]);
  }
  else if (strcmp(arg[iarg++],"delete") == 0)
      style = REMOVE_DELETE;
  else error->fix_error(FLERR,this,"expecting 'shrink' or 'delete'");

  if (strcmp(arg[iarg++],"seed") != 0)
      error->fix_error(FLERR,this,"expecting 'seed' keyword");
  seed = atoi(arg[iarg++]);
  if(seed <= 0)
      error->fix_error(FLERR,this,"seed must be > 1");

  iregion = -1;
  type_remove = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
          error->fix_error(FLERR,this,"region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"atomtype") == 0){
      iarg++;
      type_remove = atoi(arg[iarg++]);
      if (type_remove <= 0)
          error->fix_error(FLERR,this,"");
    } else error->fix_error(FLERR,this,"unknown keyword");
  }

  /*NL*///rad_mass_vary_flag = 1;
  time_depend = 1;

  mass_removed = 0.;
  mass_to_remove = 0.;
  time_origin = update->ntimestep;

  force_reneighbor = 1;
  next_reneighbor = time_origin + nevery;

  restart_global = 1; //NP modified C.K.

  // random number generator, same for all procs

  random = new RanPark(lmp,seed);
}

/* ---------------------------------------------------------------------- */

FixRemove::~FixRemove()
{

}

/* ---------------------------------------------------------------------- */

int FixRemove::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRemove::init()
{
  /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::init 0\n");

  // error checks
  if (!atom->radius_flag)
    error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");
  if (!atom->rmass_flag)
    error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");

  dt = update->dt;

  /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::init end\n");
}

/* ---------------------------------------------------------------------- */

void FixRemove::pre_exchange()
{
    int time_now = update->ntimestep;
    if(time_now != next_reneighbor) return;

    inflag = new int[atom->nlocal];

    int *type = atom->type;
    int *mask = atom->mask;
    double **x = atom->x;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *density = atom->density;

    double mass_eligible_me = 0., mass_eligible;
    double mass_shrink_me = 0.;
    double mass_removed_this_me = 0., mass_removed_this;
    double mass_to_remove_me;
    double ratio_m,ratio_r,prob;
    int nremoved_this = 0;

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 0\n");

    next_reneighbor = time_now + nevery;
    mass_to_remove = (time_now - time_origin) * dt * rate_remove - mass_removed;

    if(comm->me == 0)            fprintf( screen,"Timestep %d, removing material, amount to remove this step %f\n",time_now,mass_to_remove);
    if(comm->me == 0 && logfile) fprintf(logfile,"Timestep %d, removing material, amount to remove this step %f\n",time_now,mass_to_remove);

    if(mass_to_remove <= 0.) return;

    nremoved_this_me = 0;

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 1\n");

    //NP count total eligible mass
    for (int i = 0; i < atom->nlocal; i++)
    {
       inflag[i] = (mask[i] & groupbit && (type_remove < 0 || type[i] == type_remove) && (iregion < 0 || domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])));
       if(inflag[i])
       {
           mass_eligible_me += rmass[i];
           if(style == REMOVE_SHRINK && radius[i] >= delete_below) mass_shrink_me += rmass[i];
       }
    }

    //NP get eligible mass from all procs
    MPI_Allreduce(&mass_eligible_me,&mass_eligible,1,MPI_DOUBLE,MPI_SUM,world);
    mass_to_remove_me = mass_eligible_me / mass_eligible * mass_to_remove;

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 2\n");

    if(mass_eligible == 0.)
    {
        if(mass_to_remove > 0. && comm->me == 0)
            error->warning(FLERR,"Fix remove requested to removed mass, but no eligible particles found");
        return;
    }

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 3\n");

    //NP delete all particles if eligible mass smaller that mass to remove
    if(mass_eligible < mass_to_remove)
    {
        /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 4\n");
        int i = 0;
        while (i < atom->nlocal)
        {
            if(inflag[i])
            {
                mass_removed_this_me += rmass[i];
                delete_particle(i);
            }
            else i++;
        }
        if(comm->me == 0)
            error->warning(FLERR,"Fix remove removed less mass than requested");
    }
    //NP shrink eligible particles
    else if(style == REMOVE_SHRINK)
    {
        /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 5\n");

        //NP first remove all particles that are too small
        int i = 0;
        while (i < atom->nlocal && mass_to_remove_me > 0.)
        {
            if(inflag[i])
            {
                if(radius[i] < delete_below)
                {
                    mass_removed_this_me += rmass[i];
                    mass_to_remove_me -= rmass[i];
                    delete_particle(i);
                }
                else i++;
            }
            else i++;
        }

        //NP then shrink the rest of the particles
        //NP ratio is shrinkage factor
        if(mass_shrink_me > 0. && mass_to_remove_me > 0.)
        {
            ratio_m = 1. - mass_to_remove_me / mass_shrink_me;
            ratio_r = pow(ratio_m,0.33333);

            i = 0;
            while (i < atom->nlocal && mass_to_remove_me > 0.)
            {
                if(inflag[i])
                {
                    mass_removed_this_me += (1.-ratio_m)*rmass[i];
                    mass_to_remove_me -= (1.-ratio_m)*rmass[i];
                    rmass[i] *= ratio_m;
                    radius[i] *= ratio_r;
                }
                i++;
            }
        }
    }
    //NP delete particles
    else if(style == REMOVE_DELETE)
    {
        /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 6\n");

        double mass_to_remove_me_init = mass_to_remove_me, rand;
        int i = 0;
        while (i < atom->nlocal && mass_to_remove_me > 0.)
        {
            //NP randomize which particle to delete
            //NP this avoids a bias: large particles are inserted first

            i = static_cast<int>(random->uniform()*static_cast<double>(atom->nlocal));
            if(i == atom->nlocal) i--;

            if(inflag[i])
            {
                mass_removed_this_me += rmass[i];
                mass_to_remove_me -= rmass[i];
                delete_particle(i);
            }
        }
    }

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 7\n");

    //NP gather info about removed particles
    MPI_Allreduce(&mass_removed_this_me,&mass_removed_this,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&nremoved_this_me,&nremoved_this,1,MPI_INT,MPI_SUM,world);
    mass_removed += mass_removed_this;

    if(comm->me == 0)
        fprintf(screen,"    Ammount actually removed %f (#particles totally removed %d)\n",mass_removed_this,nremoved_this);
    if(comm->me == 0 && logfile)
        fprintf(logfile,"    Ammount actually removed %f (#particles totally removed %d)\n",mass_removed_this,nremoved_this);

    //NP tags and maps
    int i;
    if (atom->molecular == 0) {
      int *tag = atom->tag;
      for (i = 0; i < atom->nlocal; i++) tag[i] = 0;
      atom->tag_extend();
    }

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange 8\n");

    if (atom->tag_enable) {
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }

    delete []inflag;

    /*NL*/ if(DEBUG_COREX_MATERIALREMOVE) fprintf(DEBUG__OUT_COREX_MATERIALREMOVE,"FixRemove::pre_exchange end\n");
}

inline void FixRemove::delete_particle(int i)
{
    /*NL*///fprintf(screen,"proc %d deleting particle with mass %f\n",comm->me,atom->rmass[i]);
    //NP fprintf(screen,"deleting particle %d, nlocal %d, \n",i,atom->nlocal);
    atom->avec->copy(atom->nlocal-1,i,1);
    inflag[i] = inflag[atom->nlocal-1]; //NP also switch inflag
    atom->nlocal--;
    nremoved_this_me++;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixRemove::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = static_cast<double>(random->state());
  list[n++] = static_cast<double>(time_origin);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = mass_removed;
  list[n++] = rate_remove;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixRemove::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  double rate_remove_re;

  seed = static_cast<int> (list[n++]);
  time_origin = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);
  mass_removed = list[n++];
  rate_remove_re = list[n++];
  if(rate_remove_re != rate_remove)
    error->fix_error(FLERR,this,"Can not restart simulation with different removal rate");

  rate_remove = rate_remove_re;
  random->reset(seed);
}

