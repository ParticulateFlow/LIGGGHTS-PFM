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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_remove.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "region.h"
#include "domain.h"
#include "modify.h"
#include "atom_vec.h"
#include "comm.h"
#include "random_park.h"
#include "fix_multisphere.h"
#include "multisphere_parallel.h"
#include "math_const.h"
#include "input.h"
#include "variable.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{REMOVE_SHRINK,REMOVE_DELETE};


/* ---------------------------------------------------------------------- */

FixRemove::FixRemove(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
  random_(0),
  type_remove_(-1),
  iregion_(-1),
  style_(-1),            // set below
  delete_below_(0.),     // set below
  rate_remove_(0.),      // set below
  m_remove_min_(1e-6),
  seed_(0),              // set below
  integrated_error_(true),
  variable_rate_(false),
  ivar_(-1),
  mass_removed_(0.),
  mass_to_remove_(0.),
  integrated_rate_(0.),
  time_origin_(update->ntimestep),
  time_last_(update->ntimestep),
  restart_read_(true),
  restart_write_(true),
  verbose_(false),
  compress_flag_(1),
  monitor_heat_(false),
  heat_removed_(0.),
  fix_temp_(NULL),
  fix_capacity_(NULL),
  fix_internal_energy_(NULL),
  internal_energy_(false),
  fix_ms_(0),
  ms_(0)
{
  if (narg < 9) error->fix_error(FLERR,this,"");

  int iarg = 3;
  if (strcmp(arg[iarg++],"nevery") != 0)
      error->fix_error(FLERR,this,"expecting 'nevery' keyword");
  nevery = atoi(arg[iarg++]);
  if (nevery <= 0) error->fix_error(FLERR,this,"");

  if (strcmp(arg[iarg],"massrate") != 0 && strcmp(arg[iarg],"massratevariable") != 0)
      error->fix_error(FLERR,this,"expecting 'massrate' or 'massratevariable' keywords");
  if (strcmp(arg[iarg],"massrate") == 0)
  {
      rate_remove_ = atof(arg[iarg+1]);
      if (rate_remove_ <= 0.) error->fix_error(FLERR,this,"");
  }
  else if (strcmp(arg[iarg],"massratevariable") == 0)
  {
      variable_rate_ = true;
      int n = strlen(arg[iarg+1]) + 1;
      rate_name_ = new char[n];
      strcpy(rate_name_,arg[iarg+1]);
      if(comm->me == 0 && screen) fprintf(screen,"fix_remove: using variable mass rate\n");
  }
  iarg++;
  iarg++;

  if (strcmp(arg[iarg++],"style") != 0)
      error->fix_error(FLERR,this,"expecting 'style' keyword");
  if (strcmp(arg[iarg],"shrink") == 0)
  {
      style_ = REMOVE_SHRINK;
      rad_mass_vary_flag = 1;
      iarg++;
      if (strcmp(arg[iarg++],"delete_below") != 0)
          error->fix_error(FLERR,this,"expecting 'delete_below' keyword");
      delete_below_ = atof(arg[iarg++]);
  }
  else if (strcmp(arg[iarg++],"delete") == 0)
      style_ = REMOVE_DELETE;
  else error->fix_error(FLERR,this,"expecting 'shrink' or 'delete'");

  if (strcmp(arg[iarg++],"seed") != 0)
      error->fix_error(FLERR,this,"expecting 'seed' keyword");
  seed_ = atoi(arg[iarg++]);
  if(seed_ <= 0)
      error->fix_error(FLERR,this,"seed must be > 1");


  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      iregion_ = domain->find_region(arg[iarg+1]);
      if (iregion_ == -1)
          error->fix_error(FLERR,this,"region ID does not exist");
      iarg += 2;
    } else if (strcmp(arg[iarg],"atomtype") == 0){
      iarg++;
      type_remove_ = atoi(arg[iarg++]);
      if (type_remove_ <= 0)
          error->fix_error(FLERR,this,"'atomtype' > 0 required");
   } else if(strcmp(arg[iarg],"restart_read") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'restart_read'");
      if(strcmp(arg[iarg+1],"no") == 0)
        restart_read_ = false;
      else if(strcmp(arg[iarg+1],"yes"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'restart_read'");
      iarg += 2;
   } else if(strcmp(arg[iarg],"restart_write") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'restart_write'");
      if(strcmp(arg[iarg+1],"no") == 0)
        restart_write_ = false;
      else if(strcmp(arg[iarg+1],"yes"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'restart'");
      iarg += 2;
    } else if(strcmp(arg[iarg],"verbose") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'verbose'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        verbose_ = true;
      else if(strcmp(arg[iarg+1],"no"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'verbose'");
      iarg += 2;
    } else if(strcmp(arg[iarg],"integrated_error") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'integrated_error'");
      if(strcmp(arg[iarg+1],"no") == 0)
        integrated_error_ = false;
      else if(strcmp(arg[iarg+1],"yes"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'integrated_error'");
      iarg += 2;
    } else if(strcmp(arg[iarg],"monitor_heat") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'monitor_heat'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        monitor_heat_ = true;
      else if(strcmp(arg[iarg+1],"no"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'monitor_heat_'");
      iarg += 2;
    } else if (strcmp(arg[iarg],"minmass") == 0){
      iarg++;
      m_remove_min_ = atof(arg[iarg++]);
    } else if(strcmp(arg[iarg],"compress") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"Illegal compress option");
      if(strcmp(arg[iarg+1],"yes") == 0)
        compress_flag_ = 1;
      else if(strcmp(arg[iarg+1],"no") == 0)
        compress_flag_ = 0;
      else
        error->fix_error(FLERR,this,"Illegal compress option, expecing 'yes' or 'no'");
      iarg += 2;
    } else error->fix_error(FLERR,this,"unknown keyword");
  }

  /*NL*///rad_mass_vary_flag = 1;
  time_depend = 1;
  force_reneighbor = 1;
  next_reneighbor = time_origin_ + nevery;

  vector_flag = 1;
  if (monitor_heat_)
  {
    size_vector = 2;
  }
  else
  {
    size_vector = 1;
  }

  if (restart_read_ || restart_write_)
  {
    restart_global = 1; //NP modified C.K.
  }
  else
  {
    restart_global = 0;
  }

  // random number generator, same for all procs

  random_ = new RanPark(lmp,seed_);
}

/* ---------------------------------------------------------------------- */

FixRemove::~FixRemove()
{
  delete random_;
  if (variable_rate_) delete [] rate_name_;
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
  // error checks
  if (!atom->radius_flag)
    error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");
  if (!atom->rmass_flag)
    error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");
  if (!atom->tag_enable || 0 == atom->map_style)
    error->fix_error(FLERR,this,"requires atom tags and an atom map");

  dt_ = update->dt;

  fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0));

  if(modify->n_fixes_style("multisphere") > 1)
    error->fix_error(FLERR,this,"does not support more than one fix multisphere.");

  if(fix_ms_)
  {
      ms_ = &(fix_ms_->data());
      if(REMOVE_SHRINK == style_)
          error->fix_error(FLERR,this,"does not support style 'shrink' with multisphere particles");

      //NP check order

      //NP this fix must be after fix multisphere since it this fix deletes
      //NP bodies via call-back (pre_neighbor; from fix multisphere)
      //NP atoms belonging to these bodies must be deleted by fix multisphere
      //NP in pre_exchange before the next deletions are performed by fix remove
      //NP in pre_exchange; otherwise some particles would be double-counted

      //NP also, this fix must register with fix multisphere
      if(modify->index_first_fix_of_style("multisphere") > modify->my_index(this))
          error->fix_error(FLERR,this,"this fix must be defined after fix multisphere");

      fix_ms_->add_remove_callback(this);
  }
  else
    ms_ = 0;


  if (variable_rate_) ivar_ = input->variable->find(rate_name_);

  if (monitor_heat_)
  {
      fix_temp_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp", "property/atom", "scalar", 0, 0, style));

      PairGran* pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
      int max_type = pair_gran->get_properties()->max_type();
      fix_capacity_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalCapacity","property/global","peratomtype",max_type,0,style,false));

      fix_internal_energy_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("internalEnergy","property/atom","scalar",0,0,style,false));
      if (fix_internal_energy_) internal_energy_ = true;

      if (!fix_capacity_ && !fix_internal_energy_)
      {
          char errmsg[500];
          sprintf(errmsg,"Could neither locate a fix/property storing value(s) for thermalCapacity nor one for the internal energy as requested by FixRemove.");
          error->all(FLERR,errmsg);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixRemove::pre_exchange()
{
    // only perform action if neigh build triggered by this fix

    int time_now = update->ntimestep;
    if(time_now != next_reneighbor) return;

    // update if variable mass rate
    if (variable_rate_)
    {
        if(comm->me == 0 && screen) fprintf(screen,"fix_remove: updating removal rate with variable %s\n",rate_name_);
        ivar_ = input->variable->find(rate_name_);
        rate_remove_ = input->variable->compute_equal(ivar_);
    }

    // clear containers

    atom_tags_eligible_.clear();
    body_tags_eligible_.clear();
    body_tags_delete_.clear();

    // list of particles eligible for deletion
    // list of ms bodies eligible for deletion

    // schedule next event and calc mass to remove this step

    next_reneighbor = time_now + nevery;
    if (integrated_error_)
    {
        integrated_rate_ += (time_now - time_last_) * dt_ * rate_remove_;
        mass_to_remove_ = integrated_rate_ - mass_removed_;
    }
    else
    {
        mass_to_remove_ = (time_now - time_last_) * dt_ * rate_remove_;
    }
    time_last_ = time_now;

    // print to logfile

    if(comm->me == 0)
    {
        if(verbose_ && screen)
            fprintf(screen,"Timestep %d, removing material, mass to remove this step %f\n",
                        time_now,mass_to_remove_);
        if(logfile)
            fprintf(logfile,"Timestep %d, removing material, mass to remove this step %f\n",
                        time_now,mass_to_remove_);
        if(logfile)
            fprintf(logfile,"fix_remove: integrated rate = %f, rate_remove = %f\n",integrated_rate_,rate_remove_);
    }

    // return if nothing to do
    // could e.g. be because very large particle was deleted last deletion step

    if(mass_to_remove_ <= m_remove_min_) return;

    // count eligible particles

    double mass_eligible_me = 0., mass_eligible = 0.;
    double mass_to_remove_me = 0., ratio_ms_to_remove_me = 0.;
    double mass_shrink_me = 0.;

    if(!count_eligible(mass_eligible_me,mass_eligible,mass_shrink_me,
                       mass_to_remove_me,ratio_ms_to_remove_me))
        return;

    double mass_removed_this_me = 0., mass_removed_this = 0.;
    double heat_removed_this_me = 0., heat_removed_this = 0.;
    int nremoved_this = 0, nremoved_this_me = 0;

    //NP delete all particles if eligible mass smaller that mass to remove
    if(mass_eligible <= mass_to_remove_)
        delete_all(mass_eligible_me,ratio_ms_to_remove_me,mass_removed_this_me,heat_removed_this_me,nremoved_this_me);
    //NP shrink eligible particles
    else if(REMOVE_SHRINK == style_)
        shrink(mass_to_remove_me,mass_shrink_me,mass_removed_this_me,heat_removed_this_me,nremoved_this_me);
    //NP delete particles
    else if(REMOVE_DELETE == style_ && !fix_ms_)
        delete_partial_particles(mass_to_remove_me,mass_removed_this_me,heat_removed_this_me,nremoved_this_me);
    else if(REMOVE_DELETE == style_)
        delete_partial_particles_bodies(mass_to_remove_me,mass_removed_this_me,
                                        heat_removed_this_me,nremoved_this_me,ratio_ms_to_remove_me);

    // gather info about removed particles
    MPI_Allreduce(&mass_removed_this_me,&mass_removed_this,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&nremoved_this_me,&nremoved_this,1,MPI_INT,MPI_SUM,world);
    mass_removed_ += mass_removed_this;

    if(monitor_heat_)
    {
        MPI_Allreduce(&heat_removed_this_me,&heat_removed_this,1,MPI_DOUBLE,MPI_SUM,world);
        heat_removed_ += heat_removed_this;
    }

    if(comm->me == 0)
    {
        if(verbose_ && screen)
            fprintf(screen,"    fix %s: this step: mass removed %f (of %f eligible), #particles removed %d; accumulated: mass removed: %f\n",
                    id,mass_removed_this,mass_eligible,nremoved_this,mass_removed_);
        if(logfile)
            fprintf(logfile,"    fix %s: this step: mass removed %f (of %f eligible), #particles removed %d; accumulated: mass removed: %f\n",
                    id,mass_removed_this,mass_eligible,nremoved_this,mass_removed_);
    }

    // update tags and map
    int i;

    // if non-molecular system and compress flag set,
    // reset atom tags to be contiguous

    if (atom->molecular == 0 && compress_flag_) {
      int *tag = atom->tag;
      for (i = 0; i < atom->nlocal; i++) tag[i] = 0;
      atom->tag_extend();
    }


    if (atom->tag_enable) {
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }
}

/* ----------------------------------------------------------------------
   count total eligible mass
------------------------------------------------------------------------- */

bool FixRemove::count_eligible(double &mass_eligible_me,double &mass_eligible,
                               double &mass_shrink_me,double &mass_to_remove_me,
                               double &ratio_ms_to_remove_me)
{
    int nlocal = atom->nlocal;
    int *type = atom->type;
    int *mask = atom->mask;
    int *tag = atom->tag;
    double **x = atom->x;
    double *rmass = atom->rmass;
    double *radius = atom->radius;
    double mass_ms_eligible_me = 0.;
    bool fits;

    // count total eligible mass for single particles

    for (int i = 0; i < nlocal; i++)
    {
       fits = (  (mask[i] & groupbit) &&
                 (!fix_ms_ || fix_ms_->belongs_to(i) < 0) &&
                 (type_remove_ < 0 || type[i] == type_remove_) &&
                 (iregion_ < 0 || domain->regions[iregion_]->match(x[i][0],x[i][1],x[i][2]))
              );

       if(fits)
       {
           atom_tags_eligible_.push_back(tag[i]);
           mass_eligible_me += rmass[i];
           if(REMOVE_SHRINK == style_  && radius[i] >= delete_below_)
              mass_shrink_me += rmass[i];
       }
    }

    // count total eligible mass for multisphere

    if(fix_ms_)
    {
        int nbody_local = ms_->n_body();
        double xcm[3];

        for (int ibody = 0; ibody < nbody_local; ibody++)
        {
            ms_->xcm(xcm,ibody);
            fits = (  (type_remove_ < 0 || ms_->atomtype(ibody) == type_remove_) &&
                      (iregion_ < 0 || domain->regions[iregion_]->match(xcm[0],xcm[1],xcm[2]))
                   );
            if(fits)
            {
                body_tags_eligible_.push_back(ms_->tag(ibody));
                mass_ms_eligible_me += ms_->mass(ibody);
                mass_eligible_me += ms_->mass(ibody);
            }
        }
    }

    // get eligible mass from all procs

    MPI_Allreduce(&mass_eligible_me,&mass_eligible,1,MPI_DOUBLE,MPI_SUM,world);

    if(mass_eligible_me > 0.)
        ratio_ms_to_remove_me = mass_ms_eligible_me / mass_eligible_me;
    else
        ratio_ms_to_remove_me = 0.;


    if(mass_eligible == 0.)
    {
        if(verbose_ && mass_to_remove_ > 0. && comm->me == 0)
            error->warning(FLERR,"Fix remove requested to removed mass, but no eligible particles found");
        return false;
    }

    mass_to_remove_me = mass_eligible_me / mass_eligible * mass_to_remove_;

    return true;
}

/* ----------------------------------------------------------------------
   delete all particles
------------------------------------------------------------------------- */

void FixRemove::delete_all(double mass_eligible_me,double ratio_ms_to_remove_me,
                           double &mass_removed_this_me,double &heat_removed_this_me,
                           int &nremoved_this_me)
{
    double *rmass = atom->rmass;
    int *type = atom->type;
    double *T = NULL;

    if(monitor_heat_)
    {
        T = fix_temp_->vector_atom;
    }

    for(size_t ilist = 0; ilist <  atom_tags_eligible_.size(); ilist++)
    {
        int i = atom->map(atom_tags_eligible_[ilist]);
        mass_removed_this_me += rmass[i];
        if(monitor_heat_)
        {
            if (internal_energy_)
            {
                heat_removed_this_me += fix_internal_energy_->vector_atom[i];
            }
            else
            {
                double Cp = fix_capacity_->compute_vector(type[i]-1);
                heat_removed_this_me += rmass[i]*T[i]*Cp;
            }
        }
        nremoved_this_me++;
        delete_particle(i);
    }
    atom_tags_eligible_.clear();

    if(fix_ms_ && ratio_ms_to_remove_me > 0.)
    {
        body_tags_delete_ = body_tags_eligible_;
        mass_removed_this_me += ratio_ms_to_remove_me*mass_eligible_me;
        nremoved_this_me += body_tags_delete_.size();
        //if (screen) fprintf(screen,"ratio_ms_to_remove_me %f mass_eligible_me %f\n",ratio_ms_to_remove_me,mass_eligible_me);
    }
    body_tags_eligible_.clear();

    if(verbose_ && 0 == comm->me)
        error->warning(FLERR,"Fix remove removed less mass than requested");
}

/* ----------------------------------------------------------------------
   shrink particles
------------------------------------------------------------------------- */

void FixRemove::shrink(double &mass_to_remove_me,double mass_shrink_me,
                       double &mass_removed_this_me,double &heat_removed_this_me,
                       int &nremoved_this_me)
{
    double ratio_m,ratio_r;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    int *type = atom->type;
    double *T = NULL;

    if(monitor_heat_)
    {
        T = fix_temp_->vector_atom;
    }

    // multisphere case NOT handled here, cannot shrink
    // error in init() function
    if(fix_ms_)
        error->fix_error(FLERR,this,"does not support multisphere");

    // first remove all particles that are too small
    size_t ilist = 0;
    while(ilist < atom_tags_eligible_.size())
    {
        int i = atom->map(atom_tags_eligible_[ilist]);

        if(radius[i] < delete_below_)
        {
            mass_removed_this_me += rmass[i];
            if(monitor_heat_)
            {
                if (internal_energy_)
                {
                    heat_removed_this_me += fix_internal_energy_->vector_atom[i];
                }
                else
                {
                    double Cp = fix_capacity_->compute_vector(type[i]-1);
                    heat_removed_this_me += rmass[i]*T[i]*Cp;
                }
            }
            nremoved_this_me++;
            mass_to_remove_me -= rmass[i];
            delete_particle(i);
            atom_tags_eligible_.erase(atom_tags_eligible_.begin()+ilist);
        }
        else ilist++;
    }

    // then shrink the rest of the particles
    // ratio is shrinkage factor
    if(mass_shrink_me > 0. && mass_to_remove_me > 0.)
    {
        ratio_m = 1. - mass_to_remove_me / mass_shrink_me;
        ratio_r = cbrt(ratio_m);

        for(size_t ilist = 0; ilist <  atom_tags_eligible_.size(); ilist++)
        {
            int i = atom->map(atom_tags_eligible_[ilist]);
            mass_removed_this_me += (1.-ratio_m)*rmass[i];
            if(monitor_heat_)
            {
                if (internal_energy_)
                {
                    heat_removed_this_me += (1.-ratio_m)*fix_internal_energy_->vector_atom[i];
                }
                else
                {
                    double Cp = fix_capacity_->compute_vector(type[i]-1);
                    heat_removed_this_me += (1.-ratio_m)*rmass[i]*T[i]*Cp;
                }
            }
            mass_to_remove_me -= (1.-ratio_m)*rmass[i];
            rmass[i] *= ratio_m;
            radius[i] *= ratio_r;
            if(mass_to_remove_me <= 0.) break;
        }
    }

    atom_tags_eligible_.clear();
}

/* ----------------------------------------------------------------------
   delete particles
------------------------------------------------------------------------- */

void FixRemove::delete_partial_particles(double &mass_to_remove_me,
                double &mass_removed_this_me,double &heat_removed_this_me,
                int &nremoved_this_me)
{
    double *rmass = atom->rmass;
    int *type = atom->type;
    double *T = NULL;

    if(monitor_heat_)
    {
        T = fix_temp_->vector_atom;
    }

    while (atom_tags_eligible_.size() > 0 && mass_to_remove_me > m_remove_min_)
    {
        // randomize which particle to delete
        // this avoids a bias: large particles are inserted first

        size_t ilist = static_cast<int>(random_->uniform()*static_cast<double>(atom_tags_eligible_.size()));
        if(ilist == atom_tags_eligible_.size())
            ilist--;
        int i = atom->map(atom_tags_eligible_[ilist]);

        // delete particle i
        mass_removed_this_me += rmass[i];
        if(monitor_heat_)
        {
            if (internal_energy_)
            {
                heat_removed_this_me += fix_internal_energy_->vector_atom[i];
            }
            else
            {
                double Cp = fix_capacity_->compute_vector(type[i]-1);
                heat_removed_this_me += rmass[i]*T[i]*Cp;
            }
        }
        nremoved_this_me++;
        mass_to_remove_me -= rmass[i];
        delete_particle(i);

        atom_tags_eligible_.erase(atom_tags_eligible_.begin()+ilist);
    }
}

/* ----------------------------------------------------------------------
   delete particles and bodies
------------------------------------------------------------------------- */

void FixRemove::delete_partial_particles_bodies(double &mass_to_remove_me,
                double &mass_removed_this_me,double &heat_removed_this_me,
                int &nremoved_this_me,double ratio_ms_to_remove_me)
{
    double *rmass = atom->rmass;
    bool ms;
    int i,ibody;
    int *type = atom->type;
    double *T = NULL;

    if(monitor_heat_)
    {
        T = fix_temp_->vector_atom;
    }

    while ((atom_tags_eligible_.size() > 0 || body_tags_eligible_.size() > 0) && mass_to_remove_me > 0.)
    {
        if(atom_tags_eligible_.size() == 0)
            ms = true;
        else if (body_tags_eligible_.size() == 0)
            ms = false;
        else
            ms = random_->uniform() <=  ratio_ms_to_remove_me;

        if(!ms)
        {
            // randomize which particle to delete
            // this avoids a bias: large particles are inserted first

            size_t ilist = static_cast<int>(random_->uniform()*static_cast<double>(atom_tags_eligible_.size()));
            if(ilist == atom_tags_eligible_.size())
                ilist--;
            i = atom->map(atom_tags_eligible_[ilist]);
            if(i == atom->nlocal) i--;

            // delete particle i
            mass_removed_this_me += rmass[i];
            if(monitor_heat_)
            {
                if (internal_energy_)
                {
                    heat_removed_this_me += fix_internal_energy_->vector_atom[i];
                }
                else
                {
                    double Cp = fix_capacity_->compute_vector(type[i]-1);
                    heat_removed_this_me += rmass[i]*T[i]*Cp;
                }
            }
            nremoved_this_me++;
            mass_to_remove_me -= rmass[i];
            delete_particle(i);
            atom_tags_eligible_.erase(atom_tags_eligible_.begin()+ilist);
        }
        else
        {
            size_t ilist = static_cast<int>(random_->uniform()*static_cast<double>(body_tags_eligible_.size()));
            if(ilist == body_tags_eligible_.size())
                ilist--;

            ibody = ms_->map(body_tags_eligible_[ilist]);
            mass_removed_this_me += ms_->mass(ibody);
            if(monitor_heat_)
            {
                error->fix_error(FLERR,this,"Monitoring of removed heat not implemented for multispheres yet.");
            }
            nremoved_this_me++;
            mass_to_remove_me -= ms_->mass(ibody);
            body_tags_delete_.push_back(body_tags_eligible_[ilist]);
            body_tags_eligible_.erase(body_tags_eligible_.begin()+ilist);
        }
    }
}

/* ----------------------------------------------------------------------
   delete local particle i
------------------------------------------------------------------------- */

inline void FixRemove::delete_particle(int i)
{
    // if (screen) fprintf(screen,"[%d] deleting particle %d, nlocal %d, \n",comm->me,i,atom->nlocal);
    atom->avec->copy(atom->nlocal-1,i,1);
    atom->nlocal--;
}
/* ----------------------------------------------------------------------
   delete bodies
------------------------------------------------------------------------- */

void FixRemove::delete_bodies()
{
    for(size_t ilist = 0; ilist <  body_tags_delete_.size(); ilist++)
    {
        int ibody = ms_->map(body_tags_delete_[ilist]);
        ms_->remove_body(ibody);
    }
    body_tags_delete_.clear();
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixRemove::write_restart(FILE *fp)
{
  if (!restart_write_) return;
  int n = 0;
  double list[5];
  list[n++] = static_cast<double>(random_->state());
  list[n++] = static_cast<double>(time_origin_);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = mass_removed_;
  list[n++] = rate_remove_;

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
  if (!restart_read_) return;
  int n = 0;
  double *list = (double *) buf;
  double rate_remove_re;

  seed_ = static_cast<int> (list[n++]);
  time_origin_ = static_cast<int> (list[n++]);
  next_reneighbor = static_cast<int> (list[n++]);
  mass_removed_ = list[n++];
  rate_remove_re = list[n++];
  if(rate_remove_re != rate_remove_)
    error->fix_error(FLERR,this,"Can not restart simulation with different removal rate");

  rate_remove_ = rate_remove_re;
  random_->reset(seed_);
}

/* ----------------------------------------------------------------------
   provide accumulated removed mass and optionally heat
------------------------------------------------------------------------- */

double FixRemove::compute_vector(int n)
{
  if (n==0)
  {
    return mass_removed_;
  }
  else if (n==1 && monitor_heat_)
  {
    return heat_removed_;
  }

  return -1.0;
}
