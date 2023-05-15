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
#include <stdlib.h>
#include <string.h>
#include "fix_scalar_transport_equation.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "respa.h"
#include "properties.h"
#include "pair_gran.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::FixScalarTransportEquation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{

  if(strcmp(arg[2],"transportequation/scalar"))
    return;

  int iarg = 3;

  if (narg < 15)
    error->fix_error(FLERR,this,"not enough arguments");

  if(strcmp(arg[iarg++],"equation_id"))
    error->fix_error(FLERR,this,"expecting keyword 'equation_id'");
  equation_id = new char[strlen(arg[iarg])+1];
  strcpy(equation_id,arg[iarg++]);

  if(strcmp(arg[iarg++],"quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'quantity'");
  quantity_name = new char[strlen(arg[iarg])+1];
  strcpy(quantity_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"default_value"))
    error->fix_error(FLERR,this,"expecting keyword 'default_value'");
  quantity_0 = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"flux_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'flux_quantity'");
  flux_name = new char[strlen(arg[iarg])+1];
  strcpy(flux_name,arg[iarg++]);

  if(strcmp(arg[iarg++],"source_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'source_quantity'");
  source_name = new char[strlen(arg[iarg])+1];
  strcpy(source_name,arg[iarg++]);

  capacity_name = NULL;
  capacity_flag = 0;

  if(strcmp(arg[iarg++],"capacity_quantity"))
    error->fix_error(FLERR,this,"expecting keyword 'capacity_quantity'");
  if(strcmp(arg[iarg],"none"))
  {
      capacity_flag = 1;
      capacity_name = new char[strlen(arg[iarg])+1];
      strcpy(capacity_name,arg[iarg++]);
  }

  if(narg >= 17 && strcmp(arg[iarg++],"max_change") == 0)
  {
    max_change = atof(arg[iarg++]);
    limit_change = true;
  }
  else
  {
    limit_change = false;
  }

  fix_quantity = fix_flux = fix_source = NULL; fix_capacity = NULL; fix_capacity_per_atom = NULL;
  capacity = NULL;

  capacity_per_atom = false;
  time_dependent_capacity = false;
  int_flag = true;

  nevery_  = 1;
  performedIntegrationLastStep_ = true; //ensure flux is reset at the very first time step

  peratom_flag = 1;              //NP 0/1 if per-atom data is stored
  size_peratom_cols = 0;         //NP 0 = scalar, N = columns in peratom array
  peratom_freq = 1;

  scalar_flag = 1; //NP total thermal energy computed
  global_freq = 1; //NP available always
}

/* ---------------------------------------------------------------------- */

FixScalarTransportEquation::~FixScalarTransportEquation()
{
    delete []quantity_name;
    delete []flux_name;
    delete []source_name;
    delete []capacity_name;
    delete []equation_id;

    if(capacity) delete []capacity;
    //NP could delete fixes with no callbacks here since this fix has no callbacks
}

void FixScalarTransportEquation::pre_delete(bool unfixflag)
{
    //unregister property/atom fixes
    if(unfixflag)
    {
        if (fix_quantity) modify->delete_fix(quantity_name);
        if (fix_flux) modify->delete_fix(flux_name);
        if (fix_source) modify->delete_fix(source_name);
        if (fix_capacity_per_atom) modify->delete_fix(capacity_name);
    }
}

/* ---------------------------------------------------------------------- */

int FixScalarTransportEquation::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::updatePtrs()
{
  quantity = fix_quantity->vector_atom;
  flux = fix_flux->vector_atom;
  source = fix_source->vector_atom;

  vector_atom = quantity; //NP so that it is possible to access it from the input script
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::post_create()
{
  const char *fixarg[9];

  if (fix_quantity==NULL) {
    //register Temp as property/atom
    fixarg[0]=quantity_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=quantity_name;
    fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
    fixarg[5]="yes";    //NP restart yes
    fixarg[6]="yes";    //NP communicate ghost yes
    fixarg[7]="no";    //NP communicate rev no
    char arg8[30];
    sprintf(arg8,"%e",quantity_0);
    fixarg[8]=arg8;
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_quantity=static_cast<FixPropertyAtom*>(modify->find_fix_property(quantity_name,"property/atom","scalar",0,0,style));
  }

  if (fix_flux==NULL){
    //register heatFlux as property/atom
    fixarg[0]=flux_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=flux_name;
    fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
    fixarg[5]="yes";    //NP restart yes
    fixarg[6]="no";    //NP communicate ghost no
    fixarg[7]="yes";    //NP communicate rev yes
    fixarg[8]="0.";     //NP take 0 as default flux
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_flux=static_cast<FixPropertyAtom*>(modify->find_fix_property(flux_name,"property/atom","scalar",0,0,style));
  }

  if (fix_source==NULL){
    //register heatSource as property/atom
    fixarg[0]=source_name;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=source_name;
    fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
    fixarg[5]="yes";    //NP restart yes
    fixarg[6]="yes";    //NP communicate ghost yes
    fixarg[7]="no";    //NP communicate rev no
    fixarg[8]="0.";     //NP take 0 as default source
    modify->add_fix(9,const_cast<char**>(fixarg));
    fix_source=static_cast<FixPropertyAtom*>(modify->find_fix_property(source_name,"property/atom","scalar",0,0,style));
  }

  //NP Get pointer to all the fixes (also those that have the material properties)
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

double* FixScalarTransportEquation::get_capacity()
{
    return capacity;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::init()
{
  if (!atom->rmass_flag) error->all(FLERR,"Please use an atom style that defines per-particle mass for fix transportequation/scalar");

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if(capacity_flag)
  {
      if(!force->pair_match("gran", 0))
        error->fix_error(FLERR,this,"requires a granular pair style when used with capacityflag");
      PairGran* pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));

      int max_type = pair_gran->get_properties()->max_type();

      if(capacity) delete []capacity;
      capacity = new double[max_type+1];

      fix_capacity = static_cast<FixPropertyGlobal*>(modify->find_fix_property(capacity_name,"property/global","peratomtype",max_type,0,style,false));

      fix_capacity_per_atom = static_cast<FixPropertyAtom*>(modify->find_fix_property(capacity_name,"property/atom","scalar",0,0,style,false));
      if (fix_capacity_per_atom)
      {
          capacity_per_atom = true;
          time_dependent_capacity = fix_capacity_per_atom->store_old_time_values();
      }

      if (!fix_capacity && !fix_capacity_per_atom)
      {
          char errmsg[500];
          sprintf(errmsg,"Could not locate a fix/property storing value(s) for %s as requested by FixScalarTransportEquation.",capacity_name);
          error->all(FLERR,errmsg);
      }

      if (fix_capacity && fix_capacity_per_atom)
      {
          char errmsg[500];
          sprintf(errmsg,"Found both a fix/property/atom and a fix/property/global storing value(s) for %s as requested by FixScalarTransportEquation.",capacity_name);
          error->all(FLERR,errmsg);
      }

      //pre-calculate parameters for possible contact material combinations
      for(int i=1;i< max_type+1; i++)
          for(int j=1;j<max_type+1;j++)
              if (!capacity_per_atom)
              {
                  capacity[i] = fix_capacity->compute_vector(i-1);
              }
              else
              {
                  capacity[i] = 0.0;
              }
  }
}

/* ---------------------------------------------------------------------- */

int FixScalarTransportEquation::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"integrate") == 0) {
    if (narg < 2) error->fix_error(FLERR,this,"not enough arguments for fix_modify 'integrate'");

    if (strcmp(arg[1],"start") == 0) {
      int_flag = true;
    } else if (strcmp(arg[1],"stop") == 0) {
      int_flag = false;
    } else
      error->fix_error(FLERR,this,"wrong argument for fix_modify 'integrate'");
    return 2;
  }

  if (strcmp(arg[0],"every") == 0) {
    if (narg < 2) error->fix_error(FLERR,this,"not enough arguments for fix_modify 'every'");

    nevery_ = force->inumeric(FLERR,arg[1]);
//    printf("FixScalarTransportEquation: will perform update every %d time steps. \n", nevery_);

    return 1;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::initial_integrate(int vflag)
{
  //NP update because re-allocation might have taken place
  updatePtrs();

  //Skip in case there was NO Integration the last time step (to keep flux in mem)
  if(!performedIntegrationLastStep_)
        return;

  //reset heat flux
  //sources are not reset
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  {
       if (mask[i] & groupbit)
           flux[i]=0.;
  }

  /*NL*/ //if (screen) fprintf(screen,"executing FixScalarTransportEquation::initial_integrate, flux[0] %f [1] %f\n",flux[0],flux[1]);

  //NP do forward comm to send quantity (calculated from last time-step) to ghost particles
  fix_quantity->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::pre_force(int vflag)
{
    //NP do forward comm to send quantity (calculated from last time-step) to ghost particles
    //NP need to re-do this here on reneighboring steps since forward comm in
    //NP initial_integrate() was on 'old' ghost list
    if(neighbor->ago == 0)
        fix_quantity->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::final_integrate()
{
    double dt = update->dt;
    int nlocal = atom->nlocal;
    double *rmass = atom->rmass;
    int *type = atom->type;
    int *mask = atom->mask;

    // skip if integration turned off
    if(!int_flag)
        return;

    // skip if integration not wanted at this timestep
    if (update->ntimestep % nevery_)
    {
        performedIntegrationLastStep_ = false;
        return;
    }

    //NP update because re-allocation might have taken place
    updatePtrs();

    //NP do forward comm to send sources (the ones added this time-step) to ghost particles
    fix_source->do_forward_comm();

    if(capacity_flag)
    {
        double capacity, capacity_prev, quantity_new, quantity_prev;
        for (int i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit)
            {
                if (!capacity_per_atom)
                {
                    capacity = fix_capacity->compute_vector(type[i]-1);
                    capacity_prev = capacity;
                }
                else
                {
                    capacity = fix_capacity_per_atom->vector_atom[i];
                    if (time_dependent_capacity)
                    {
                        capacity_prev = fix_capacity_per_atom->old_time_values()->vector_atom[i];
                    }
                    else
                    {
                        capacity_prev = capacity;
                    }
                }
                quantity_prev = quantity[i];
                if(fabs(capacity) > SMALL)
                {
                    quantity_new = (capacity_prev * quantity_prev + (flux[i] + source[i]*double(nevery_)) * dt / rmass[i]) / capacity;
                    if (limit_change)
                    {
                        if (quantity_new - quantity_prev > max_change) quantity_new = quantity_prev + max_change;
                        else if (quantity_new - quantity_prev < -max_change) quantity_new = quantity_prev - max_change;
                    }
                    quantity[i] = quantity_new;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit)
            {
                quantity[i] += (
                                  flux[i]
                                + source[i]*double(nevery_) //multiply source to account for missing steps
                             ) * dt ;
            }
        }
    }
    performedIntegrationLastStep_ = true;
}

/* ---------------------------------------------------------------------- */

void FixScalarTransportEquation::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  // outermost level - update
  // all other levels - nothing

  if (ilevel == nlevels_respa-1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

double FixScalarTransportEquation::compute_scalar()
{
    double *rmass = atom->rmass;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    double capacity;

    //NP update because re-allocation might have taken place
    updatePtrs();

    double quantity_sum = 0.;

    if(capacity_flag)
    {
        for (int i = 0; i < nlocal; i++)
        {
            if (!capacity_per_atom)
            {
                capacity = fix_capacity->compute_vector(type[i]-1);
            }
            else
            {
                capacity = fix_capacity_per_atom->vector_atom[i];
            }
            quantity_sum += capacity * rmass[i] * quantity[i];
            /*NL*///if (screen) fprintf(screen,"step %d, proc %d, i %d quantity %f\n",update->ntimestep, comm->me,i,quantity[i]);
        }
    }
    else
    {
        for (int i = 0; i < nlocal; i++)
        {
           quantity_sum += quantity[i];
        }
    }

    MPI_Sum_Scalar(quantity_sum,world);

    /*NL*///if (screen) fprintf(screen,"step %d, proc %d, nlocal %d quantity_sum %f\n",update->ntimestep, comm->me,nlocal,quantity_sum);
    /*NL*/ //if(nlocal) error->all(FLERR,"end");
    return quantity_sum;
}

/* ---------------------------------------------------------------------- */

bool FixScalarTransportEquation::match_equation_id(const char* id)
{
    if(strcmp(id,equation_id)) return false;
    return true;
}
