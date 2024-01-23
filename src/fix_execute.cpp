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
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "fix_execute.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include <math.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixExecute::FixExecute(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  execution_point(1),
  conditional(false),
  file(false),
  var(NULL),
  var_valid(0.0),
  var_threshold(0.0),
  once(false),
  execution_step(0)
{
  if (narg < 5) error->all(FLERR,"Illegal fix execute command");
  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix execute command");

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  string = new char[n];
  strcpy(string,arg[4]);


  int iarg = 5;

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;

    if(strcmp(arg[iarg],"once") == 0)
    {
      once = true;
      execution_step = nevery;
      nevery = 1;
      iarg++;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"conditional") == 0)
    {
      conditional = true;
      if(narg < iarg+4)
        error->fix_error(FLERR,this,"not enough arguments for 'conditional'");
      iarg++;
      int n = strlen(arg[iarg]) + 1;
      var = new char[n];
      strcpy(var,arg[iarg]);
      iarg++;
      var_valid = atof(arg[iarg]);
      iarg++;
      var_threshold = atof(arg[iarg]);
      iarg++;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"file") == 0)
    {
      file = true;
      iarg++;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"start_of_step") == 0)
    {
      execution_point = 0;
      iarg++;
      hasargs = true;
    }
    else 
    {
      error->fix_error(FLERR,this,"unknown keyword");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixExecute::~FixExecute()
{
  delete [] string;
}

/* ---------------------------------------------------------------------- */

int FixExecute::setmask()
{
  int mask = 0;
  if (execution_point == 0) mask |= INITIAL_INTEGRATE;
  if (execution_point == 1) mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExecute::initial_integrate(int)
{
  execution_command();
}

/* ---------------------------------------------------------------------- */

void FixExecute::end_of_step()
{
  execution_command();
}

/* ---------------------------------------------------------------------- */

void FixExecute::execution_command()
{
  if (conditional)
  {
    int ivar = input->variable->find(var);
    if (ivar < 0) error->fix_error(FLERR,this,"target variable not found");
    double value = input->variable->compute_equal(ivar);
    if (fabs(value-var_valid) > var_threshold) return;
  }

  if (!once || update->ntimestep == execution_step)
  {
    if (file)
    {
      input->file(string);
    }
    else
    {
      input->one(string);
    }
  }
}
