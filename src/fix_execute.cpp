/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
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
  conditional(false),
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

  if (narg == 6 && strcmp(arg[5],"once")==0)
  {
    once = true;
    execution_step = nevery;
    nevery = 1;
  }
  else if (narg == 9 && strcmp(arg[5],"conditional") == 0) {
    conditional = true;
    int n = strlen(arg[6]) + 1;
    var = new char[n];
    strcpy(var,arg[6]);
    var_valid = atof(arg[7]);
    var_threshold = atof(arg[8]);
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
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExecute::end_of_step()
{
  if (conditional)
  {
    int ivar = input->variable->find(var);
    if (ivar < 0) error->fix_error(FLERR,this,"target variable not found");
    double value = input->variable->compute_equal(ivar);
    if (fabs(value-var_valid) > var_threshold) return;
  }

  if (!once || update->ntimestep == execution_step )
  input->one(string);
}
