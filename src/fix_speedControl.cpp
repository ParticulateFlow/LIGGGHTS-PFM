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

#include <string.h>
#include <stdlib.h>
#include "fix_speedControl.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixSpeedControl::FixSpeedControl(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix speedControl command");

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else if (strcmp(arg[3],"NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else if (strcmp(arg[4],"NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else if (strcmp(arg[5],"NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  if (xstyle == NONE && ystyle == NONE && zstyle == NONE) {
    error->fix_error(FLERR,this,"no set-velocity specified; all velocity components are left floating, remove this fix!");
  }

  // default controller gain
  K = 0.1;

  // optional args

  iregion = -1;
  idregion = NULL;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix speedControl command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix speedControl does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"gain") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix speedControl command");
      K = force->numeric(FLERR,arg[iarg+1]);
      if(K <= 0.0)
        error->fix_error(FLERR,this,"gain K must be in range ]0;inf[");
      iarg += 2;

    } else error->all(FLERR,"Illegal fix speedControl command");
  }

  maxatom = 0;
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

FixSpeedControl::~FixSpeedControl()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSpeedControl::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix speedControl does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix speedControl is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix speedControl does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix speedControl is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix speedControl does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix speedControl is invalid style");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix speedControl does not exist");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;


  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  //NP modified C.K.
  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce,maxatom,4,"speedControl:sforce");
  }


  // constant force

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;

        if (xstyle) f[i][0] -= K*(v[i][0] - xvalue);
        if (ystyle) f[i][1] -= K*(v[i][1] - yvalue);
        if (zstyle) f[i][2] -= K*(v[i][2] - zvalue);
      }

  // variable force, wrap with clear/add
  // wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM && sforce)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM && sforce)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM && sforce)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (iregion >= 0 &&
            !domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
          continue;


        if (xstyle == ATOM) f[i][0] += sforce[i][0];
        else if (xstyle) f[i][0] -= K*(v[i][0] - xvalue);
        if (ystyle == ATOM) f[i][1] += sforce[i][1];
        else if (ystyle) f[i][1] -= K*(v[i][1] - yvalue);
        if (zstyle == ATOM) f[i][2] += sforce[i][2];
        else if (zstyle) f[i][2] -= K*(v[i][2] - zvalue);
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSpeedControl::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSpeedControl::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}
