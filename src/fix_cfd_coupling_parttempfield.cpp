/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particulate Flow Modelling
   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "fix_cfd_coupling_parttempfield.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingPartTempField::FixCfdCouplingPartTempField(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_temp = NULL;

    int iarg = 3;

    if(narg < iarg + 2) error->all(FLERR,"Fix couple/cfd/parttempfield: Wrong number of arguments");
    if(strcmp(arg[iarg++],"T0") != 0) error->all(FLERR,"Fix couple/cfd/parttempfield: Expecting keyword 'T0'");
    T0 = atof(arg[iarg++]);

    if(T0 < 0.) error->all(FLERR,"Fix couple/cfd/parttempfield: T0 must be >= 0");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingPartTempField::~FixCfdCouplingPartTempField()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingPartTempField::pre_delete(bool unfixflag)
{
    if(fix_temp) modify->delete_fix("Temp");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingPartTempField::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingPartTempField::post_create()
{
  fix_temp=static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  if(!fix_temp)
  {
        const char* fixarg[9];
        fixarg[0]="Temp";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Temp";
        fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
        fixarg[5]="yes";    //NP restart yes
        fixarg[6]="yes";    //NP communicate ghost no
        fixarg[7]="no";    //NP communicate rev yes
        char arg8[30];
        sprintf(arg8,"%e",T0);
        fixarg[8]=arg8;

        modify->add_fix(9,const_cast<char**>(fixarg));
        fix_temp=static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0,style));
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingPartTempField::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    //values to come from OF
    fix_coupling->add_pull_property("Temp","scalar-atom");
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingPartTempField::post_force(int)
{
}
