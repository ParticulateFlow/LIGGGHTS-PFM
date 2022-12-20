/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2021- Eindhoven University of Technology

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Tim Nijssen (TU/e)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "atom_vec.h"
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
#include "math_const.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_fluidproperties.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingFluidproperties::FixCfdCouplingFluidproperties(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg),
    fix_fluid_density_(0),
    fix_fluid_viscosity_(0)
{
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingFluidproperties::~FixCfdCouplingFluidproperties()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingFluidproperties::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_fluid_density_) modify->delete_fix("fluidDensity");
        if (fix_fluid_viscosity_) modify->delete_fix("fluidViscosity");
    }
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingFluidproperties::setmask()
{
    int mask = 0;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingFluidproperties::post_create()
{
    // register fluid density
    fix_fluid_density_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("fluidDensity","property/atom","scalar",0,0,style,false));
    if(!fix_fluid_density_)
    {
        const char* fixarg[9];
        fixarg[0]="fluidDensity";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="fluidDensity";
        fixarg[4]="scalar"; // 1 value per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_fluid_density_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register fluid viscosity
    fix_fluid_viscosity_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("fluidViscosity","property/atom","scalar",0,0,style,false));
    if(!fix_fluid_viscosity_)
    {
        const char* fixarg[9];
        fixarg[0]="fluidViscosity";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="fluidViscosity";
        fixarg[4]="scalar"; // 1 value per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_fluid_viscosity_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingFluidproperties::init()
{
    // make sure there is only one fix of this style
    if (modify->n_fixes_style(style) != 1)
        error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if (!fix_coupling_)
        error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    // values to come from OF
    fix_coupling_->add_pull_property("fluidDensity","scalar-atom");
    fix_coupling_->add_pull_property("fluidViscosity","scalar-atom");
}

/* ---------------------------------------------------------------------- */
