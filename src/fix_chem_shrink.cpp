/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

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
   Thomas Lichtenegger (JKU Linz)
   M.Efe Kinaci (JKU Linz)
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_chem_shrink.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrink::FixChemShrink(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg)
{
    if (strncmp(style,"chem/shrink",14) == 0 && (!atom->radius_flag))
            error -> all (FLERR,"Fix chem/shrink needs particle radius");

    // defaults
    fix_concA_   =   NULL;
    fix_concC_   =   NULL;
    fix_changeOfA_   =   NULL;
    fix_changeOfC_   =   NULL;
    // fix_tpart_       =   NULL;
    // fix_rhogas_      =   NULL;
    // fix_reactionheat_    =   0;

    iarg_ = 3;
    if (narg < 14)
        error -> all (FLERR,"not enough arguments");

    // check and define species A
    if (strcmp(arg[iarg_++],"speciesA") != 0)
        error -> all (FLERR, "missing keyword speciesA");
    speciesA = new char [strlen(arg[iarg_])+1];
    strcpy(speciesA, arg[iarg_++]);
    if (speciesA == 0)
        error -> all (FLERR, "speciesA not defined");

    // check and define molecular mass of A
    if (strcmp(arg[iarg_++],"molMassA") != 0)
        error -> all (FLERR, "keyword molMassA for species A is missing");
    molMass_A_  =   atof(arg[iarg_++]);
    if (molMass_A_ < 1)
        error -> all (FLERR,"molMass species A is not defined");

    // check and define Species C
    if (strcmp(arg[iarg_++],"speciesC") != 0)
        error -> all (FLERR, "missing keyword speciesC");
    speciesC = new char [strlen(arg[iarg_])+1];
    strcpy(speciesC, arg[iarg_++]);
    if (speciesC == 0)
        error -> all (FLERR, "speciesC not defined");

    // check and define molecular mass of C
    if (strcmp(arg[iarg_++],"molMassC") != 0)
        error -> all (FLERR, "keyword molMassC for species C is missing");
    molMass_C_  =   atof(arg[iarg_++]);
    if (molMass_C_ < 1)
        error -> all (FLERR,"molMass species C is not defined");

    // define the molecular mass of solid molecule
    if (strcmp(arg[iarg_++],"molMassB") != 0)
        error -> all (FLERR, "keyword molMassB for species B is missing");
    molMass_B_  =   atof(arg[iarg_++]);
    if (molMass_B_ < 1)
        error -> all (FLERR,"molMass species solid is not defined");

    // define reaction rate coefficient
    if (strcmp(arg[iarg_++], "k") != 0)
        error -> all (FLERR,"keyword k for reaction rate is missing");
    k   =   atof(arg[iarg_++]);

    // flags for vector output
    scalar_flag =   1;
    global_freq =   1;
    extscalar   =   1;

}

/* ---------------------------------------------------------------------- */

FixChemShrink::~FixChemShrink()
{
    if(fix_changeOfA_)  delete  []fix_changeOfA_;
    if(fix_changeOfC_)  delete  []fix_changeOfC_;

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_concA_) modify  ->  delete_fix(speciesA);
    if(unfixflag && fix_concC_) modify  ->  delete_fix(speciesC);
    // if(unfixflag && fix_rhogas_)    modify  -> delete_fix("partRho");
    // if(unfixflag && fix_tgas_)  modify  ->  delete_fix("partTemp");

    if(unfixflag && fix_changeOfA_) modify  ->  delete_fix("changeOfSpeciessMass_A");
    if(unfixflag && fix_changeOfC_) modify  ->  delete_fix("changeOfSpeciesMass_C");


}

/* ---------------------------------------------------------------------- */

int FixChemShrink::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_create()
{
    // register concentration species A
    if (!fix_concA_)
    {
        const char* fixarg[9];
        fixarg[0]=speciesA;
        fixarg[2]="property/atom";
        fixarg[3]=speciesA;
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_concA_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register concentration species C
    if (!fix_concC_)
    {
        const char* fixarg[9];
        fixarg[0]=speciesC;
        fixarg[2]="property/atom";
        fixarg[3]=speciesC;
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_concC_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register change of species mass A
    if (!fix_changeOfA_)
    {
        const char* fixarg[9];
        fixarg[0]="changeOfSpeciessMass_A";
        fixarg[2]="property/atom";
        fixarg[3]="changeOfSpeciessMass_A";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register change of species mass C
    if (!fix_changeOfA_)
    {
        const char* fixarg[9];
        fixarg[0]="changeOfSpeciessMass_C";
        fixarg[2]="property/atom";
        fixarg[3]="changeOfSpeciessMass_C";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::init()
{
    if (!atom -> radius_flag || !atom -> density_flag)
        error -> all(FLERR,"Fix chem/shrink cna only be used with sphere atom style");

    // references
    fix_concA_            =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesA,"property/atom","scalar",0,0,style));
    fix_concC_          =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesC,"property/atom","scalar",0,0,style));
    fix_changeOfA_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("changeOfSpeciessMass_A","property/atom","scalar",0,0,style));
    fix_changeOfC_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("changeOfSpeciessMass_C","property/atom","scalar",0,0,style));

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    // radius_ -> atom -> radius;

    // forAll,i
    // reaction()


    // void reaction
}

