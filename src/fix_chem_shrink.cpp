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
    //if (strncmp(style,"chem/shrink",15) == 0 && (!atom->radius_flag)||(!atom->rmass_flag))
    //        error -> all (FLERR,"Fix chem/shrink needs particle radius and mass");

    // defaults
    fix_concA_          =   NULL;
    fix_concC_          =   NULL;
    fix_changeOfA_      =   NULL;
    fix_changeOfC_      =   NULL;
    fix_rhogas_         =   NULL;
    fix_tgas_           =   NULL;
    fix_reactionheat_   =   0;

    int n;
    char cha[30];

    iarg_ = 3;
    if (narg < 15)
        error -> all (FLERR,"not enough arguments");

    // check and define species A
    if (strcmp(arg[iarg_++],"speciesA") != 0)
        error -> all (FLERR, "missing keyword speciesA");

    printf("arg 3: %s \n",arg[3]); // check

    speciesA = new char [strlen(arg[iarg_])+1];
    strcpy(speciesA, arg[iarg_++]);

    printf("arg 4: %s \n",arg[4]); // check
    if (speciesA == 0)
        error -> all (FLERR, "speciesA not defined");

    // check and define molecular mass of A
    if (strcmp(arg[iarg_++],"molMassA") != 0)
        error -> all (FLERR, "keyword molMassA for species A is missing");
    printf("arg 5: %s \n",arg[5]); // check

    molMass_A_  =   atof(arg[iarg_++]);
    printf("arg 6: %f \n",molMass_A_); // check

    if (molMass_A_ == 0)
        error -> all (FLERR,"molMass species A is not defined");

    // check and define Species C
    if (strcmp(arg[iarg_++],"speciesC") != 0)
        error -> all (FLERR, "missing keyword speciesC");
    printf("arg 7: %s \n",arg[7]); // check

    speciesC = new char [strlen(arg[iarg_])+1];

    strcpy(speciesC, arg[iarg_++]);
    if (speciesC == 0)
        error -> all (FLERR, "speciesC not defined");
    printf("arg 8: %s \n",arg[8]); // check

    // check and define molecular mass of C
    if (strcmp(arg[iarg_++],"molMassC") != 0)
        error -> all (FLERR, "keyword molMassC for species C is missing");
    printf("arg 9: %s \n",arg[9]); // check

    molMass_C_  =   atof(arg[iarg_++]);
    printf("arg 10: %f \n",molMass_C_ ); // check

    if (molMass_C_ == 0)
        error -> all (FLERR,"molMass species C is not defined");

    // define the molecular mass of solid molecule
    if (strcmp(arg[iarg_++],"molMassB") != 0)
        error -> all (FLERR, "keyword molMassB for species B is missing");
    printf("arg 11: %s \n",arg[11]); // check

    molMass_B_  =   atof(arg[iarg_++]);
    printf("arg 12: %f \n",molMass_B_); // check

    if (molMass_B_ == 0)
        error -> all (FLERR,"molMass species solid is not defined");

    // define reaction rate coefficient
    if (strcmp(arg[iarg_++], "k") != 0)
        error -> all (FLERR,"keyword k for reaction rate is missing");
    printf("arg 13: %s \n",arg[13]); // check

    k   =   atof(arg[iarg_++]);
    printf("arg 14: %.4e \n",k); // check
    printf("arg 15: %s \n",arg[15]); // check

    // define changed species mass A
    massA = new char [n];
    n = strlen(arg[5]+1);
    strcpy(cha,"Changed_");
    strcat(cha,speciesA);
    strcpy(massA,cha);

    // define changed species mass A
    massC = new char [n];
    n = strlen(arg[9]+1);
    strcpy(cha,"Changed_");
    strcat(cha,speciesC);
    strcpy(massC,cha);

    // flags for vector output
    /*peratom_flag =  1;  //0/1 per-atom data is stored
    peratom_freq =  1;
    scalar_flag =   1;
    global_freq =   1;
    extscalar   =   1;*/

    vector_flag =   1;
    size_vector =   3;
    global_freq =   1;
    extvector   =   1;

}

/* ---------------------------------------------------------------------- */

FixChemShrink::~FixChemShrink()
{
    if (fix_concA_)     delete  []fix_concA_;
    if (fix_concC_)     delete  []fix_concC_;
    delete massA;
    delete massC;
    if(fix_changeOfA_)  delete  []fix_changeOfA_;
    if(fix_changeOfC_)  delete  []fix_changeOfC_;

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_concA_)         modify  ->  delete_fix(speciesA);
    if(unfixflag && fix_concC_)         modify  ->  delete_fix(speciesC);
    if(unfixflag && fix_rhogas_)        modify  ->  delete_fix("partRho");
    if(unfixflag && fix_tgas_)          modify  ->  delete_fix("partTemp");
    if(unfixflag && fix_reactionheat_)  modify  ->  delete_fix("reactionHeat");

    if(unfixflag && fix_changeOfA_) modify  ->  delete_fix(massA);
    if(unfixflag && fix_changeOfC_) modify  ->  delete_fix(massC);
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
   //register concentration species A
   if (!fix_concA_)
    {
        const char* fixarg[9];
        fixarg[0]= speciesA;
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= speciesA;
        fixarg[4]= "scalar"; // 1 scalar per particle to be registered
        fixarg[5]= "yes";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

   // register concentration species C
    if (!fix_concC_)
    {
        const char* fixarg[9];
        fixarg[0]= speciesC;
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= speciesC;
        fixarg[4]= "scalar"; // 1 scalar per particle to be registered
        fixarg[5]= "yes";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register change of species mass A
    if (!fix_changeOfA_)
    {
        const char* fixarg[9];
        fixarg[0]= massA;
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= massA;
        fixarg[4]= "scalar"; // 1 scalar per particle to be registered
        fixarg[5]= "yes";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register change of species mass C
    if (!fix_changeOfC_)
    {
        const char* fixarg[9];
        fixarg[0]= massC;
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= massC;
        fixarg[4]= "scalar"; // 1 scalar per particle to be registered
        fixarg[5]= "yes";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register rhogas (partRho)
    if (!fix_rhogas_)
    {
        const char* fixarg[9];
        fixarg[0]= "partRho";
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= "partRho";
        fixarg[4]= "scalar";              // 1 scalar per particle to be registered
        fixarg[5]= "no";                  // restart yes
        fixarg[6]= "yes";                 // communicate ghost no
        fixarg[7]= "no";                  // communicate rev
        fixarg[8]= "0.";
        fix_rhogas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register tgas (partTemp)
    if (!fix_tgas_)
    {
        const char* fixarg[9];
        fixarg[0]= "partTemp";
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= "partTemp";
        fixarg[4]= "scalar";             // 1 vector per particle to be registered
        fixarg[5]= "yes";                // restart yes
        fixarg[6]= "no";                 // communicate ghost no
        fixarg[7]= "no";                 // communicate rev
        fixarg[8]= "0.";
        fix_tgas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register reactionheat
    if (!fix_reactionheat_)
    {
        const char* fixarg[9];
        fixarg[0]= "reactionHeat";
        fixarg[1]= "all";
        fixarg[2]= "property/atom";
        fixarg[3]= "reactionHeat";
        fixarg[4]= "scalar"; // 1 vector per particle to be registered
        fixarg[5]= "yes";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        fix_reactionheat_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::updatePtrs()
{
    concA_      =   fix_concA_      -> vector_atom;
    concC_      =   fix_concC_      -> vector_atom;
    changeOfA_  =   fix_changeOfA_  -> vector_atom;
    changeOfC_  =   fix_changeOfC_  -> vector_atom;
    rhogas_     =   fix_rhogas_     -> vector_atom;
    //tgas_       =   fix_tgas_       -> vector_atom;

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::init()
{
   if (!atom -> radius_flag || !atom -> density_flag)
        error -> all(FLERR,"Fix chem/shrink can only be used with sphere atom style");

    // references
    fix_concA_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesA,"property/atom","scalar",0,0,style));
    fix_concC_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesC,"property/atom","scalar",0,0,style));
    fix_changeOfA_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massA,"property/atom","scalar",0,0,style));
    fix_changeOfC_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massC,"property/atom","scalar",0,0,style));

    fix_tgas_        =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_=   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));

    //updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;

    reaction();
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::reaction()
    {
        updatePtrs();
        int nlocal  =   atom->nlocal;

        for (int i = 0; i<nlocal; i++)
        {
            // double dA   = -k*fix_rhogas_*fix_concA_*partSurfArea(radius_);
            double dA   =   -k*rhogas_[i]*concA_[i]*partSurfArea(radius_[i]);
            double dC   =   -dA*(molMass_C_/molMass_A_);

            // rate of change of the total number of moles of A
            // fix_changeOfA_  +=  dA;
            changeOfA_[i]       +=  dA;

            // rate of change of the total number of moles of C
            // fix_changeOfC_  +=  dC;
            changeOfC_[i]       +=  dC;

            // rate of change of the total number of moles of B
            pmass_[i]          +=  dA*(molMass_B_/molMass_A_);

            // change of radius of particle -assumption: density of particle is constant
            radius_[i]         =   pow(0.75*pmass_[i]/(M_PI*pdensity_[i]),0.333333);
        }
     }

/* ----------------- compute particle surface area ------------------------ */

double FixChemShrink::partSurfArea(double radius)
    {
        double A_p =   4*M_PI*radius*radius;
        return (A_p);
    }

