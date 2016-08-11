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
#include "atom_vec.h"
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
    if (strncmp(style,"chem/shrink",14) == 0 && (!atom->radius_flag)||(!atom->rmass_flag))
            error -> all (FLERR,"Fix chem/shrink needs particle radius and mass");

    // defaults
    fix_concA_          =   NULL;
    fix_concC_          =   NULL;
    fix_changeOfA_      =   NULL;
    fix_changeOfC_      =   NULL;
    fix_rhogas_         =   NULL;
    fix_tgas_           =   NULL;
    fix_reactionheat_   =   0;


    iarg_ = 3;
    molMass_A_ = 0;
    molMass_C_ = 0;
    molMass_B_ = 0;
    k = 0;
    rmin = 0;
    int n = 16;
    char cha[30];
    spcA = 0;
    spcC = 0;
    pmass_ = NULL;

    //if (narg < 15)
    //    error -> all (FLERR,"not enough arguments");

    bool hasargs = true;
    while (iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            spcA = 1;
            iarg_++;
            speciesA = new char [strlen(arg[iarg_]) + 1];
            strcpy(speciesA,arg[iarg_]);
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            if (spcA != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesA' before 'molMassA'");
            if (strlen(speciesA) < 1)
                error -> fix_error(FLERR, this, "speciesA not defined");
            iarg_++;
            molMass_A_ = atof(arg[iarg_]);
            if (molMass_A_ < 1)
                error -> fix_error(FLERR, this, "molar mass of A is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            spcC = 1;
            iarg_++;
            speciesC = new char [strlen(arg[iarg_])+1];
            strcpy(speciesC,arg[iarg_]);
            iarg_ ++;
            hasargs = true;
        }
       else if (strcmp(arg[iarg_],"molMassC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            if (spcC != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesC' before 'molMassC'");
            if (strlen(speciesC) < 1)
                error -> fix_error(FLERR, this, "speciesC not defined");
            iarg_++;
            molMass_C_ = atof(arg[iarg_]);
            if (molMass_C_ < 1)
                error -> fix_error(FLERR, this, "molar mass of C is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassB") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_++;
            molMass_B_ = atof(arg[iarg_]);
            if (molMass_B_ < 1)
                error -> fix_error(FLERR, this, "molar mass of B is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"k") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_ ++;
            k = atof(arg[iarg_]);
            if (k == 0)
                error -> fix_error(FLERR, this, "k is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"rmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_++;
            rmin = atof(arg[iarg_]);
            if (rmin == 0)
                error -> fix_error(FLERR, this, "rmin is not defined");
            iarg_++;
            hasargs = true;
        }

        // print the arguments on screen
        if (comm -> me == 0 && screen)
        {
            fprintf(screen,"arg 3: %s \n",arg[3]); // check
            fprintf(screen,"arg 4: %s \n",arg[4]); // check
            fprintf(screen,"arg 5: %s \n",arg[5]); // check
            fprintf(screen,"arg 6: %f \n",molMass_A_); // check
            fprintf(screen,"arg 7: %s \n",arg[7]); // check
            fprintf(screen,"arg 8: %s \n",arg[8]); // check
            fprintf(screen,"arg 9: %s \n",arg[9]); // check
            fprintf(screen,"arg 10: %f \n",molMass_C_ ); // check
            fprintf(screen,"arg 11: %s \n",arg[11]); // check
            fprintf(screen,"arg 12: %f \n",molMass_B_); // check
            fprintf(screen,"arg 13: %s \n",arg[13]); // check
            fprintf(screen,"arg 14: %e \n",k); // check
            fprintf(screen,"arg 15: %s \n",arg[15]); // check
            fprintf(screen,"arg 16: %f \n",rmin); // check
        }
    }

    // define changed species mass A
    massA = new char [n];
    strcpy(cha,"Modified_");
    strcat(cha,speciesA);
    strcpy(massA,cha);

    // define changed species mass C
    massC = new char [n];
    strcpy(cha,"Modified_");
    strcat(cha,speciesC);
    strcpy(massC,cha);


    vector_flag =   1;
    size_vector =   3;
    global_freq =   1;
    extvector   =   1;
    atom -> rmass_flag = 1;

    if (comm -> me == 0 && screen)
    {
        fprintf(screen,"massA: %s \n", massA);
        fprintf(screen,"massC: %s \n", massC);
        fprintf(screen,"constructor successfully completed \n");
    }
}


/* ---------------------------------------------------------------------- */

FixChemShrink::~FixChemShrink()
{
    delete [] massA;
    delete [] massC;

    delete [] speciesA;
    delete [] speciesC;

    if (comm -> me == 0 && screen)
        fprintf(screen,"deconstruct successfully completed \n");
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if(fix_concA_)         modify  ->  delete_fix(speciesA);
        if(fix_concC_)         modify  ->  delete_fix(speciesC);
        if(fix_rhogas_)        modify  ->  delete_fix("partRho");
        if(fix_tgas_)          modify  ->  delete_fix("partTemp");
        if(fix_reactionheat_)  modify  ->  delete_fix("reactionHeat");
        if(fix_changeOfA_)     modify  ->  delete_fix(massA);
        if(fix_changeOfC_)     modify  ->  delete_fix(massC);

        if (comm -> me == 0 && screen)
            fprintf(screen,"pre-delete successfully completed \n");
    }
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
        fixarg[5]= "no";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        fix_concA_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
        fixarg[5]= "no";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        fix_concC_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
        fixarg[5]= "no";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        fix_changeOfA_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
        fixarg[5]= "no";    // restart
        fixarg[6]= "no";     // communicate ghost
        fixarg[7]= "no";     // communicate rev
        fixarg[8]= "0.";
        fix_changeOfC_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register rhogas (partRho)
    if (!fix_rhogas_)
    {
        const char* fixarg[9];
        fixarg[0]="partRho";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partRho";
        fixarg[4]="scalar";              // 1 vector per particle to be registered
        fixarg[5]="no";                  // restart yes
        fixarg[6]="yes";                 // communicate ghost no
        fixarg[7]="no";                  // communicate rev
        fixarg[8]="0.";
        fix_rhogas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register tgas (partTemp)
    if (!fix_tgas_)
    {
        const char* fixarg[9];
        fixarg[0]="partTemp";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partTemp";
        fixarg[4]="scalar";             // 1 vector per particle to be registered
        fixarg[5]="yes";                // restart yes
        fixarg[6]="no";                 // communicate ghost no
        fixarg[7]="no";                 // communicate rev
        fixarg[8]="0.";
        fix_tgas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register reactionheat
    if (!fix_reactionheat_)
    {
        const char* fixarg[9];
        fixarg[0]="reactionHeat";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="reactionHeat";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_reactionheat_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    if (comm -> me == 0 && screen)
        fprintf(screen,"post-create successfully completed \n");

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

    if (comm -> me == 0 && screen)
            fprintf(screen,"updatePtrs successfully \n");

}

/* ---------------------------------------------------------------------- */

void FixChemShrink::reaction()
    {
        updatePtrs();
        int nlocal  =   atom -> nlocal;
        // int natoms  =   atom -> natoms;
        // int nall = nlocal + atom -> nghost;
        double TimeStep = update -> dt;

        for (int i = 0; i<nlocal; i++)
        {
            if (screen)
            {
                fprintf(screen, " double k : %e \n", k);
                fprintf(screen, " double rhogas_[i] : %f \n", rhogas_[i]);
                fprintf(screen, " double concA : %f  \n", concA_[i]);
                fprintf(screen, " double radius_[i] : %f \n", radius_[i]);
                fprintf(screen, " timestep:%f \n", TimeStep);
            }
            // double dA   = -k*fix_rhogas_*fix_concA_*partSurfArea(radius_);
            double dA   =   -k*rhogas_[i]*concA_[i]*partSurfArea(radius_[i])*TimeStep;
            double dB   =    dA*(molMass_B_/molMass_A_);
            double dC   =   -dA*(molMass_C_/molMass_A_);

            // rate of change of the total number of moles of A
            // fix_changeOfA_  +=  dA;
            changeOfA_[i]       +=  dA;

            // rate of change of the total number of moles of C
            // fix_changeOfC_  +=  dC;
            changeOfC_[i]       +=  dC;

            // mass of particle subtracted by mass change of particle (delta m)
            pmass_[i]           +=  dB;

            // change of radius of particle -assumption: density of particle is constant
            radius_[i]         =   pow(0.75*pmass_[i]/(M_PI*pdensity_[i]),0.333333);

            if (screen)
            {
                fprintf(screen, " double dA: %e \n", dA);
                fprintf(screen, " double dC: %e \n", dC);
                fprintf(screen, " double changeOfA: %f  \n", changeOfA_[i]);
                fprintf(screen, " double changeOfC: %f \n", changeOfC_[i]);
                fprintf(screen, " double pmass: %f \n", pmass_[i]);
                fprintf(screen, " density: %f \n", pdensity_[i]);
                fprintf(screen, " testing reaction inside for loop \n");

            }
        }
        if (screen)
            fprintf(screen,"nlocal number is = %i \n", nlocal);
    }

/* ----------------- compute particle surface area ------------------------ */

double FixChemShrink::partSurfArea(double radius)
    {
        double A_p =   4*M_PI*radius*radius;
        return (A_p);
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

    // updatePtrs();

    if (comm -> me == 0 && screen)
        fprintf(screen,"init succesfully completed  \n");
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;
    int nlocal = atom -> nlocal;

    for (int i = 0; i < nlocal; i++)
    {
        if (radius_[i] >= rmin)
        {
            reaction();
        }
        else if (radius_[i] < rmin)
        {
            dlist[i] = radius_[i] < rmin;
            delete_atoms();
        }
    }

    if (comm -> me == 0 && screen)
        fprintf(screen,"post_force succesfully completed \n");
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::delete_atoms()
{
    AtomVec *avec = atom->avec;
    int nlocal = atom -> nlocal;

    int i = 0;
    while (i < nlocal) {
        if (dlist[i]) {
        avec->copy(nlocal-1,i,1);
        dlist[i] = dlist[nlocal-1];
        nlocal--;
        } else i++;
      }

    atom->nlocal = nlocal;
    memory->destroy(dlist);
    int comp_flag = 0;

    if (atom->molecular == 0 && comp_flag)
    {
        int *tag = atom -> tag;
        for (i = 0; i < nlocal; i++) tag[i] = 0;
        atom->tag_extend();
    }

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
}
