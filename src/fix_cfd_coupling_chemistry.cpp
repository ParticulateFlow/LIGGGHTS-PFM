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
#include "fix_cfd_coupling_chemistry.h"
#include "fix_property_atom.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */


FixCfdCouplingChemistry::FixCfdCouplingChemistry(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg)
    // use_Re_(false)
{
    // defaults
    fix_coupling        =  NULL;
    fix_tgas            =  NULL;
    fix_rhogas          =  NULL;
    fix_massfrac_       =  NULL;
    fix_masschange_     =  NULL;
    fix_reactionheat_   =  NULL;

    iarg_ = 3;
    num_species = 0;
    int n;
    char mod[30];
    n_species = 0;

   if (narg < iarg_ + 3)
        error -> all (FLERR,"Fix couple/cfd/chemistry: Wrong number of arguments");

   bool hasargs = true;
   while (iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"n_species") == 0)
        {
            n_species = 1;
            iarg_++;
            num_species = atoi(arg[iarg_]);
            hasargs = true;
            iarg_ ++;
        }
        else if (strcmp(arg[iarg_],"species_names") == 0)
        {
            if (n_species != 1)
                error -> fix_error (FLERR,this, "have to define keyword 'n_species' before 'species_names'");
            if (iarg_ + num_species > narg)
                error -> fix_error(FLERR,this, "Wrong number of arguments");
            if (num_species < 1)
                error -> fix_error(FLERR,this,"n_species > 0 is required");
            iarg_++;

            species_names_ = new char*[num_species];
            mod_spec_names_ = new char*[num_species];
            for (int i = 0; i < num_species; i++)
            {
                n = strlen(arg[iarg_+i]) + 1;
                species_names_[i] = new char [n];
                mod_spec_names_[i] = new char [n];
                strcpy(species_names_[i], arg[iarg_+i]);
                strcpy(mod,"Modified_");
                strcat(mod,species_names_[i]);
                strcpy(mod_spec_names_[i],mod);

                // Fix_cfd_coupling_chemistry_control (are the speciesdetermined correctly)
                if (comm->me == 0 && screen)
                {
                    fprintf(screen,"Fix_cfd_coupling_chemistry_control (are the speciesdetermined correctly) \n");
                    fprintf(screen,"num species: %i \n",num_species);
                    fprintf(screen,"modified species is: %s \n", mod_spec_names_[i]);
                }
                // ----------------------------------------------------------------------
            }

            fix_massfrac_ = new FixPropertyAtom*[num_species];
            fix_masschange_ =   new FixPropertyAtom*[num_species];
            iarg_++;
            hasargs = true;
        }
        /* else if
         * {
           if(strcmp(arg[iarg],"transfer_Re") == 0)
           {
            if(narg < iarg+2)
                 error->fix_error(FLERR,this,"not enough arguments for 'transfer_Re'");
             iarg++;
             if(strcmp(arg[iarg],"yes") == 0)
                 use_Re_ = true;
             else if(strcmp(arg[iarg],"no") == 0)
                 use_Re_ = false;
             else
                 error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_Re'");
             iarg++;
             hasargs = true;
         }
         else if (strcmp(this->style,"couple/cfd/chemistry") == 0)
         {
             error->fix_error(FLERR,this,"unknown keyword");
         }*/

    }
    // flags for vector output
    vector_flag = 1;
    size_vector = 3;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingChemistry::~FixCfdCouplingChemistry()
{
    for (int i=0;i<num_species;i++)
    {
        if (species_names_[i]) delete [] species_names_[i];
        if (mod_spec_names_[i]) delete [] mod_spec_names_[i];
    }

    delete [] species_names_;
    delete [] mod_spec_names_;

    if(fix_massfrac_)       delete []fix_massfrac_;
    if(fix_masschange_)     delete []fix_masschange_;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_tgas)    modify -> delete_fix("partTemp");
    if(unfixflag && fix_rhogas)  modify -> delete_fix("partRho");

    if(unfixflag && fix_reactionheat_)  modify -> delete_fix("reactionHeat");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingChemistry::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_create()
{
    // register tgas (partTemp)
    if (!fix_tgas)
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
        fix_tgas = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register rhogas (partRho)
    if (!fix_rhogas)
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
        fix_rhogas = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
        fixarg[5]="no";     // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_reactionheat_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    //register massfractions
    for (int i=0; i<num_species;i++)
    {
        if(!fix_massfrac_[i])
        {
            char **fixarg = new char*[9];
            fixarg[0]=species_names_[i];                //fixarg[0]="concentrations";
            fixarg[1]=(char*) "all";
            fixarg[2]=(char*) "property/atom";
            fixarg[3]=species_names_[i];                //fixarg[3]="concentrations";
            fixarg[4]=(char*) "scalar";                 // 1 vector per particle to be registered
            fixarg[5]=(char*) "no";                     // restart
            fixarg[6]=(char*) "no";                     // communicate ghost
            fixarg[7]=(char*) "no";                     // communicate rev
            fixarg[8]=(char*) "0.";
            modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }

        // register masschange/changeOfSpeciesMass
        if(!fix_masschange_[i])
        {
            const char* fixarg[9];
            fixarg[0]=mod_spec_names_[i];       // fixarg[0]="concentrations";
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=mod_spec_names_[i];       // fixarg[3]="concentrations";
            fixarg[4]="scalar";                 // 1 vector per particle to be registered
            fixarg[5]="no";                     // restart
            fixarg[6]="no";                     // communicate ghost
            fixarg[7]="no";                     // communicate rev
            fixarg[8]="0.";
            modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::init()
{
    // make sure there is only one fix of thfor(int i = 1; i <= n_FixMesh_; i++)is style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->fix_error(FLERR,this,"Fix couple/cfd/chemistry needs a fix of type couple/cfd");

    // values to come from OF
    // if(use_Re_) fix_coupling->add_pull_property("Re","scalar-atom");
    fix_coupling    ->  add_pull_property("partTemp","scalar-atom");
    fix_coupling    ->  add_pull_property("partRho","scalar-atom");
    for (int i=0; i<num_species; i++)
    {
        fix_coupling -> add_pull_property(species_names_[i],"scalar-atom");
    }

    //  values to be transfered to OF
    fix_coupling    ->  add_push_property("reactionHeat","scalar-atom");
    for (int i = 0; i<num_species; i++)
    {
        fix_coupling    ->  add_push_property(mod_spec_names_[i],"scalar-atom");
    }

    // reference to partTemp, partRho and reactionheat
    fix_tgas            =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas          =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));

    // get reference for concentrations and changeOfSpeciesMass
    for (int i=0; i<num_species;i++)
    {
        fix_massfrac_[i]    =   static_cast<FixPropertyAtom*>(modify->find_fix_property(species_names_[i],"property/atom","scalar",0,0,style));
        fix_masschange_[i]  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(mod_spec_names_[i],"property/atom","scalar",0,0,style));
    }

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::initial_integrate(bigint)
{
    // for all species, reaction heat
    // if current timestep - 1 == latestpush(species name)
    // reset fix_masschange_(species name)
    // -1 is needed because time step is advanced before this function is called
    int nlocal = atom -> nlocal;
    double *reactionHeat_ = fix_reactionheat_ -> vector_atom;

    for (int i = 0; i < num_species; i++)
    {
        for (int j=0;j<nlocal;j++)
        {    
            if (update -> ntimestep - 1 == fix_coupling -> latestpush(species_names_[i]))
            {
                fix_masschange_[i] = NULL;
            }
           /* if (fix_masschange_[i] == NULL)
            {
                reactionHeat_[j] = NULL;
            }*/
            reactionHeat_[i] = NULL;
        }
    }

}

void FixCfdCouplingChemistry::post_force(int)
{
    // for all species names i and reaction heat
    // if dc_->pushednow(i)
    // clear masschange(i),reactionheat
    // pushednow(i)=false
}

