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
------------------------------------------------------------------------- */

#include "string.h"
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
    fix_coupling_    =  NULL;
    fix_tgas_        =  0;
    fix_rhogas_      =  0;
    fix_massfrac_    =  NULL;
    fix_masschange_  =  NULL;
    fix_reactionheat_=  0;

    iarg_ = 3;
    narg_ = narg;
    num_species_ = 0;
    int n;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        // checking in liggghts_run if arguments are set up correctly
        if(strcmp(arg[iarg_],"n_species") == 0)
        {
            if(iarg_ + 2 > narg)
                error -> fix_error (FLERR, this, "Fix couple/cfd/chemistry: Wrong number of arguments");
            num_species_ = atoi(arg[iarg_++]);
            if(num_species_ < 1)
                error -> fix_error(FLERR, this, "'n_species' > 0 is required");
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_],"species_names") == 0)
        {
            if (num_species_ == 0)
                error -> fix_error(FLERR, this, "have to define # of species at 'n_species' first");
            if (narg < iarg_+1+num_species_)
                error -> fix_error (FLERR, this, "not enough arguments/all species are not defined");

            species_names_ = new char* [num_species_];
            for (int i = 0; i < num_species_; i++)
            {
                n = strlen(arg[iarg_] + 1);
                species_names_[i] = new char[n];
                strcpy(species_names_[i], arg[iarg_+i]);
            }

            fix_massfrac_ = new FixPropertyAtom*[num_species_];
            for (int i=0; i<num_species_; i++)
            {
                int n_i = modify -> find_fix(arg[iarg_+i]);
                if (n_i == -1)
                    error -> fix_error (FLERR, this, "could not find fix you provided");
                fix_massfrac_[i-1] = static_cast<FixPropertyAtom*>(modify->fix[n_i]);
            }
            hasargs = true;
            iarg_ = 1 + num_species_;
        }
        /*else  if(strcmp(arg[iarg],"transfer_Re") == 0)
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
    for (int i=0;i<num_species_;i++) delete [] species_names_[i];
    delete [] species_names_;
    if(fix_massfrac_) delete []fix_massfrac_;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_create()
{
    // register tgas (partTemp)
    if (!fix_tgas_)
    {
        const char* fixarg[9];
        fixarg[0]="partTemp";
        fixarg[2]="property/atom";
        fixarg[3]="partTemp";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_tgas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register rhogas (partRho)
    if (!fix_rhogas_)
    {
        const char* fixarg[9];
        fixarg[0]="partRho";
        fixarg[2]="property/atom";
        fixarg[3]="partRho";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_rhogas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    //register massfractions // concentrations
    for( int i=0; i<num_species_;i++)
    {
        fix_massfrac_[i] = static_cast<FixPropertyAtom*>(modify->find_fix_property(species_names_[i],"property/atom","scalar",0,0,style,false));
        if(!fix_massfrac_[i])
        {
            const char* fixarg[9];
            fixarg[0]=species_names_[i]; //fixarg[0]="concentrations";
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=species_names_[i]; //fixarg[3]="concentrations";
            fixarg[4]="scalar"; // 1 vector per particle to be registered
            fixarg[5]="no";    // restart
            fixarg[6]="no";     // communicate ghost
            fixarg[7]="no";     // communicate rev
            fixarg[8]="1.";
            fix_massfrac_[i] = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }
    }

    // register reactionheat
    if (!fix_reactionheat_)
    {
        const char* fixarg[9];
        fixarg[0]="reactionHeat";
        fixarg[2]="property/atom";
        fixarg[3]="reactionHeat";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_reactionheat_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register masschange/changeOfSpeciesMass
    for( int i=0; i<num_species_;i++)
    {
        fix_masschange_[i] = static_cast<FixPropertyAtom*>(modify->find_fix_property(species_names_[i],"property/atom","scalar",0,0,style,false));
        if(!fix_masschange_[i])
        {
            const char* fixarg[9];
            fixarg[0]=species_names_[i]; //fixarg[0]="concentrations";
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=species_names_[i]; //fixarg[3]="concentrations";
            fixarg[4]="scalar"; // 1 vector per particle to be registered
            fixarg[5]="no";    // restart
            fixarg[6]="no";     // communicate ghost
            fixarg[7]="no";     // communicate rev
            fixarg[8]="1.";
            fix_masschange_[i] = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }
    }


}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_tgas_) modify -> delete_fix("partTemp");
    if(unfixflag && fix_rhogas_) modify -> delete_fix("partRho");
    if(unfixflag)
    {
        for (int i=0; i<num_species_; i++)
        {
          if (fix_massfrac_[i]) modify -> delete_fix(species_names_[i]);//("concentrations");        // should be deleted for every single species
          //if (fix_masschange_[i]) modify -> delete_fix("changeOfSpeciesMass"); // should be deleted for every single species
        }
    }
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingChemistry::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::init()
{
    // make sure there is only one fix of thfor(int i = 1; i <= n_FixMesh_; i++)is style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/chemistry needs a fix of type couple/cfd");

    //  values to be transfered to OF
    fix_coupling_->add_push_property("reactionheat","scalar-atom");
    for (int i = 0; i<num_species_; i++)
    {
        fix_coupling_->add_push_property("changeOfSpeciesMass","scalar-atom");
    }

    // values to come from OF
    // if(use_Re_) fix_coupling_->add_pull_property("Re","scalar-atom");
    fix_coupling_->add_pull_property("partTemp","scalar-atom");
    fix_coupling_->add_pull_property("partRho","scalar-atom");
    for (int i = 0; i<num_species_; i++)
    {
        fix_coupling_->add_pull_property(species_names_[i],"scalar-atom");
    }


}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_force(int)
{
  
// for all species names i and reaction heat
// if dc_->pushednow(i)
// clear masschange(i),reactionheat
// pushednow(i)=false
  
  
  
}

