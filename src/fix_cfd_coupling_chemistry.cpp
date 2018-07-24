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
#include "fix_property_global.h"
#include "group.h"
#include "pair.h"
#include "properties.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingChemistry::FixCfdCouplingChemistry(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg),
  fix_coupling_(0),
  fix_tgas_(0),
  fix_rhogas_(0),
  fix_masschange_(0),
  fix_reactionheat_(0),
  fix_totalmole_(0),
  fix_diffusionCoeff_(0),
  fix_nufField_(0),
  fix_partReynolds_(0),
  fix_molarfraction_(0),
  fix_partPressure_(0),
  use_reactant_(false)
{
    num_species = 0;
    n_species = 0;
    n_diff = 0;
    iarg_ = 3;
    mod_spec_names_ = diffusant_names_ = molarfraction_names = NULL;
    num_diffusant = 0;

   if (narg < iarg_ + 3)
        error -> all (FLERR,"Fix couple/cfd/chemistry: Wrong number of arguments");

   bool hasargs = true;

   while (iarg_ < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg_],"n_species") == 0)
        {
            n_species = 1;
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR,this,"Wrong number of arguments");
            iarg_++;
            num_species = atoi(arg[iarg_]);
            if (num_species < 1)
                error -> fix_error(FLERR,this,"n_species > 0 is required");
            hasargs = true;
            iarg_ ++;
        }
        else if (strcmp(arg[iarg_],"species_names") == 0)
        {
            if (n_species != 1)
                error -> fix_error (FLERR,this, "have to define n_species before 'species_names'");
            if (iarg_ + num_species > narg)
                error -> fix_error(FLERR,this, "Wrong number of arguments");
            species_names_ = new char*[num_species];
            iarg_++;
            for (int i = 0; i < num_species; ++i)
            {
                species_names_[i] = new char [strlen(arg[iarg_])+1];
                strcpy(species_names_[i], arg[iarg_]);
                iarg_ += 1;

                /*if(screen)
                {
                    fprintf(screen,"species_names in cfd/coupling/chemistry: %s \n",species_names_[i]);
                }*/
            }
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"n_diff") == 0)
        {
            n_diff = 1;
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR,this,"Wrong number of arguments");
            iarg_++;
            num_diffusant = atoi(arg[iarg_]);
            if (num_diffusant < 1)
                error -> fix_error(FLERR,this,"num_diffusant > 0 is required");
            hasargs = true;
            iarg_ ++;
        }
        else if (strcmp(arg[iarg_],"diffusant_names") == 0)
        {
            if (n_diff != 1)
                error -> fix_error (FLERR,this, "have to define n_diff before 'diffusant_names'");
            if (iarg_ + num_diffusant > narg)
                error -> fix_error(FLERR,this, "Wrong number of arguments");
            diffusant_names_ = new char*[num_diffusant];
            iarg_++;
            for (int i = 0; i < num_diffusant; ++i)
            {
                diffusant_names_[i] = new char [strlen(arg[iarg_])+strlen("_diffCoeff")+1];
                strcpy(diffusant_names_[i], arg[iarg_]);
                strcat(diffusant_names_[i],"_diffCoeff");
                iarg_ += 1;

                /* if(screen)
                {
                    fprintf(screen, "diffusant_names in cfd/coupling/chemistry: %s \n",diffusant_names_[i]);
                }*/
            }
            hasargs = true;
        }
        else if(strcmp(arg[iarg_],"transfer_reactant") == 0)
        {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_reactant'");
            iarg_++;
            if(strcmp(arg[iarg_],"yes") == 0)
                use_reactant_ = true;
            else if(strcmp(arg[iarg_],"no") == 0)
                use_reactant_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_reactant'");
            iarg_++;
            hasargs = true;
        }
        else if (strcmp(this->style,"couple/cfd/chemistry") == 0)
        {
            error->fix_error(FLERR,this,"unknown keyword");
        }
   }

   mod_spec_names_ = new char*[num_species];
   molarfraction_names = new char *[num_species];
   for (int i = 0; i < num_species; ++i)
   {
       mod_spec_names_[i] = new char [strlen("Modified_")+strlen(species_names_[i])+1];
       strcpy(mod_spec_names_[i],"Modified_");
       strcat(mod_spec_names_[i],species_names_[i]);

       molarfraction_names[i] = new char [strlen("X_")+strlen(species_names_[i])+1];
       strcpy(molarfraction_names[i],"X_");
       strcat(molarfraction_names[i],species_names_[i]);
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
    for (int i=0;i<num_species;++i)
    {
        if (species_names_[i]) delete [] species_names_[i];
        if (mod_spec_names_[i]) delete [] mod_spec_names_[i];
        if (molarfraction_names[i]) delete [] molarfraction_names[i];
    }

    for (int j=0; j<num_diffusant;++j)
        if (diffusant_names_[j]) delete [] diffusant_names_[j];

    delete [] species_names_;
    delete [] mod_spec_names_;
    delete [] diffusant_names_;
    delete [] molarfraction_names;

    if(fix_masschange_)     delete []fix_masschange_;
    if(fix_diffusionCoeff_) delete []fix_diffusionCoeff_;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_tgas_)          modify -> delete_fix("partTemp");
    if(unfixflag && fix_rhogas_)        modify -> delete_fix("partRho");
    if(unfixflag && fix_reactionheat_)  modify -> delete_fix("reactionHeat");
    if(unfixflag && fix_totalmole_)     modify -> delete_fix("partMolarConc");
    if(unfixflag && fix_nufField_)      modify -> delete_fix("partNu");
    if(unfixflag && fix_partReynolds_)  modify -> delete_fix("partRe");
    if(unfixflag && fix_partPressure_)  modify -> delete_fix("partP");

    for (int i = 0; i < num_species; ++i)
    {
        if (unfixflag && fix_masschange_[i])        modify -> delete_fix(mod_spec_names_[i]);
        if (unfixflag && fix_molarfraction_[i])     modify -> delete_fix(molarfraction_names[i]);
    }

    for (int j = 0; j < num_diffusant; ++j)
        if (unfixflag && fix_diffusionCoeff_[j])    modify -> delete_fix(diffusant_names_[j]);
    
    if(unfixflag && fix_reactantPerParticle_) modify->delete_fix("reactantPerParticle");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingChemistry::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_create()
{
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
        fixarg[6]="yes";                // communicate ghost no
        fixarg[7]="no";                 // communicate rev
        fixarg[8]="0.";
        fix_tgas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
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
        fixarg[5]="yes";                 // restart yes
        fixarg[6]="no";                  // communicate ghost no
        fixarg[7]="no";                  // communicate rev
        fixarg[8]="0.";
        fix_rhogas_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register reactionheat
    if (!fix_reactionheat_)
    {
        const char* fixarg[9];
        fixarg[0]="reactionHeat";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="reactionHeat";
        fixarg[4]="scalar";     // 1 vector per particle to be registered
        fixarg[5]="yes";        // restart
        fixarg[6]="yes";         // communicate ghost
        fixarg[7]="no";         // communicate rev
        fixarg[8]="0.";
        fix_reactionheat_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register partMolarConc
    if (!fix_totalmole_)
    {
        const char* fixarg[9];
        fixarg[0]="partMolarConc";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partMolarConc";
        fixarg[4]="scalar";     // 1 vector per particle to be registered
        fixarg[5]="no";        // restart
        fixarg[6]="yes";         // communicate ghost
        fixarg[7]="no";         // communicate rev
        fixarg[8]="0.";
        fix_totalmole_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register nu field
    if (!fix_nufField_)
    {
        const char* fixarg[9];
        fixarg[0]="partNu";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partNu";
        fixarg[4]="scalar";     // 1 vector per particle to be registered
        fixarg[5]="yes";        // restart
        fixarg[6]="no";         // communicate ghost
        fixarg[7]="no";         // communicate rev
        fixarg[8]="0.";
        fix_nufField_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register part Reynolds number
    if (!fix_partReynolds_)
    {
        const char* fixarg[9];
        fixarg[0]="partRe";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partRe";
        fixarg[4]="scalar";     // 1 vector per particle to be registered
        fixarg[5]="yes";        // restart
        fixarg[6]="no";         // communicate ghost
        fixarg[7]="no";         // communicate rev
        fixarg[8]="0.";
        fix_partReynolds_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    fix_masschange_ = new FixPropertyAtom*[num_species];
    fix_molarfraction_  =   new FixPropertyAtom*[num_species];
    for (int i=0; i<num_species;++i)
    {
        fix_masschange_[i]      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(mod_spec_names_[i],"property/atom","scalar",0,0,style,false));
        if (fix_masschange_[i] == NULL)
        // register masschange/changeOfSpeciesMass
        {
            const char* fixarg[9];
            fixarg[0]=mod_spec_names_[i];
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=mod_spec_names_[i];
            fixarg[4]="scalar";      // 1 vector per particle to be registered
            fixarg[5]="yes";         // restart
            fixarg[6]="yes";          // communicate ghost
            fixarg[7]="no";          // communicate rev
            fixarg[8]="0.";
            fix_masschange_[i] = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }

        fix_molarfraction_[i]   =   static_cast<FixPropertyAtom*>(modify->find_fix_property(molarfraction_names[i],"property/atom","scalar",0,0,style,false));
        if (fix_molarfraction_[i] == NULL)
        // register fix for molar fraction
        {
            const char* fixarg[9];
            fixarg[0]=molarfraction_names[i];
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=molarfraction_names[i];
            fixarg[4]="scalar";        // 1 vector per particle to be registered
            fixarg[5]="yes";           // restart
            fixarg[6]="yes";            // communicate ghost
            fixarg[7]="no";            // communicate rev
            fixarg[8]="0.";
            fix_molarfraction_[i] = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }
    }

    fix_diffusionCoeff_ = new FixPropertyAtom*[num_diffusant];
    for (int j=0; j<num_diffusant;++j)
    {

        fix_diffusionCoeff_[j]  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffusant_names_[j],"property/atom","scalar",0,0,style,false));
        // register diffusion coefficient variables
        if(fix_diffusionCoeff_[j] == NULL)
        {
            const char* fixarg[9];
            fixarg[0]=diffusant_names_[j];
            fixarg[1]="all";
            fixarg[2]="property/atom";
            fixarg[3]=diffusant_names_[j];
            fixarg[4]="scalar";        // 1 vector per particle to be registered
            fixarg[5]="no";           // restart
            fixarg[6]="yes";            // communicate ghost
            fixarg[7]="yes";            // communicate rev
            fixarg[8]="0.";
            fix_diffusionCoeff_[j] = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        }
    }

    fix_reactantPerParticle_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("reactantPerParticle","property/atom","scalar",0,0,style,false));
    if(!fix_reactantPerParticle_ && use_reactant_)
    {
        const char* fixarg[9];
        fixarg[0]="reactantPerParticle";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="reactantPerParticle";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_reactantPerParticle_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register pressure (partP)
    if (!fix_partPressure_)
    {
        const char* fixarg[9];
        fixarg[0]="partP";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partP";
        fixarg[4]="scalar";              // 1 vector per particle to be registered
        fixarg[5]="yes";                 // restart yes
        fixarg[6]="no";                  // communicate ghost no
        fixarg[7]="no";                  // communicate rev
        fixarg[8]="0.";
        fix_partPressure_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::init()
{
    // make sure there is only one fix of thfor(int i = 1; i <= n_FixMesh_; ++i)is style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/chemistry needs a fix of type couple/cfd");

    // reference to partTemp, partRho and reactionheat
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas_         =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));
    fix_totalmole_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partMolarConc","property/atom","scalar",0,0,style));
    fix_nufField_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partNu","property/atom","scalar",0,0,style));
    fix_partReynolds_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRe","property/atom","scalar",0,0,style));
    fix_partPressure_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partP","property/atom","scalar",0,0,style));

    /*for (int i = 0; i < num_species; ++i)
    {
        //fix_masschange_[i]      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(mod_spec_names_[i],"property/atom","scalar",0,0,style));
        //fix_diffusionCoeff_[i]  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffusant_names_[i],"property/atom","scalar",0,0,style));
        //fix_molarfraction_[i]   =   static_cast<FixPropertyAtom*>(modify->find_fix_property(molarfraction_names[i],"property/atom","scalar",0,0,style));
    }*/

    // values to come from OF
    fix_coupling_->add_pull_property("partTemp","scalar-atom");
    fix_coupling_->add_pull_property("partRho","scalar-atom");
    fix_coupling_->add_pull_property("partMolarConc","scalar-atom");
    fix_coupling_->add_pull_property("partNu","scalar-atom");
    fix_coupling_->add_pull_property("partRe","scalar-atom");
    fix_coupling_->add_pull_property("partP","scalar-atom");

    for (int i=0; i<num_species; ++i)
        fix_coupling_->add_pull_property(molarfraction_names[i],"scalar-atom");

    for (int j=0; j < num_diffusant; ++j)
        fix_coupling_ -> add_pull_property(diffusant_names_[j],"scalar-atom");

    
    if(use_reactant_) fix_coupling_->add_pull_property("reactantPerParticle","scalar-atom");

    //  values to be transfered to OF
    fix_coupling_->add_push_property("reactionHeat","scalar-atom");

    for (int i = 0; i<num_species; ++i)
    {
        fix_coupling_->add_push_property(mod_spec_names_[i],"scalar-atom");
    }

     /*for (int j=0; j < num_diffusant; ++j)
        {
            double *diff_Coeff = fix_diffusionCoeff_[j] -> vector_atom;
            for (int i =0; i < atom->nlocal; i++)
                if (screen)
                    fprintf(screen,"diffCoeff from cfd/coupling/chemistry: %6.15f \n", diff_Coeff[i]);
        }*/
    }

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::initial_integrate(int)
{
    /* if (comm -> me == 0 && screen)
        fprintf(screen,"activate initial integrate \n"); */
    // for all species, reaction heat
    // if current timestep - 1 == latestpush(species name)
    // reset fix_masschange_(species name)
    // -1 is needed because time step is advanced before this function is called
    bigint prev_time = update->ntimestep - 1;
    int *mask   = atom -> mask;
    int  nlocal = atom -> nlocal;

    for (int k = 0; k < num_species; ++k)
    {
        if (prev_time == fix_coupling_ -> latestpush(mod_spec_names_[k]))
        {
            for (int i = 0; i < nlocal; ++i)
            {
                if (mask[i] & groupbit)
                     fix_masschange_[k] -> vector_atom[i] = 0.;
            }
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

