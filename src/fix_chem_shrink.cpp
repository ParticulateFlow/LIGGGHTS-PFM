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
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "fix_chem_shrink.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrink::FixChemShrink(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg),
    screenflag_ (0),
    speciesA(NULL),
    speciesC(NULL),
    massA(NULL),
    massC(NULL),
    moleFracA(NULL),
    moleFracC(NULL),
    fix_changeOfA_(NULL),
    fix_changeOfC_(NULL),
    fix_rhogas(NULL),
    fix_tgas(NULL),
    fix_reactionheat_(NULL),
    fix_totalMole_(NULL),
    fix_nuField_(NULL),
    fix_partRe_(NULL),
    fix_reactantPerParticle_(NULL),
    fix_moleFractionA_(NULL),
    fix_moleFractionC_(NULL),
    xA_(NULL),
    xC_(NULL),
    nuf_(NULL),
    Rep_(NULL),
    iarg_(3),
    k0(0),
    molMass_A_(0),
    molMass_B_(0),
    molMass_C_(0),
    relaxFac_(0.5),
    nu_A_(1),
    nu_B_(1),
    nu_C_(1),
    radius_(NULL),
    pmass_(NULL),
    pdensity_(NULL),
    rmin(0.),
    minMolarFrac_(1e-3),
    changeOfA_(NULL),
    changeOfC_(NULL),
    rhogas_(NULL),
    tgas_(NULL),
    molarConc_(NULL),
    reactantPerParticle_(NULL),
    spcA(0),
    spcC(0),
    TimeStep(0.0),
    current_timestep(0),
    use_reactant_(false)
{
    if ((strncmp(style,"chem/shrink",11) == 0) && (!atom->radius_flag||!atom->rmass_flag))
            error -> all (FLERR,"Fix chem/shrink needs per particle radius and mass");

    bool hasargs = true;

    if (narg < 24) error -> all (FLERR,"Fix chem/shrink: Wrong number of arguments");

    while (iarg_ < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg_],"speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            spcA = 1;
            speciesA = new char [strlen(arg[iarg_+1])];
            strcpy(speciesA,arg[iarg_+1]);
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_A_ = atof(arg[iarg_+1]);
            if (molMass_A_ < 0)
                error -> fix_error(FLERR, this, "molar mass of A is not defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nuA") == 0)
        {
            nu_A_ = atoi(arg[iarg_+1]);
            if (nu_A_ < 1)
                error -> fix_error(FLERR, this, "nuA is not defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            spcC = 1;
            speciesC = new char [strlen(arg[iarg_+1])];
            strcpy(speciesC,arg[iarg_+1]);
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"molMassC") == 0)
        {
            if (spcC != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesC' before 'molMassC'");
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_C_ = atof(arg[iarg_+1]);
            if (molMass_C_ < 0)
                error -> fix_error(FLERR, this, "molMassC > 0 is required");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nuC") == 0)
        {
            nu_C_ = atoi(arg[iarg_+1]);
            if (nu_C_ < 1)
                error -> fix_error(FLERR, this, "nuC is not defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"molMassB") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_B_ = atof(arg[iarg_+1]);
            if (molMass_B_ < 0)
                error -> fix_error(FLERR, this, "molMassB > 0 is required");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nuB") == 0)
        {
            nu_B_ = atoi(arg[iarg_+1]);
            if (nu_B_ < 1)
                error -> fix_error(FLERR, this, "nuB is not defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"k") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            k0 = atof(arg[iarg_+1]);
            if (k0 <= 0)
                error -> fix_error(FLERR, this, "k is not (well-)defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"rmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            rmin = atof(arg[iarg_+1]);
            if (rmin == 0)
                error -> fix_error(FLERR, this, "rmin is not defined");
            hasargs = true;
            iarg_+=2;
        }
        else if (strcmp(arg[iarg_],"nevery") == 0)
        {
            nevery = atoi(arg[iarg_+1]);
            if (nevery <= 0) error->fix_error(FLERR,this,"");
            hasargs = true;
            iarg_+=2;
        }
        else if(strcmp(arg[iarg_],"use_reactant") == 0)
        {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'use_reactant'");
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
        else if (strcmp(arg[iarg_],"relaxFactor") == 0)
        {
            relaxFac_ = atof(arg[iarg_+1]);
            if (relaxFac_ < 0)
                error -> fix_error(FLERR, this, "relaxFactor is not well-defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"screen") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) screenflag_ = 1;
            else if (strcmp(arg[iarg_+1],"no") == 0) screenflag_ = 0;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(this->style,"chem/shrink") == 0)
        {
            error->fix_error(FLERR,this,"necessary keyword is missing");
        }
    }

    // define changed species mass A
    massA = new char [strlen("Modified_")+strlen(speciesA)+1];
    strcpy(massA,"Modified_");
    strcat(massA,speciesA);

    // define changed species mass C
    massC = new char [strlen("Modified_")+strlen(speciesC)+1];
    strcpy(massC,"Modified_");
    strcat(massC,speciesC);

    // reactant species molar fraction
    moleFracA = new char [strlen("X_")+strlen(speciesA)+1];
    strcpy(moleFracA,"X_");
    strcat(moleFracA,speciesA);

    // product species molar fraction
    moleFracC = new char [strlen("X_")+strlen(speciesC)+1];
    strcpy(moleFracC,"X_");
    strcat(moleFracC,speciesC);

    time_depend = 1;
    force_reneighbor =1;
    next_reneighbor = update ->ntimestep + nevery;
    restart_global = 0;
}


/* ---------------------------------------------------------------------- */

FixChemShrink::~FixChemShrink()
{
    delete []massA;
    delete []massC;

    delete []speciesA;
    delete []speciesC;

    delete []moleFracA;
    delete []moleFracC;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if(fix_moleFractionA_)       modify  ->  delete_fix(moleFracA);
        if(fix_moleFractionC_)       modify  ->  delete_fix(moleFracC);
        if(fix_rhogas)               modify  ->  delete_fix("partRho");
        if(fix_tgas)                 modify  ->  delete_fix("partTemp");
        if(fix_reactionheat_)        modify  ->  delete_fix("reactionHeat");
        if(fix_changeOfA_)           modify  ->  delete_fix(massA);
        if(fix_changeOfC_)           modify  ->  delete_fix(massC);
        if(fix_totalMole_)           modify  ->  delete_fix("partMolarConc");
        if(fix_reactantPerParticle_) modify  ->  delete_fix("reactantPerParticle");
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

void FixChemShrink::updatePtrs()
{
    xA_         =   fix_moleFractionA_  ->  vector_atom;
    xC_         =   fix_moleFractionC_  ->  vector_atom;
    changeOfA_  =   fix_changeOfA_  ->  vector_atom;
    changeOfC_  =   fix_changeOfC_  ->  vector_atom;
    rhogas_     =   fix_rhogas      ->  vector_atom;
    molarConc_  =   fix_totalMole_  ->  vector_atom;

    if(use_reactant_)
    {
        reactantPerParticle_ = fix_reactantPerParticle_ -> vector_atom;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::reaction()
{
    updatePtrs();
    int nlocal  =   atom -> nlocal;
    double molarConc = 0.0;
    double molarFrac = 0.0;

    for (int i = 0 ; i < nlocal; i++)
    {
        if(radius_[i] > rmin)
        {
            if (screenflag_ && screen)
            {
                fprintf(screen,"check TimeStep %f \n",TimeStep);
                fprintf(screen,"check current timestep %d \n", current_timestep);
            }

            // Current time step concentration change of reactant A and product C
            double dA   =   0.0;
            molarConc = molarConc_[i];
            molarFrac = xA_[i];
            if(molarFrac < minMolarFrac_)
            {
                continue;
            }
            double k = reactionRatConst(i);
            if(nu_A_ < 2)
            {
                dA = -k*molarFrac*molarConc*molMass_A_*partSurfArea(radius_[i])*TimeStep*nevery;
            }
            else
            {
                dA = -k*pow((molarFrac*molarConc),nu_A_)*partSurfArea(radius_[i])*molMass_A_*nu_A_*TimeStep*nevery;
            }

            // limit mass change - can't remove more than present in cell
            // limit it with species mass per volume x voidfraction x cell volume / particles in cell x relaxation factor (0.8)
            if(use_reactant_)
            {
                if (screenflag_ && screen) fprintf(screen,"checking reactant limitation\n");

                double dAmax = molarFrac * molarConc * molMass_A_ * reactantPerParticle_[i] * relaxFac_;

                if(-dA > dAmax) dA = -dAmax;
            }

            double dC   =   -dA*molMass_C_*nu_C_/(molMass_A_*nu_A_);
            // mass removed from particle
            double dB   =    dA*molMass_B_*nu_B_/(molMass_A_*nu_A_);

            // total rate of change for species A
            changeOfA_[i]       +=  dA;

            // total rate of change for species C
            changeOfC_[i]       +=  dC;

            // Mass of single particle
            // never remove more than half the particle's mass at once
            if(-dB > 0.5*pmass_[i])
            {
                pmass_[i] *= 0.5;
            }
            else
            {
                pmass_[i] += dB;
            }


            // change of radius of particle -assumption: density of particle is constant
            radius_[i]           =   cbrt(0.75*pmass_[i]/(M_PI*pdensity_[i]));

            // uncomment if postproc (verification)
            if (screenflag_ && screen)
            {
                fprintf(screen,"Co2 Mass %f \n",changeOfC_[i]);
                fprintf(screen,"O2 Mass %f \n",changeOfA_[i]);
                fprintf(screen,"Particle Mass %f \n",pmass_[i]);
                fprintf(screen,"Gas Density %f \n",rhogas_[i]);
                fprintf(screen,"total number of moles in chem/shrink: %f \n",molarConc_[i]);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

double FixChemShrink::reactionRatConst(int i)
{
    return k0;
}

/* ----------------- compute particle surface area ------------------------ */

double FixChemShrink::partSurfArea(double radius)
{
    return 4.*M_PI*radius*radius;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::init()
{
    // error checks
    if (!atom->radius_flag)
      error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");
    if (!atom->rmass_flag)
      error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");
    if (!atom->tag_enable || 0 == atom->map_style)
      error->fix_error(FLERR,this,"requires atom tags and an atom map");

    // references
    fix_moleFractionA_  =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(moleFracA,"property/atom","scalar",0,0,style));
    fix_moleFractionC_  =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(moleFracC,"property/atom","scalar",0,0,style));
    fix_changeOfA_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massA,"property/atom","scalar",0,0,style));
    fix_changeOfC_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massC,"property/atom","scalar",0,0,style));
    fix_tgas            =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas          =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));
    fix_totalMole_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partMolarConc","property/atom","scalar",0,0,style));

    if(use_reactant_)
    {
        fix_reactantPerParticle_ = static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactantPerParticle","property/atom","scalar",0,0,style));
    }

    updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;
    TimeStep = update -> dt;
    current_timestep = update->ntimestep;

    // skip if integration not wanted at this timestep
    if (current_timestep % nevery)
    {
        return;
    }

    reaction();
}
