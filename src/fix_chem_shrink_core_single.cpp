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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "fix_chem_shrink_core_single.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_property_atom_polydispparcel.h"
#include "force.h"
#include "group.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define SMALL   1e-10

const double FixChemShrinkCoreSingle::Runiv = 8.3144;

const double FixChemShrinkCoreSingle::a_coeff_CO[]    = { 200.00, 6000., 1000.,
                                                    3.048486E+00,  1.351728E-03, -4.857941E-07,  7.885364E-11, -4.698075E-15, -1.426612E+04,  6.017098E+00,
                                                    3.579534E+00, -6.103537E-04,  1.016814E-06,  9.070059E-10, -9.044245E-13, -1.434409E+04,  3.508409E+00 };
const double FixChemShrinkCoreSingle::a_coeff_CO2[]   = { 200.00, 6000., 1000.,
                                                    4.636511E+00,  2.741457E-03, -9.958976E-07,  1.603867E-10, -9.161986E-15, -4.902490E+04, -1.934896E+00,
                                                    2.356813E+00,  8.984130E-03, -7.122063E-06,  2.457301E-09, -1.428855E-13, -4.837197E+04,  9.900904E+00 };
const double FixChemShrinkCoreSingle::a_coeff_O2[]    = { 200.00, 6000., 1000.,
                                                    3.660960E+00,  6.563655E-04, -1.411494E-07,  2.057976E-11, -1.299132E-15, -1.215977E+03,  3.415362E+00,
                                                    3.782456E+00, -2.996734E-03,  9.847302E-06, -9.681295E-09,  3.243728E-12, -1.063943E+03,  3.657676E+00 };
// coke data are a linear fit to Patisson and Hanrot, Metallurgical and Materials Transactions B 31.2 (2000): 381-390.
// sixth coefficient is needed to get heat of formation at room temperature
const double FixChemShrinkCoreSingle::a_coeff_coke[]  = { 298.00, 1300., 1000.,
                                                    4.040000E-01,  2.020000E-03,  0.000000E+00,  0.000000E+00,  0.000000E+00, -2.100800E+02,  0.000000E+00,
                                                    4.040000E-01,  2.020000E-03,  0.000000E+00,  0.000000E+00,  0.000000E+00, -2.100800E+02,  0.000000E+00 };
// no reasonable data for ash available; due to low mass content, it should not have any significant impact
const double FixChemShrinkCoreSingle::a_coeff_ash[]  = { 298.00, 1300., 1000.,
                                                    4.040000E-01,  2.020000E-03,  0.000000E+00,  0.000000E+00,  0.000000E+00, -2.100800E+02,  0.000000E+00,
                                                    4.040000E-01,  2.020000E-03,  0.000000E+00,  0.000000E+00,  0.000000E+00, -2.100800E+02,  0.000000E+00 };


/* ---------------------------------------------------------------------- */

FixChemShrinkCoreSingle::FixChemShrinkCoreSingle(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    layers_(1),
    minMolarFrac_(1e-3),
    rmin_(1e-5),      //  [m]
    layerDiffusion_(true),
    heatToFluid_(false),
    heatToParticle_(false),
    shrink_(true),
    reactionHeatIndex_(-1),
    KeqIndex_(0),
    k0_(-1.0),
    T0_(-1.0),
    Tmin_(0.0),
    nPreFactor_(0),
    Cp_coke_(10.2),
    T_room_(298)
{
    if ((strncmp(style, "chem/shrink/core/single", 15) == 0) && ((!atom->radius_flag) || (!atom->rmass_flag)))
        error->all(FLERR, "Fix chem/shrink/core/single needs per particle radius and mass");

    // set defaults
    init_defaults();
    screenflag_ = 0;
    cg_ = -1.0;
    int iarg_ = 3;
    bool hasargs = true;

    while (iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_], "speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error->fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            speciesA = new char[strlen(arg[iarg_+1])+1];
            strcpy(speciesA, arg[iarg_+1]);
            if (screenflag_ && screen) fprintf(screen, "speciesA: %s", speciesA);
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_], "molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_A_ = atof(arg[iarg_ + 1]);
            if (molMass_A_ < 0.0)
                error->fix_error(FLERR, this, "molar mass of A is not defined");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nuA") == 0)
        {
            nu_A_ = atoi(arg[iarg_+1]);
            if (nu_A_ < 1)
                error -> fix_error(FLERR, this, "nuA is not well-defined");
            if (nu_A_ > 1)
                error -> fix_error(FLERR, this, "nuA>1 not implemented");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"molMassB") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_B_ = atof(arg[iarg_+1]);
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nuB") == 0)
        {
            nu_B_ = atoi(arg[iarg_+1]);
            if (nu_B_ < 1)
                error -> fix_error(FLERR, this, "nuB is not well-defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_], "speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            speciesC = new char[strlen(arg[iarg_+1])+1];
            strcpy(speciesC, arg[iarg_+1]);
            if (screenflag_ && screen) fprintf(screen, "speciesC: %s", speciesC);
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_], "molMassC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_C_ = atof(arg[iarg_ + 1]);
            if (molMass_C_ < 0.0)
                error->fix_error(FLERR, this, "molar mass of C is not defined");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nuC") == 0)
        {
            nu_C_ = atoi(arg[iarg_+1]);
            if (nu_C_ < 1)
                error -> fix_error(FLERR, this, "nuC is not well-defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"molMassD") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_D_ = atof(arg[iarg_+1]);
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nuD") == 0)
        {
            nu_D_ = atoi(arg[iarg_+1]);
            if (nu_D_ < 1)
                error -> fix_error(FLERR, this, "nuB is not well-defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"layerDiffusion") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) layerDiffusion_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) layerDiffusion_ = false;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"heatToParticle") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) heatToParticle_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) heatToParticle_ = false;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"heatToFluid") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) heatToFluid_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) heatToFluid_ = false;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"shrink") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) shrink_ = true;
            else if (strcmp(arg[iarg_+1],"no") == 0) shrink_ = false;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"k") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            k0_ = atof(arg[iarg_+1]);
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"T") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            T0_ = atof(arg[iarg_+1]);
            if (T0_ <= 0.)
                error -> fix_error(FLERR, this, "T is not (well-)defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"Tmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            Tmin_ = atof(arg[iarg_+1]);
            if (Tmin_ <= 0.)
                error -> fix_error(FLERR, this, "Tmin is not (well-)defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nPreFactor") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            nPreFactor_ = atoi(arg[iarg_+1]);
            if (nPreFactor_ < 0)
                error -> fix_error(FLERR, this, "nPreFactor is not (well-)defined");
            hasargs = true;
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"rmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            rmin_ = atof(arg[iarg_+1]);
            if (rmin_ < SMALL)
                error -> fix_error(FLERR, this, "rmin is not well-defined");
            hasargs = true;
            iarg_+=2;
        }
        else if (strcmp(arg[iarg_],"nevery") == 0)
        {
            nevery = atoi(arg[iarg_+1]);
            if (nevery <= 0) error->fix_error(FLERR,this,"nevery must be larger than 0");
            iarg_+=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"cg") == 0)
        {
            cg_ = atof(arg[iarg_+1]);
            if (cg_ <= 0.0) error->fix_error(FLERR,this,"cg must be larger than 0");
            iarg_+=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"scale_reduction_rate") == 0 )
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR,this,"Wrong number of arguments");
            scale_reduction_rate = atof(arg[iarg_ + 1]);
            if (scale_reduction_rate < 0.0)
                error->fix_error(FLERR, this, "scale_reduction_rate cannot be smaller or less then 0");
            iarg_ += 2;
            hasargs = true;
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
        else if(strcmp(arg[iarg_],"limit_reactant_consumption") == 0)
        {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'limit_reactant_consumption'");
            iarg_++;
            if(strcmp(arg[iarg_],"yes") == 0)
                limit_reactant_consumption_ = true;
            else if(strcmp(arg[iarg_],"no") == 0)
                limit_reactant_consumption_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'limit_reactant_consumption'");
            iarg_++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"maxReactantConsumptionFrac") == 0)
        {
            maxReactantConsumptionFrac_ = atof(arg[iarg_+1]);
            if (maxReactantConsumptionFrac_ < 0)
                error -> fix_error(FLERR, this, "maxReactantConsumptionFractor is not well-defined");
            hasargs = true;
            iarg_ +=2;
        }
        else
        {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    if (k0_ < 0.0 || T0_ < 0.0)
    {
        error->fix_error(FLERR,this,"specify positive values for k0 and T0");
    }

    if (molMass_D_ < SMALL * molMass_B_)
    {
        error->fix_error(FLERR,this,"set sensible value for molar mass of solid reaction product");
    }

    if (heatToParticle_ && heatToFluid_)
    {
        error->fix_error(FLERR,this,"trying to apply reaction heat to both particles and fluid");
    }

    // check if one of the predefined reaction enthalpies matches current reaction
    if (heatToParticle_ || heatToFluid_)
    {
        if (strcmp(speciesA,"O2") == 0 && strcmp(speciesC,"CO") == 0) reactionHeatIndex_ = 0;
        else if (strcmp(speciesA,"CO2") == 0 && strcmp(speciesC,"CO") == 0) reactionHeatIndex_ = 1;
        else error->fix_error(FLERR,this,"trying to use reaction heat for unknown reaction");
    }
    // check if current reaction has a non-trivial, predefined equilibrium constant
    if (strcmp(speciesA,"CO2") == 0 && strcmp(speciesC,"CO") == 0) KeqIndex_ = 1;


    // define changed species mass A
    massA = new char [strlen("Modified_")+strlen(speciesA)+1];
    strcpy(massA,"Modified_");
    strcat(massA,speciesA);

    // define changed species mass C
    massC = new char [strlen("Modified_")+strlen(speciesC)+1];
    strcpy(massC,"Modified_");
    strcat(massC,speciesC);

    // define diffusant species
    diffA = new char [strlen(speciesA)+strlen("_diffCoeff")+1];
    strcpy(diffA,speciesA);
    strcat(diffA,"_diffCoeff");

    // reacting species bulk molar fraction
    moleFracA = new char [strlen("X_")+strlen(speciesA)+1];
    strcpy(moleFracA,"X_");
    strcat(moleFracA,speciesA);

    // product species bulk molar fraction
    moleFracC = new char [strlen("X_")+strlen(speciesC)+1];
    strcpy(moleFracC,"X_");
    strcat(moleFracC,speciesC);

    time_depend = 1;
    if (cg_ < 0.0) cg_ = force->cg();
    peratom_flag = 1;
    peratom_freq = 1;
    global_freq = 1;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::post_create()
{
    if (fix_Aterm == NULL) {
        char *fixname = new char [strlen("Aterm_")+strlen(id)+1];
        strcpy(fixname,"Aterm_");
        strcat(fixname,id);

        const char* fixarg[9];
        // register property/atom for chemical reaction resistance
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";          // vector with 3 values (for the different layers)
        fixarg[5]="yes";             // restart yes
        fixarg[6]="no";              // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";             // tale 0 as default value
        fix_Aterm = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete [] fixname;
    }

    if (fix_Bterm == NULL) {
        char *fixname = new char[strlen("Bterm_")+strlen(id)+1];
        strcpy(fixname,"Bterm_");
        strcat(fixname,id);

        const char* fixarg[9];
        // register property/atom for diffusion resistance term
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";          // vector with 3 values (for the different layers)
        fixarg[5]="yes";             // restart yes
        fixarg[6]="no";              // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";
        fix_Bterm = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_Massterm == NULL) {
        char *fixname = new char[strlen("Massterm_")+strlen(id)+1];
        strcpy(fixname,"Massterm_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";             // restart yes
        fixarg[6]="yes";             // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";
        fix_Massterm = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_effDiffBinary == NULL) {
        char *fixname = new char[strlen("effDiffBinary_")+strlen(id)+1];
        strcpy(fixname,"effDiffBinary_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup]; //"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fix_effDiffBinary = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete [] fixname;
    }

    if (fix_effDiffKnud == NULL) {
        char *fixname = new char[strlen("effDiffKnud_")+strlen(id)+1];
        strcpy(fixname,"effDiffKnud_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fix_effDiffKnud = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_dY_ == NULL) {
        char* fixname = new char[strlen("dY_")+strlen(id)+1];
        strcpy(fixname,"dY_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fix_dY_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_dmA_ == NULL) {
        char* fixname = new char[strlen("dmA_")+strlen(id)+1];
        strcpy(fixname,"dmA_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fix_dmA_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        delete []fixname;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_Aterm)          { modify->delete_fix(fix_Aterm->id); Aterm = NULL; }
        if (fix_Bterm)          { modify->delete_fix(fix_Bterm->id); Bterm = NULL; }
        if (fix_Massterm)       { modify->delete_fix(fix_Massterm->id); Massterm = NULL; }
        if (fix_effDiffBinary)  { modify->delete_fix(fix_effDiffBinary->id); effDiffBinary = NULL; }
        if (fix_effDiffKnud)    { modify->delete_fix(fix_effDiffKnud->id); effDiffKnud = NULL; }
        if (fix_dY_)            { modify->delete_fix(fix_dY_->id); dY = NULL; }
        if (fix_dmA_)           { modify->delete_fix(fix_dmA_->id); dmA_f_ = NULL; }
    }
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCoreSingle::~FixChemShrinkCoreSingle()
{
    delete [] massA;
    delete [] massC;
    delete [] diffA;
    delete [] moleFracA;
    delete [] moleFracC;

    delete [] speciesA;
    delete [] speciesC;
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCoreSingle::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::updatePtrs()
{
    changeOfA_      =   fix_changeOfA_      ->  vector_atom;
    changeOfC_      =   fix_changeOfC_      ->  vector_atom;
    T_              =   fix_tgas_           ->  vector_atom;
    Tpart_          =   fix_tpart_          ->  vector_atom;
    heatFlux_       =   fix_heatFlux_       ->  vector_atom;
    molecularDiffusion_  = fix_diffcoeff_   ->  vector_atom;
    nuf_            =   fix_nuField_        ->  vector_atom;
    Rep_            =   fix_partRe_         ->  vector_atom;
    xA_             =   fix_moleFractionA_  ->  vector_atom;
    xC_             =   fix_moleFractionC_  ->  vector_atom;
    relRadii_       =   fix_layerRelRad_    ->  array_atom;
    massLayer_      =   fix_layerMass_      ->  array_atom;
    effvolfactors_  =   fix_polydisp_       ->  vector_atom;

#ifdef PER_ATOM_LAYER_DENSITIES
    layerDensities_ =   fix_layerDens_      ->  array_atom;
#else
    layerDensities_ =   fix_layerDens_      ->  values;
#endif

    porosity_       =   fix_porosity_       ->  values;
    rhoeff_         =   fix_rhoeff_         ->  array_atom;
    tortuosity_     =   fix_tortuosity_     ->  compute_scalar();
    pore_diameter_  =   fix_pore_diameter_  ->  compute_scalar();
    fracRed_        =   fix_fracRed         ->  array_atom;
    Aterm           =   fix_Aterm           ->  vector_atom;
    Bterm           =   fix_Bterm           ->  vector_atom;
    Massterm        =   fix_Massterm        ->  vector_atom;
    effDiffBinary   =   fix_effDiffBinary   ->  vector_atom;
    effDiffKnud     =   fix_effDiffKnud     ->  vector_atom;
    partP_          =   fix_partPressure_   ->  vector_atom;

    TimeStep        =   update  -> dt;
    radius_         =   atom    -> radius;
    pmass_          =   atom    -> rmass;
    pdensity_       =   atom    -> density;

    dY          =   fix_dY_         ->  vector_atom;
    dmA_f_      =   fix_dmA_ -> vector_atom;
    reactionheat_ = fix_reactionheat -> vector_atom;

    if(limit_reactant_consumption_)
    {
        reactantPerParticle_ = fix_reactantPerParticle_ -> vector_atom;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::init()
{
    // error checks
    if (!atom->radius_flag)
      error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");
    if (!atom->rmass_flag)
      error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");
    if (!atom->tag_enable || 0 == atom->map_style)
      error->fix_error(FLERR,this,"requires atom tags and an atom map");

    char* propertyname = new char [strlen("density_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"density_");
    strcat(propertyname,group->names[igroup]);
#ifdef PER_ATOM_LAYER_DENSITIES
    fix_layerDens_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname,"property/atom","vector",0,0,style));
#else
    fix_layerDens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname,"property/global","vector",2,0,style));
#endif
    delete [] propertyname;

    // references for per atom properties.
    fix_changeOfA_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massA, "property/atom", "scalar", 0, 0, style));
    fix_changeOfC_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massC, "property/atom", "scalar", 0, 0, style));
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, style)); // gas temperature at location of particle
    fix_tpart_          =   static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp", "property/atom", "scalar", 0, 0, style));     // particle temperature
    fix_heatFlux_       =   static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux", "property/atom", "scalar", 0, 0, style));
    fix_diffcoeff_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffA, "property/atom", "scalar", 0, 0, style));
    fix_nuField_        =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partNu", "property/atom", "scalar", 0, 0, style));
    fix_partRe_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partRe", "property/atom", "scalar", 0, 0, style));
    fix_moleFractionA_  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(moleFracA, "property/atom", "scalar", 0, 0, style));
    fix_moleFractionC_  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(moleFracC, "property/atom", "scalar", 0, 0, style));
    fix_partPressure_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partP", "property/atom", "scalar", 0, 0, style));
    fix_layerRelRad_    =   static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii", "property/atom", "vector", 0, 0, style));
    fix_layerMass_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer","property/atom","vector",0,0,style));
    fix_thermal_capacity_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("thermalCapacity","property/atom","scalar",0,0,style));
    fix_rhoeff_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("rhoeff", "property/atom", "vector", 0, 0, style));
    fix_polydisp_       =   static_cast<FixPropertyAtomPolydispParcel*>(modify->find_fix_property("effvolfactor", "property/atom","scalar",0,0,style));

    if (screenflag_ && screen)
    {
        if (fix_polydisp_ == NULL) fprintf(screen,"fix chem/shrink/core/single: found no fix for particle effective volumina\n");
        else fprintf(screen,"fix chem/shrink/core/single: found fix for particle effective volumina\n");
    }

    propertyname = new char [strlen("porosity_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"porosity_");
    strcat(propertyname,group->names[igroup]);
    fix_porosity_       =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "vector", 0, 0, style));
    delete [] propertyname;

    // references for global properties - valid for every particle of same type equally
    propertyname = new char [strlen("tortuosity_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"tortuosity_");
    strcat(propertyname,group->names[igroup]);
    fix_tortuosity_     =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char [strlen("pore_diameter_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"pore_diameter_");
    strcat(propertyname,group->names[igroup]);
    fix_pore_diameter_  =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "scalar", 0, 0,style));
    delete [] propertyname;

    propertyname = new char [strlen("Aterm_")+strlen(id)+1];
    strcpy (propertyname,"Aterm_");
    strcat(propertyname,id);
    fix_Aterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("Bterm_")+strlen(id)+1];
    strcpy (propertyname,"Bterm_");
    strcat(propertyname,id);
    fix_Bterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("Massterm_")+strlen(id)+1];
    strcpy (propertyname,"Massterm_");
    strcat(propertyname,id);
    fix_Massterm        =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("effDiffBinary_")+strlen(id)+1];
    strcpy (propertyname,"effDiffBinary_");
    strcat(propertyname,id);
    fix_effDiffBinary = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("effDiffKnud_")+strlen(id)+1];
    strcpy (propertyname,"effDiffKnud_");
    strcat(propertyname,id);
    fix_effDiffKnud = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("dY_")+strlen(id)+1];
    strcpy (propertyname,"dY_");
    strcat(propertyname,id);
    fix_dY_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("dmA_")+strlen(id)+1];
    strcpy (propertyname,"dmA_");
    strcat(propertyname,id);
    fix_dmA_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    fix_fracRed         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("fracRed", "property/atom", "vector", 0, 0, style));
    fix_reactionheat   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("reactionHeat", "property/atom", "scalar", 0, 0, style));

    if(limit_reactant_consumption_)
    {
        fix_reactantPerParticle_ = static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactantPerParticle","property/atom","scalar",0,0,style));
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::setup(int)
{
    updatePtrs();
    int *mask   =   atom->mask;
    // set initial values for rhoeff, and use them to calculate mass of layers
    for (int i = 0; i < atom->nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            active_layers(i);

            for (int layer=0; layer <= layers_; layer++)
            {
#ifdef PER_ATOM_LAYER_DENSITIES
                rhoeff_[i][layer] = (1.0 - porosity_[layer])*layerDensities_[i][layer];
#else
                rhoeff_[i][layer] = (1.0 - porosity_[layer])*layerDensities_[layer];
#endif
            }
            calcMassLayer(i);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::post_force(int)
{
    if (update->ntimestep % nevery)
        return;

    updatePtrs();
    int i = 0;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    double x0_eq_= 0.; // molar fraction of reactant gas
    double dmA_ = 0.;   // mass flow rate of reactant gas species for each layer at w->fe, m->w & h->m interfaces

    for (i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (xA_[i] < minMolarFrac_)
            {
                continue;
            }

            if (T_[i] < Tmin_)
            {
                continue;
            }

            // 1st recalculate masses of layers if layer has reduced
            // is ignored if there is no change in layers
            if (active_layers(i) > 0)
            {
                // calculate values for fractional reduction f_i = (1-relRadii_i^3)
                // or with mass ratio provides simplicity for calculations of A & B terms.
                FractionalReduction(i);
                // get values for equilibrium molar fraction of reactant gas species,
                // this value is calculated from the Equilibrium constants function Keq(layer,T).
                // and used in the reaction rate determination.
                getXi(i,x0_eq_);
                // calculate the reaction resistance term
                getA(i);
                // calculate the diffusion resistance term
                getB(i);
                // calculate mass transfer resistance term
                getMassT(i);
                // do the reaction calculation with pre-calculated values of A, B and ß (massT)
                // the USCM model chemical reaction rate with gaseous species model
                // based on the works of Philbrook, Spitzer and Manning
                reaction(i, dmA_, x0_eq_);
                // the results of reaction gives us the mass change of reactant species gas
                // in the usual case that means the CO gas mass species change is given
                // this information is used then to calculate mass changes of particle layers
                update_atom_properties(i, dmA_);
                // also the results of reaction function is used to calculate
                // the changes in gas species
                update_gas_properties(i, dmA_);
                // calculate delta_h, and dot_delta_h for heat of reaction
                heat_of_reaction(i, dmA_);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCoreSingle::active_layers(int i)
{
    layers_ = 1;
    if (relRadii_[i][1]*(radius_[i]/cg_) < rmin_)
    {
        layers_=0;
    }

    if (screenflag_ && screen)
        fprintf(screen, "active layers: %i \n", layers_);
    return layers_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::calcMassLayer(int i)
{
    double rad[2] = {0.};
    rad[1] = (radius_[i]/cg_)*relRadii_[i][1];
    rad[0] = (radius_[i]/cg_)*relRadii_[i][0];

    massLayer_[i][1] = MY_4PI3*rhoeff_[i][1]*rad[1]*rad[1]*rad[1];
    massLayer_[i][0] = MY_4PI3*rhoeff_[i][0]*(rad[0]*rad[0]*rad[0]
                                                         -rad[1]*rad[1]*rad[1]);

    if (fix_polydisp_)
    {
        massLayer_[i][1] *= effvolfactors_[i];
        massLayer_[i][0] *= effvolfactors_[i];
    }
}

/* ---------------------------------------------------------------------- */

double FixChemShrinkCoreSingle::K_eq(int i)
{
    double Keq_ = 1e9; // simple forward reaction
    double T = T_[i];
    if (KeqIndex_ == 1) // CO2 + C -> 2CO
    {
        double exp = 9141/T + 0.000224*T - 9.595;
        Keq_ = 1.0/pow(10,exp);
    }
    return Keq_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::getXi(int i, double &x0_eq_)
{
    const double kch2_ = xA_[i] + xC_[i];

    x0_eq_  =   kch2_/(1.0+K_eq(i));

    if (screenflag_ && screen)
        fprintf(screen, "x0_eq: %f \n", x0_eq_);
}

/* ---------------------------------------------------------------------- */

// calculate A_[i] [s/m] - the chemical reaction resistance term of generalized Arrhenius form k0*T^n*exp(-T0/T)*...
void FixChemShrinkCoreSingle::getA(int i)
{
    double T = T_[i];
    Aterm[i] = k0_ * exp(-T0_ / T)
                        * cbrt((1.0 - fracRed_[i][0]) * (1.0 - fracRed_[i][0]))
                        * (1.0 + 1.0 / K_eq(i));
    for (int j=0;j < nPreFactor_;j++)
    {
        Aterm[i] *= T;
    }
    Aterm[i] = 1.0 / Aterm[i];
}

/* ---------------------------------------------------------------------- */

// Calculate B [s/m] - the diffusion resistance term
// Use binary diffusion for mixture, and Knudsen diffusion to determine the effective diffusion term
void FixChemShrinkCoreSingle::getB(int i)
{
    if (!layerDiffusion_)
    {
        Bterm[i]= 0.0;
        return;
    }

    double fracRedThird = 0.;
    double diffEff = 0.;


    // calculate fractional reduction to the power of 1/3 for simpler use
    fracRedThird = cbrt(1.0-fracRed_[i][0]);

    // Calculate the effective molecular diffusion
    effDiffBinary[i] = molecularDiffusion_[i]*(porosity_[0]/tortuosity_) + SMALL;

    // Calculate the Knudsen diffusion
    effDiffKnud[i]  =  (pore_diameter_/6.0)*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_))*(porosity_[0]/tortuosity_) + SMALL;

    // total effective diffusivity
    diffEff =   effDiffKnud[i]*effDiffBinary[i]/(effDiffBinary[i]+effDiffKnud[i]) + SMALL;


    // calculation of diffusion term
    if (molecularDiffusion_[i] < SMALL)
    {
        Bterm[i]= 0.0;
    }
    else
    {
        // diffusion resistance
        Bterm[i]   =   ((1-fracRedThird)/fracRedThird)*((radius_[i]/cg_)/(diffEff));
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::getMassT(int i)
{
    // initialize sherwood & schmidt numbers for every particle
    double Sc_i = 0.;
    double Sh_i = 0.;
    // if molecular diffusion is around 0, overwrite to avoid numerical errors.
    if (molecularDiffusion_[i] < SMALL)
        molecularDiffusion_[i] += SMALL;

    if (nuf_[i] == 0.0 || Rep_[i] == 0.0)
    {
        Massterm[i] = 0.0;
    }
    else
    {
        Sc_i  =   nuf_[i]/molecularDiffusion_[i];
        Sh_i  =   2.0+0.6*sqrt(Rep_[i])*cbrt(Sc_i);

        Massterm[i] = Sh_i*molecularDiffusion_[i]/(2.0*(radius_[i]/cg_)+SMALL);
        Massterm[i] = 1.0/(Massterm[i]);
    }

    if (screenflag_ && screen)
        fprintf(screen, "Schmidt number: %f \n",Sc_i);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::reaction(int i, double &dmA_, const double x0_eq_)
{

    double p_eq_ = 0.;

    p_eq_ = x0_eq_ * partP_[i];

    const double p_A = xA_[i] * partP_[i];

    if (screenflag_ && screen)
    {
        fprintf(screen, "p_eq_I: %f, p_A: %f \n", p_eq_,p_A);
    }

    // rate of chemical reaction for 1 active layer
    const double W = Aterm[i] + Bterm[i] + Massterm[i];

    dY[i]   =   (p_A - p_eq_) / W;

    if (dY[i] < 0.0)
        dY[i] = 0.0;


    // mass flow rate for reactant gas species
    // dmA is a positive value
    dmA_ =   dY[i] * (1.0 / (Runiv * T_[i])) * molMass_A_
               * (MY_4PI * ((radius_[i] * radius_[i]) / (cg_ * cg_)))
               * TimeStep * nevery;

    // limit mass change - can't remove more than present in cell
    // limit it with species mass per volume x voidfraction x cell volume / particles in cell x relaxation factor
    if(limit_reactant_consumption_)
    {
        if (screenflag_ && screen) fprintf(screen,"checking reactant limitation\n");

        double dAmax = p_A / (Runiv * T_[i]) * molMass_A_ * reactantPerParticle_[i] * maxReactantConsumptionFrac_;

        if(dmA_ > dAmax) dmA_ = dAmax;
    }

    // fix property added so values are outputted to file
    dmA_f_[i] = dmA_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::update_atom_properties(int i, const double dmA_)
{
    if (!shrink_) return;
    if (screenflag_ && screen)
        fprintf(screen,"run update atom props \n");

    // based on material change: update relative radii, average density and mass of atom i

    // initialize radius, mass change of layer and sum of particle
    double rad[2] = {0.};
    double dmL_[2] = {0.};
    double sum_mass_p_new = 0.0;
    double Cp = 0.0;
    double layer_Cp[2] = {0.};
    double Tpart = Tpart_[i];

    if (shrink_)
    {
        // Mass change of inner layer
        dmL_[1] = dmA_ * nu_B_ * (molMass_B_ / molMass_A_);

        // New layer mass
        massLayer_[i][1] -= dmL_[1]*scale_reduction_rate;
        if (massLayer_[i][1] < 0.0) massLayer_[i][1] = 0.0;

        dmL_[0] = dmL_[1] * molMass_D_ * nu_D_ / (molMass_B_ * nu_B_);
        massLayer_[i][0] += dmL_[0]*scale_reduction_rate;
    }

    layer_Cp[0] = spec_heat(a_coeff_ash,Tpart);
    layer_Cp[1] = spec_heat(a_coeff_coke,Tpart);
    Cp += massLayer_[i][0] * layer_Cp[0] / molMass_D_;
    if (layers_ > 0) Cp += massLayer_[i][1] * layer_Cp[1] / molMass_B_; // only include this contribution if core is still present

    for (int j = 0; j <= layers_; j++)
    {
        // calculate total mass of particle
        // since there is a minimum radius for layers, there is always a
        // non-zero contribution (at least from the innermost layer)
        sum_mass_p_new += massLayer_[i][j];
    }
    Cp /= sum_mass_p_new;
    fix_thermal_capacity_->vector_atom[i] = Cp;

    if (!shrink_) return;

    // if (screen) fprintf(screen,"total mass of particle = %f \n", sum_mass_p_new);

    // Total mass of particle with coarse-graining
    pmass_[i] = sum_mass_p_new*cg_*cg_*cg_;

    // if (screen) fprintf(screen, "pmass = %f \n",pmass_[i]);
    rad[1] = cbrt((0.75*massLayer_[i][1])/(rhoeff_[i][1]*M_PI));
    // NOTE: This keeps the outer radius of the particle constant. This induces a slight mass conservation error
    //       unless the inner (B) and outer (D) layer satisfy molMass_D_ * nu_D_ / (porosity_D_ * rho_D_) = molMass_B_ * nu_B_ / (porosity_B_ * rho_B_) 
    rad[0] = radius_[i]/cg_;

    if (fix_polydisp_)
    {
        rad[1] /= cbrt(effvolfactors_[i]);
    }

    // Determine new relative radii after reduction

    relRadii_[i][1] = rad[1]/rad[0];
    relRadii_[i][1] = std::min(0.9999, relRadii_[i][1]);

    // total particle effective density
    pdensity_[i] = 0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);
    if (fix_polydisp_)
    {
        pdensity_[i] /= effvolfactors_[i];
    }

    if (screenflag_ && screen)
        fprintf(screen, "radius_: %f, pmass_: %f, pdensity_: %f\n ", radius_[i], pmass_[i], pdensity_[i]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::update_gas_properties(int i, const double dmA_)
{
    /*double kch2_ = 0.0;
    kch2_ = xA_[i] + xC_[i];*/

    // based on material change: update gas-phase source terms for mass and heat
    double dmA = dmA_*cg_*cg_*cg_;
    // Reactant gas mass change
    changeOfA_[i]   -=  dmA;
    // Limit maximum reactant gas
    if (changeOfA_[i] > 0.0) changeOfA_[i] = 0.0;
    // Product gas mass change
    changeOfC_[i]   +=  dmA*nu_C_*molMass_C_/molMass_A_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::FractionalReduction(int i)
{
    const double f = 1.0 - relRadii_[i][1]*relRadii_[i][1]*relRadii_[i][1];

    fracRed_[i][0] = f;
}


/* ---------------------------------------------------------------------- */

/* heat of reaction calcualtion using the NASA 7 coefficient polynomials */
void FixChemShrinkCoreSingle::heat_of_reaction(int i, const double dmA_)
{

    if (!heatToParticle_ && !heatToFluid_) return;

    double T = T_[i];
    double Tpart = Tpart_[i];
    double dH = 0.0;


    // reaction enthalpy in J/mol
    if (reactionHeatIndex_ == 0)  // O2 + 2C -> 2C0
    {
        dH -= conv_enthalpy(a_coeff_O2,T);
        dH += 2 * conv_enthalpy(a_coeff_CO,T);
        dH -= 2 * conv_enthalpy(a_coeff_coke,Tpart);
        // neglect any contribution to ash
    }
    else if (reactionHeatIndex_ == 1) // CO2 + C -> 2CO
    {
        dH -= conv_enthalpy(a_coeff_CO2,T);
        dH += 2 * conv_enthalpy(a_coeff_CO,T);
        dH -= conv_enthalpy(a_coeff_coke,Tpart);
        // neglect any contribution to ash
    }

    double cg3 = cg_ * cg_ * cg_;

    if (heatToParticle_)
    {
        // heat flux = released heat / (nevery * dt), but it is reset to 0 in scalar transport equation
        // at the beginning of each step, hence an additional factor nevery which cancels that in the denominator 
        heatFlux_[i] -= dmA_ / molMass_A_ * dH / TimeStep * cg3;
    }
    else
    {
        // CFDEM takes accumulated heat, resets it after pull, and divides by time since last pull
        reactionheat_[i] -= dmA_ / molMass_A_ * dH * cg3;
    }
}
/* ---------------------------------------------------------------------- */

/* Calculate conventional enthalpies of species */

double FixChemShrinkCoreSingle::conv_enthalpy (const double *a, double Ti)
{
    double value = 0.;

    if (Ti < SMALL)
        error->fix_error(FLERR, this, "Error T <= ZERO");

    if (Ti < a[0]) { // Temperature smaller than lower bound
        const double Tbound_low = a[0];
        const double Tbound_low_sq = Tbound_low*Tbound_low;
        const double Tbound_low_cb = Tbound_low_sq*Tbound_low;
        value =   a[10]*Tbound_low
                + a[11]*Tbound_low_sq*0.5
                + a[12]*Tbound_low_cb/3.0
                + a[13]*Tbound_low_sq*Tbound_low_sq*0.25
                + a[14]*Tbound_low_sq*Tbound_low_cb*0.20
                + a[15];
    } else if (Ti < a[2]) {
        const double Ti_sq = Ti*Ti;
        const double Ti_cb = Ti_sq*Ti;
        value =   a[10]*Ti
                + a[11]*Ti_sq*0.5
                + a[12]*Ti_cb/3.0
                + a[13]*Ti_sq*Ti_sq*0.25
                + a[14]*Ti_sq*Ti_cb*0.20
                + a[15];
    } else if (Ti < a[1]) {
        const double Ti_sq = Ti*Ti;
        const double Ti_cb = Ti_sq*Ti;
        value =   a[3]*Ti
                + a[4]*Ti_sq*0.5
                + a[5]*Ti_cb/3.0
                + a[6]*Ti_sq*Ti_sq*0.25
                + a[7]*Ti_sq*Ti_cb*0.20
                + a[8];
    } else {
        const double Tbound_high = a[1];
        const double Tbound_high_sq = Tbound_high*Tbound_high;
        const double Tbound_high_cb = Tbound_high_sq*Tbound_high;
        value =   a[3]*Tbound_high
                + a[4]*Tbound_high_sq*0.5
                + a[5]*Tbound_high_cb/3.0
                + a[6]*Tbound_high_sq*Tbound_high_sq*0.25
                + a[7]*Tbound_high_sq*Tbound_high_cb*0.20
                + a[8];
    }

    return value*Runiv;
}

/* ---------------------------------------------------------------------- */


/* Calculate specific heat capacity of species [J / (mol K)] */

double FixChemShrinkCoreSingle::spec_heat (const double *a, double Ti)
{
    double value = 0.;

    if (Ti < SMALL)
        error->fix_error(FLERR, this, "Error T <= ZERO");

    if (Ti < a[0]) { // Temperature smaller than lower bound
        const double Tbound_low = a[0];
        const double Tbound_low_sq = Tbound_low*Tbound_low;
        const double Tbound_low_cb = Tbound_low_sq*Tbound_low;
        value =   a[10]
                + a[11]*Tbound_low
                + a[12]*Tbound_low_sq
                + a[13]*Tbound_low_cb
                + a[14]*Tbound_low_sq*Tbound_low_sq;
    } else if (Ti < a[2]) {
        const double Ti_sq = Ti*Ti;
        const double Ti_cb = Ti_sq*Ti;
        value =   a[10]
                + a[11]*Ti
                + a[12]*Ti_sq
                + a[13]*Ti_cb
                + a[14]*Ti_sq*Ti_sq;
    } else if (Ti < a[1]) {
        const double Ti_sq = Ti*Ti;
        const double Ti_cb = Ti_sq*Ti;
        value =   a[3]
                + a[4]*Ti
                + a[5]*Ti_sq
                + a[6]*Ti_cb
                + a[7]*Ti_sq*Ti_sq;
    } else {
        const double Tbound_high = a[1];
        const double Tbound_high_sq = Tbound_high*Tbound_high;
        const double Tbound_high_cb = Tbound_high_sq*Tbound_high;
        value =   a[3]
                + a[4]*Tbound_high
                + a[5]*Tbound_high_sq
                + a[6]*Tbound_high_cb
                + a[7]*Tbound_high_sq*Tbound_high_sq;
    }

    return value*Runiv;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::init_defaults()
{
    molMass_A_ = molMass_B_ = molMass_C_ = molMass_D_ = 0.0;
    rhoeff_ = NULL;
    porosity_ = NULL;
    pore_diameter_ = tortuosity_ = 0.0;
    relRadii_ = massLayer_ = NULL;
    xA_ = xC_ = NULL;
    scale_reduction_rate = 1.;
    layerDensities_ = NULL;

    nu_A_ = nu_B_ = nu_C_ = nu_D_ = 1;

    // particle properties total
    radius_ = pmass_ = pdensity_ = NULL;

    // initialise fix handles
    changeOfA_ = changeOfC_ = T_ = Tpart_ = molecularDiffusion_ = nuf_ = Rep_ = partP_ = Massterm = heatFlux_ = NULL;
    Aterm = Bterm = effDiffBinary = effDiffKnud = NULL;
    fracRed_ = NULL;

    dY = dmA_f_ = NULL;

    TimeStep = 0.0;

    fix_changeOfA_ = NULL;
    fix_changeOfC_ = NULL;
    fix_tgas_ = NULL;       // [K]
    fix_tpart_ = NULL;       // [K]
    fix_heatFlux_= NULL;
    fix_diffcoeff_ = NULL;  // [m^2/s]
    fix_nuField_ = NULL;    // [m^2/s]
    fix_partRe_ = NULL;

    fix_moleFractionA_ = NULL;
    fix_moleFractionC_ = NULL;

    fix_fracRed = NULL;
    fix_Aterm = NULL;
    fix_Bterm = NULL;
    fix_Massterm = NULL;
    fix_effDiffBinary = NULL;
    fix_effDiffKnud = NULL;
    fix_partPressure_ = NULL;   // Pascal
    fix_layerRelRad_ = NULL;
    fix_layerMass_ = NULL;      //  [kg]
    fix_thermal_capacity_ = NULL;
    fix_layerDens_ = NULL;      //  [kg/m^3]
    fix_porosity_ = NULL;       //  [%]
    fix_rhoeff_ = NULL;
    fix_tortuosity_ = NULL;
    fix_pore_diameter_ = NULL;  // [m]
    fix_dY_ = NULL;
    fix_dmA_ = NULL;

    fix_reactionheat = NULL;

    fix_reactantPerParticle_ = NULL;
    reactantPerParticle_ = NULL;
    limit_reactant_consumption_ = false;
    maxReactantConsumptionFrac_ = 0.5;

    massA = massC = NULL;
    diffA = moleFracA = moleFracC = NULL;
    speciesA = speciesC = NULL;
}

void FixChemShrinkCoreSingle::update_fix(int narg, char **arg)
{
    setup(0);
    int nlocal = atom->nlocal;
    int *mask   =   atom->mask;
    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            pmass_[i] = (massLayer_[i][0] + massLayer_[i][1])*cg_*cg_*cg_;
            pdensity_[i] = 0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);
            if (fix_polydisp_)
            {
                pdensity_[i] /= effvolfactors_[i];
            }

            double Cp = 0.0;
            double layer_Cp[2] = {0.};
            double Tpart = Tpart_[i];
            layer_Cp[0] = spec_heat(a_coeff_ash,Tpart);
            layer_Cp[1] = spec_heat(a_coeff_coke,Tpart);
            Cp += massLayer_[i][0] * layer_Cp[0] / molMass_D_;
            Cp += massLayer_[i][1] * layer_Cp[1] / molMass_B_;
            Cp /= (massLayer_[i][0] + massLayer_[i][1]);
            fix_thermal_capacity_->vector_atom[i] = Cp;
        }
    }
}
