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
#include "fix_chem_shrink_core.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_property_atom_polydispparcel.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL   1e-10

const double FixChemShrinkCore::Runiv = 8.3144;

// 7-coefficient NASA polynomials, see
// Burcat, A., Ruscic, B. (2005). Third millenium ideal gas and condensed phase thermochemical database for combustion
// (with update from active thermochemical tables) (No. ANL-05/20).
// T_low, T_high, T_middle,
// The first  set of 7 constants belongs to the T_middle - T_high K polynomial,
// the second set of 7 constants belongs to the T_low  - T_middle K polynomial
const double FixChemShrinkCore::a_coeff_nasa_Fe2O3[] = { 298.15, 1700., 1000.,
                                                        -3.240851E+02,  1.152686E+00, -1.413588E-03,  7.496435E-07, -1.455880E-10, -2.647718E+04,  1.609668E+03,
                                                         1.066786E+01, -4.869774E-03,  5.056287E-05, -4.500105E-08,  7.709213E-12, -1.025006E+05, -5.068585E+01 };
const double FixChemShrinkCore::a_coeff_nasa_Fe3O4[] = { 298.15, 1870., 1000.,
                                                         1.948880E+02, -4.225771E-01,  3.843026E-04, -1.536471E-07,  2.301151E-11, -1.944678E+05, -1.023246E+03,
                                                         7.781798E+01, -4.602735E-01,  1.199812E-03, -1.203296E-06,  4.119178E-10, -1.457550E+05, -3.319564E+02 };
const double FixChemShrinkCore::a_coeff_nasa_FeO[]   = { 298.15, 1650., 1000.,
                                                         5.762730E+00,  1.442345E-03,  1.309698E-07, -2.476342E-10,  4.756658E-14, -3.453377E+04, -2.603715E+01,
                                                         5.108889E+00,  3.732856E-03, -2.817518E-06,  1.393173E-09, -2.814222E-13, -3.438676E+04, -2.280153E+01 };
const double FixChemShrinkCore::a_coeff_nasa_Fe[]    = { 298.15, 1809., 1000.,
                                                        -2.674281E+02,  8.564543E-01, -9.799567E-04,  4.844980E-07, -8.760776E-11,  6.533737E+04,  1.349404E+03,
                                                        -8.123461E+00,  8.207546E-02, -2.126480E-04,  2.357984E-07, -9.114276E-11,  3.344653E+02,  3.269875E+01 };
const double FixChemShrinkCore::a_coeff_nasa_CO[]    = { 200.00, 6000., 1000.,
                                                         3.048486E+00,  1.351728E-03, -4.857941E-07,  7.885364E-11, -4.698075E-15, -1.426612E+04,  6.017098E+00,
                                                         3.579534E+00, -6.103537E-04,  1.016814E-06,  9.070059E-10, -9.044245E-13, -1.434409E+04,  3.508409E+00 };
const double FixChemShrinkCore::a_coeff_nasa_CO2[]   = { 200.00, 6000., 1000.,
                                                         4.636511E+00,  2.741457E-03, -9.958976E-07,  1.603867E-10, -9.161986E-15, -4.902490E+04, -1.934896E+00,
                                                         2.356813E+00,  8.984130E-03, -7.122063E-06,  2.457301E-09, -1.428855E-13, -4.837197E+04,  9.900904E+00 };
const double FixChemShrinkCore::a_coeff_nasa_H2[]    = { 200.00, 6000., 1000.,
                                                         2.932831E+00,  8.265980E-04, -1.464006E-07,  1.540985E-11, -6.887962E-16, -8.130558E+02, -1.024316E+00,
                                                         2.344303E+00,  7.980425E-03, -1.947792E-05,  2.015697E-08, -7.376029E-12, -9.179241E+02,  6.830022E-01 };
const double FixChemShrinkCore::a_coeff_nasa_H2O[]   = { 200.00, 6000., 1000.,
                                                         2.677039E+00,  2.973182E-03, -7.737689E-07,  9.443351E-11, -4.268999E-15, -2.988589E+04,  6.882550E+00,
                                                         4.198635E+00, -2.036402E-03,  6.520342E-06, -5.487927E-09,  1.771968E-12, -3.029373E+04, -8.490090E-01 };

const double FixChemShrinkCore::v_reac_[] = { 1.0, 1.0, 3.0 }; // reaction of w, m, h
const double FixChemShrinkCore::v_prod_[] = { 1.0, 3.0, 2.0 }; // production of Fe, w, m

#define SWITCH_LOW_HIGH_TEMPERATURE 843.15

#ifdef TWO_LAYERS
#define MAX_LAYERS 2
const double FixChemShrinkCore::v_reac_low_[] = { 0.25, 3.0, 0.0 };
const double FixChemShrinkCore::v_prod_low_[] = { 0.75, 2.0, 0.0 };
const double FixChemShrinkCore::k0_low_CO[] = { 150., 150. };
const double FixChemShrinkCore::k0_low_H2[] = {  50.,  25. };
const double FixChemShrinkCore::Ea_low_CO[] = { 70000., 75000. };
const double FixChemShrinkCore::Ea_low_H2[] = { 75000., 75000. };
//                                                  {       Fe,    Fe3O4,     Fe2O3 }
const double FixChemShrinkCore::layerMolMasses_[] = { 0.055845, 0.231532, 0.1596882 };
#else
#define MAX_LAYERS 3
#define PSEUDO_THREE_LAYERS
#ifdef PSEUDO_THREE_LAYERS // ignore wustite layer; TODO! beware of CO transformation!!!
const double FixChemShrinkCore::v_reac_low_[] = { 0.0, 0.25, 3.0 }; // reaction of Fe, m, h
const double FixChemShrinkCore::v_prod_low_[] = { 0.75, 0.0, 2.0 }; // production of Fe, Fe, m
const double FixChemShrinkCore::k0_low_CO[] = { 100., 100., 150. };
const double FixChemShrinkCore::k0_low_H2[] = {  33.,  33.,  25. };
const double FixChemShrinkCore::Ea_low_CO[] = { 70000., 70000., 75000. };
const double FixChemShrinkCore::Ea_low_H2[] = { 70000., 70000., 75000. };
#else

#endif
//                                                  {       Fe,      FeO,    Fe3O4,     Fe2O3 }
const double FixChemShrinkCore::layerMolMasses_[] = { 0.055845, 0.071844, 0.231532, 0.1596882 }; // [kg/mol]
#endif

enum {
    LAYER_WUSTITE = 0,
    LAYER_MAGNETITE,
    LAYER_HEMATITE,
    LAYER_MAX
};

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    layers_(MAX_LAYERS),
    minMolarFrac_(1e-3),
    rmin_(1e-5),      //  [m]
    created_fix_layerMass_(false),
    created_fix_rhoeff_(false),
    created_fix_fracRed(false),
    created_fix_internal_energy_(false)
{
    if ((strncmp(style, "chem/shrink/core", 15) == 0) && ((!atom->radius_flag) || (!atom->rmass_flag)))
        error->all(FLERR, "Fix chem/shrink/core needs per particle radius and mass");

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
            if (screenflag_ && screen) fprintf(screen, "spreciesA: %s", speciesA);
            iarg_ += 2; // iarg  = 5
            hasargs = true;
        }
        else if (strcmp(arg[iarg_], "molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_A_ = atof(arg[iarg_ + 1]);
            if (molMass_A_ < 0.0)
                error->fix_error(FLERR, this, "molar mass of A is not defined");
            iarg_ += 2; // iarg = 7
            hasargs = true;
        }
        else if (strcmp(arg[iarg_], "speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            speciesC = new char[strlen(arg[iarg_+1])+1];
            strcpy(speciesC, arg[iarg_+1]);
            if (screenflag_ && screen) fprintf(screen, "spreciesC: %s", speciesC);
            iarg_ += 2; // iarg = 9
            hasargs = true;
        }
        else if (strcmp(arg[iarg_], "molMassC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_C_ = atof(arg[iarg_ + 1]);
            if (molMass_C_ < 0.0)
                error->fix_error(FLERR, this, "molar mass of C is not defined");
            iarg_ += 2; // iarg = 11
            hasargs = true;
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
            hasargs = true; //iarg = 13
        }
        else if (strcmp(arg[iarg_],"screen") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) screenflag_ = 1;
            else if (strcmp(arg[iarg_+1],"no") == 0) screenflag_ = 0;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2; // iarg = 15
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

void FixChemShrinkCore::post_create()
{
    if (fix_Aterm == NULL) {
        char *fixname = new char [strlen("Aterm_")+strlen(id)+1];
        strcpy(fixname,"Aterm_");
        strcat(fixname,id);

        const char* fixarg[11];
        // register property/atom for chemical reaction resistance
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";          // vector with 3 values (for the different layers)
        fixarg[5]="yes";             // restart yes
        fixarg[6]="no";              // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";             // tale 0 as default value
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_Aterm = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete [] fixname;
    }

    if (fix_Bterm == NULL) {
        char *fixname = new char[strlen("Bterm_")+strlen(id)+1];
        strcpy(fixname,"Bterm_");
        strcat(fixname,id);

        const char* fixarg[11];
        // register property/atom for diffusion resistance term
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";          // vector with 3 values (for the different layers)
        fixarg[5]="yes";             // restart yes
        fixarg[6]="no";              // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_Bterm = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
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

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup]; //"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_effDiffBinary = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete [] fixname;
    }

    if (fix_effDiffKnud == NULL) {
        char *fixname = new char[strlen("effDiffKnud_")+strlen(id)+1];
        strcpy(fixname,"effDiffKnud_");
        strcat(fixname,id);

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_effDiffKnud = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_dY_ == NULL) {
        char* fixname = new char[strlen("dY_")+strlen(id)+1];
        strcpy(fixname,"dY_");
        strcat(fixname,id);

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_dY_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    if (fix_dmA_ == NULL) {
        char* fixname = new char[strlen("dmA_")+strlen(id)+1];
        strcpy(fixname,"dmA_");
        strcat(fixname,id);

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]=group->names[igroup];//"all";
        fixarg[2]="property/atom";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="yes";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_dmA_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete []fixname;
    }

    // shared fixes
    fix_layerMass_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer", "property/atom", "vector", MAX_LAYERS+1, 0, style, false));
    if (fix_layerMass_ == NULL) {
        const char* fixarg[12];
        fixarg[0]="LayerMasses";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="massLayer";
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fixarg[11]="0.0";
        fix_layerMass_ = modify->add_fix_property_atom(12,const_cast<char**>(fixarg),style);
        created_fix_layerMass_ = true;
    }

    fix_rhoeff_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("rhoeff", "property/atom", "vector", MAX_LAYERS+1, 0, style, false));
    if (fix_rhoeff_ == NULL) {
        const char* fixarg[12];
        fixarg[0]="rhoeff";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="rhoeff";
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fixarg[11]="0.0";
        fix_rhoeff_ = modify->add_fix_property_atom(12,const_cast<char**>(fixarg),style);
        created_fix_rhoeff_ = true;
    }

    fix_fracRed = static_cast<FixPropertyAtom*>(modify->find_fix_property("fracRed", "property/atom", "vector", MAX_LAYERS, 0, style, false));
    if (fix_fracRed == NULL) {
        const char* fixarg[11];
        fixarg[0]="fracRed";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="fracRed";
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.0";
        fixarg[9]="0.0";
        fixarg[10]="0.0";
        fix_fracRed = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        created_fix_fracRed = true;
    }

    fix_internal_energy_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("internalEnergy","property/atom","scalar",0,0,style,false));
    if (fix_internal_energy_ == NULL) {
        const char* fixarg[9];
        fixarg[0]="internalEnergy";            // fixid
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="internalEnergy";           // propertyid
        fixarg[4]="scalar";
        fixarg[5]="yes";             // restart yes
        fixarg[6]="yes";             // communicate ghost - yes
        fixarg[7]="no";              // communicate rev no
        fixarg[8]="0.0";
        fix_internal_energy_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
        created_fix_internal_energy_ = true;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
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
        if (fix_layerMass_ && created_fix_layerMass_)   { modify->delete_fix(fix_layerMass_->id); massLayer_ = NULL; }
        if (fix_rhoeff_    && created_fix_rhoeff_)      { modify->delete_fix(fix_rhoeff_->id); rhoeff_ = NULL; }
        if (fix_fracRed    && created_fix_fracRed)      { modify->delete_fix(fix_fracRed->id); fracRed_ = NULL; }
        if (fix_internal_energy_ && created_fix_internal_energy_) {modify->delete_fix(fix_internal_energy_->id);}
    }
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
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

int FixChemShrinkCore::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::updatePtrs()
{
    changeOfA_      =   fix_changeOfA_      ->  vector_atom;
    changeOfC_      =   fix_changeOfC_      ->  vector_atom;
    T_              =   fix_tgas_           ->  vector_atom;
    Tpart_          =   fix_tpart_          ->  vector_atom;
    reactionHeat_   =   fix_reactionHeat_   ->  vector_atom;
    molecularDiffusion_  = fix_diffcoeff_   ->  vector_atom;
    nuf_            =   fix_nuField_        ->  vector_atom;
    Rep_            =   fix_partRe_         ->  vector_atom;
    xA_             =   fix_moleFractionA_  ->  vector_atom;
    xC_             =   fix_moleFractionC_  ->  vector_atom;
    relRadii_       =   fix_layerRelRad_    ->  array_atom;
    massLayer_      =   fix_layerMass_      ->  array_atom;
    if (fix_polydisp_)
        effvolfactors_  =   fix_polydisp_   ->  vector_atom;

#ifdef PER_ATOM_LAYER_DENSITIES
    layerDensities_ =   fix_layerDens_      ->  array_atom;
#else
    layerDensities_ =   fix_layerDens_      ->  values;
#endif

    k0_             =   fix_k0_             ->  values;
    Ea_             =   fix_Ea_             ->  values;
    rhoeff_         =   fix_rhoeff_         ->  array_atom;
    //
    porosity_       =   fix_porosity_       ->  values;
    tortuosity_     =   fix_tortuosity_     ->  compute_scalar();
    pore_diameter_  =   fix_pore_diameter_  ->  values;
    //
    fracRed_        =   fix_fracRed         ->  array_atom;
    Aterm           =   fix_Aterm           ->  array_atom;
    Bterm           =   fix_Bterm           ->  array_atom;
    Massterm        =   fix_Massterm        ->  vector_atom;
    effDiffBinary   =   fix_effDiffBinary   ->  array_atom;
    effDiffKnud     =   fix_effDiffKnud     ->  array_atom;
    partP_          =   fix_partPressure_   ->  vector_atom;

    TimeStep        =   update  -> dt;
    radius_         =   atom    -> radius;
    pmass_          =   atom    -> rmass;
    pdensity_       =   atom    -> density;

    dY          =   fix_dY_         ->  array_atom;
    dmA_f_      =   fix_dmA_ -> array_atom;

    if(limit_reactant_consumption_)
    {
        reactantPerParticle_ = fix_reactantPerParticle_ -> vector_atom;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::init()
{
    // error checks
    if (!atom->radius_flag)
      error->fix_error(FLERR,this,"requires atom attribute radius (per-particle)");
    if (!atom->rmass_flag)
      error->fix_error(FLERR,this,"requires atom attribute mass (per-particle)");
    if (!atom->tag_enable || 0 == atom->map_style)
      error->fix_error(FLERR,this,"requires atom tags and an atom map");

    char* propertyname = new char[strlen("k0_")+strlen(id)+1];
    strcpy (propertyname,"k0_");
    strcat(propertyname,id);
    fix_k0_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname,"property/global","vector",MAX_LAYERS,0,style));
    delete [] propertyname;

    // look up activation energies Ea
    // differs for every ore id
    propertyname = new char [strlen("Ea_")+strlen(id)+1];
    strcpy(propertyname, "Ea_");
    strcat(propertyname, id);
    fix_Ea_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname,"property/global","vector",MAX_LAYERS,0,style));
    delete [] propertyname;

    propertyname = new char [strlen("density_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"density_");
    strcat(propertyname,group->names[igroup]);
#ifdef PER_ATOM_LAYER_DENSITIES
    fix_layerDens_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname,"property/atom","vector",0,0,style));
#else
    fix_layerDens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname,"property/global","vector",MAX_LAYERS+1,0,style));
#endif
    delete [] propertyname;

    // references for per atom properties.
    fix_changeOfA_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massA, "property/atom", "scalar", 0, 0, style));
    fix_changeOfC_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massC, "property/atom", "scalar", 0, 0, style));
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, style));
    fix_tpart_          =   static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp", "property/atom", "scalar", 0, 0, style));     // particle temperature
    fix_reactionHeat_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("reactionHeat", "property/atom", "scalar", 0, 0, style));
    fix_diffcoeff_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffA, "property/atom", "scalar", 0, 0, style));
    fix_nuField_        =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partNu", "property/atom", "scalar", 0, 0, style));
    fix_partRe_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partRe", "property/atom", "scalar", 0, 0, style));
    fix_moleFractionA_  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(moleFracA, "property/atom", "scalar", 0, 0, style));
    fix_moleFractionC_  =   static_cast<FixPropertyAtom*>(modify->find_fix_property(moleFracC, "property/atom", "scalar", 0, 0, style));

    fix_partPressure_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partP", "property/atom", "scalar", 0, 0, style));
    fix_layerRelRad_    =   static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii", "property/atom", "vector", 0, 0, style));
    fix_layerMass_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer","property/atom","vector",0,0,style));
    fix_rhoeff_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("rhoeff", "property/atom", "vector", 0, 0, style));
    fix_thermal_capacity_=  static_cast<FixPropertyAtom*>(modify->find_fix_property("thermalCapacity", "property/atom", "scalar", 0, 0, style,false));
    fix_polydisp_       =   static_cast<FixPropertyAtomPolydispParcel*>(modify->find_fix_property("effvolfactor", "property/atom","scalar",0,0,style,false));

    // references for global properties - valid for every particle equally
    propertyname = new char [strlen("porosity_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"porosity_");
    strcat(propertyname,group->names[igroup]);
    fix_porosity_       =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "vector", 0, 0, style));
    delete [] propertyname;

    propertyname = new char [strlen("tortuosity_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"tortuosity_");
    strcat(propertyname,group->names[igroup]);
    fix_tortuosity_     =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char [strlen("pore_diameter_")+strlen(group->names[igroup])+1];
    strcpy(propertyname,"pore_diameter_");
    strcat(propertyname,group->names[igroup]);
    fix_pore_diameter_  =   static_cast<FixPropertyGlobal*>(modify->find_fix_property(propertyname, "property/global", "vector", 0, 0,style));
    delete [] propertyname;

    propertyname = new char [strlen("Aterm_")+strlen(id)+1];
    strcpy (propertyname,"Aterm_");
    strcat(propertyname,id);
    fix_Aterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("Bterm_")+strlen(id)+1];
    strcpy (propertyname,"Bterm_");
    strcat(propertyname,id);
    fix_Bterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", MAX_LAYERS, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("Massterm_")+strlen(id)+1];
    strcpy (propertyname,"Massterm_");
    strcat(propertyname,id);
    fix_Massterm        =   static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "scalar", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("effDiffBinary_")+strlen(id)+1];
    strcpy (propertyname,"effDiffBinary_");
    strcat(propertyname,id);
    fix_effDiffBinary = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", MAX_LAYERS, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("effDiffKnud_")+strlen(id)+1];
    strcpy (propertyname,"effDiffKnud_");
    strcat(propertyname,id);
    fix_effDiffKnud = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", MAX_LAYERS, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("dY_")+strlen(id)+1];
    strcpy (propertyname,"dY_");
    strcat(propertyname,id);
    fix_dY_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", 0, 0, style));
    delete [] propertyname;

    propertyname = new char[strlen("dmA_")+strlen(id)+1];
    strcpy (propertyname,"dmA_");
    strcat(propertyname,id);
    fix_dmA_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(propertyname, "property/atom", "vector", 0, 0, style));
    delete [] propertyname;

    fix_fracRed         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("fracRed", "property/atom", "vector", 0, 0, style));

    if(limit_reactant_consumption_)
    {
        fix_reactantPerParticle_ = static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactantPerParticle","property/atom","scalar",0,0,style));
    }

    if(fix_thermal_capacity_)
    {
         variableCp_ = true;
    }
    else
    {
         variableCp_ = false;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::setup(int)
{
    updatePtrs();

    int nlocal = atom->nlocal;
    int *mask  = atom->mask;

    // set initial values for rhoeff, and use them to calculate mass of layers
    for (int i = 0; i < nlocal; ++i)
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

void FixChemShrinkCore::pre_neighbor()
{
    setup(0);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    if (update->ntimestep % nevery)
        return;

    updatePtrs();
    int i = 0;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    double x0_eq_[MAX_LAYERS] = {0.}; // molar fraction of reactant gas
    double dmA_[MAX_LAYERS] = {0.};   // mass flow rate of reactant gas species for each layer at w->fe, m->w & h->m interfaces

    for (i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (xA_[i] < minMolarFrac_)
            {
                continue;
            }

            // 1st recalculate masses of layers if layer has reduced
            // is ignored if there is no change in layers
            if (active_layers(i) > 0)
            {
                if (T_[i] < 573.15) // Kelvin / 300 Celsius
                {
                    // do nothing -- no reaction takes place
                    if (screenflag_) error->warning(FLERR, "The temperature is too low for reduction to take place!");
                }
                else if (T_[i] < SWITCH_LOW_HIGH_TEMPERATURE) // T_[i] between 573.15 and 843.15 Kelvin (300 and 570 Celsius)
                {
#ifdef TWO_LAYERS
                    FractionalReduction_low(i);
#else
#ifdef PSEUDO_THREE_LAYERS
                    // treating wustite layer as if it were Fe, thus no more reaction at this point
                    if (layers_ <= 1)
                        continue;
#endif

                    FractionalReduction(i);
#endif
                    getXi_low(i,x0_eq_);
                    getA_low(i);
                    getB(i);
                    getMassT(i);
                    reaction_low(i, dmA_, x0_eq_);
                    update_atom_properties(i, dmA_,v_reac_low_,v_prod_low_);
                    update_gas_properties(i, dmA_);
                    heat_of_reaction(i, dmA_,v_reac_low_,v_prod_low_);
                }
                else // T_[i] > 843.15 K
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
                    update_atom_properties(i, dmA_,v_reac_,v_prod_);
                    // also the results of reaction function is used to calculate
                    // the changes in gas species
                    update_gas_properties(i, dmA_);
                    // calculate delta_h, and dot_delta_h for heat of reaction
                    heat_of_reaction(i, dmA_,v_reac_,v_prod_);
                }
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCore::active_layers(int i)
{
    layers_ = MAX_LAYERS;

    for (int j = layers_; j > 0; j--)
    {
        if (relRadii_[i][j]*(radius_[i]/cg_) < rmin_)
        {
            --layers_;
        }
    }

    if (screenflag_ && screen)
        fprintf(screen, "active layers: %i \n", layers_);
    return layers_;
}

/* ---------------------------------------------------------------------- */

// 0 = iron shell, 1 = wüstite layer, 2 = magnetite layer, 3 = hematite layer
void FixChemShrinkCore::calcMassLayer(int i)
{
    double rad[MAX_LAYERS+1] = {0.};
    for (int layer = 0; layer <= layers_ ; layer++)
        rad[layer] = (radius_[i]/cg_)*relRadii_[i][layer];

    massLayer_[i][layers_] = MY_4PI3*rhoeff_[i][layers_]*rad[layers_]*rad[layers_]*rad[layers_];

    for (int layer = 0 ; layer < layers_; layer++)
    {
        massLayer_[i][layer] = MY_4PI3*rhoeff_[i][layer]*(rad[layer  ]*rad[layer  ]*rad[layer  ]
                                                         -rad[layer+1]*rad[layer+1]*rad[layer+1]);
    }

    if (fix_polydisp_)
    {
        for (int layer = 0 ; layer <= layers_; layer++)
        {
            massLayer_[i][layer] *= effvolfactors_[i];
        }
    }
}

/* ---------------------------------------------------------------------- */

double FixChemShrinkCore::K_eq(int layer, int i)
{
    // 0 = FeO (wüstite)     -> Fe (iron)
    // 1 = Fe3O4 (magnetite) -> FeO (wüstite)
    // 2 = Fe2O3 (hematite)  -> Fe3O4 (magnetite)
    double Keq_ = 0.;

    if (strcmp(speciesA, "CO") == 0)
    {
        if (layer == 2)
            Keq_ = exp(3968.37/T_[i]+3.94);        // Valipour 2009, Nietrost 2012
        else if (layer == 1)
            Keq_ = pow(10.0,(-1834.0/T_[i]+2.17)); // Nietrost 2012
            // Keq_ = exp(-3585.64/T_[i]+4.58);    // Valipour 2006 // exp(-3585.64/T_[i]+8.98); // Valipour 2009
        else if (layer == 0)
            Keq_ = pow(10.0,(914.0/T_[i]-1.097));  // Nietrost 2012
            // Keq_ = exp(2744.63/T_[i]-2.946);    // Valipour 2009
     }
     else if(strcmp(speciesA,"H2")==0)
     {
        if (layer == 2)
            Keq_   =   exp(-362.6/T_[i] + 10.334); // Valipour 2009, Nietrost 2012
        else if (layer == 1)
            Keq_    =   pow(10.0,(-3577.0/T_[i]+3.74)); // Nietrost 2012
            // Keq_   =   exp(-7916.6/T_[i] + 8.46); // Valipour 2009
        else if (layer == 0)
            // Keq_    =   pow(10.0,(-827.0/T_[i]+0.468)); // Nietrost --> With this equilibrium constant R1 is occuring with all temperatures for H2
            Keq_    =   pow(10.0,(-856.66/T_[i]+0.4387)); // Equilibrium constant from Turkdogan "Physical Chemistry of High Temperature Technology"
            // Keq_   =   exp(-1586.9/T_[i] + 0.9317); // Valipour 2009
     }
     else
     {
         error->one(FLERR, "Undefined Reaction\n");
     }

     if(screenflag_ && screen)
         fprintf(screen,"Keq_ : %f \n",Keq_);

    return Keq_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getXi(int i, double *x0_eq_)
{
    const double kch2_ = xA_[i] + xC_[i];

    for (int j = 0; j < layers_; j++)
    {
        x0_eq_[j]  =   kch2_/(1.0+K_eq(j,i));
    }

    if (screenflag_ && screen)
        fprintf(screen, "x0_eq_0: %f, x0_eq_1: %f, x0_eq_2: %f \n", x0_eq_[0], x0_eq_[1], x0_eq_[2]);
}

/* ---------------------------------------------------------------------- */

// calculate A_[j] [s/m] - the chemical reaction resistance term
// Equation available in literature. (Valipour, Natsui, Nietrost...)
// 0 = wüstite interface, 1 = magnetite interface, 2 = hematite interface
void FixChemShrinkCore::getA(int i)
{
    for (int j = 0; j < layers_ ; j++)
    {
            Aterm[i][j] = (k0_[j] * exp(-Ea_[j] / (Runiv*T_[i])))
                        * cbrt(square(1.0 - fracRed_[i][j]))
                        * (1.0 + 1.0 / K_eq(j,i));
            Aterm[i][j] = 1.0 / Aterm[i][j];
    }
}

/* ---------------------------------------------------------------------- */

// Calculate B [s/m] - the diffusion resistance term
// Use binary diffusion for mixture, and knudsen diffusion to determine the effective diffusion term
// 0 : diffusion through iron layer, 1 : diffusion through wüstite, 2 : diffusion through magnetite.
// there is no diffusion through the hematite layer
void FixChemShrinkCore::getB(int i)
{
    double fracRedThird_[MAX_LAYERS] = {0.};
    double diffEff_[MAX_LAYERS] = {0.};

    for (int layer = 0; layer < layers_; layer++)
    {
        // calculate fractional reduction to the power of 1/3 for simpler use
        fracRedThird_[layer] = cbrt(1.0-fracRed_[i][layer]);

        // Calculate the effective molecular diffusion
        effDiffBinary[i][layer] = molecularDiffusion_[i]*(porosity_[layer]/tortuosity_) + SMALL;

        // Calculate the knudsen diffusion
        effDiffKnud[i][layer]  =  (pore_diameter_[layer]/6.0)*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_))*(porosity_[layer]/tortuosity_) + SMALL;

        // total effective diffusivity
        diffEff_[layer] =   effDiffKnud[i][layer]*effDiffBinary[i][layer]/(effDiffBinary[i][layer]+effDiffKnud[i][layer]) + SMALL;
    }

    if (screenflag_ && screen)
    {
#ifdef TWO_LAYERS
        fprintf(screen, "diffEff_[0]: %f, diffEff_[1]: %f \n"
                    "fracRedThird_[0]: %f, fracRedThird_[1]: %f \n"
                    "fracRed_[0]: %f, fracRed_[1]: %f \n",
                    diffEff_[0], diffEff_[1],
                    fracRedThird_[0], fracRedThird_[1],
                    fracRed_[0][0], fracRed_[0][1]);
#else
        fprintf(screen, "diffEff_[0]: %f, diffEff_[1]: %f, diffEff_[2]: %f \n"
                    "fracRedThird_[0]: %f, fracRedThird_[1]: %f, fracRedThird_[2] : %f \n"
                    "fracRed_[0]: %f, fracRed_[1]: %f, fracRed_[2]: %f \n",
                    diffEff_[0], diffEff_[1], diffEff_[2],
                    fracRedThird_[0], fracRedThird_[1] , fracRedThird_[2],
                    fracRed_[0][0], fracRed_[0][1], fracRed_[0][2]);
#endif
    }

    // calculation of diffusion term
    if (molecularDiffusion_[i] < SMALL)
    {
        for (int layer = 0; layer < layers_; layer++)
            Bterm[i][layer] = 0.0;
    }
    else
    {
        // diffusion resistance through Fe
        Bterm[i][0]   =   ((1-fracRedThird_[0])/fracRedThird_[0])*((radius_[i]/cg_)/(diffEff_[0]));
        for (int layer = 1; layer < layers_; layer++)
        {
            Bterm[i][layer] = (fracRedThird_[layer-1]-fracRedThird_[layer])/(fracRedThird_[layer-1]*fracRedThird_[layer])*((radius_[i]/cg_)/(diffEff_[layer]));
        }
    }

#ifndef TWO_LAYER
    if (T_[i] < SWITCH_LOW_HIGH_TEMPERATURE) {
#ifdef PSEUDO_THREE_LAYERS
        Bterm[i][0] = 1.0; // not used, set to something that does not disturb plotting of Bterms
        if (layers_ > 1)
            Bterm[i][1] = ((1-fracRedThird_[1])/fracRedThird_[1])*((radius_[i]/cg_)/(diffEff_[0]));
#else

#endif
    }
#endif

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getMassT(int i)
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

void FixChemShrinkCore::reaction(int i, double *dmA_, const double *x0_eq_)
{

#ifdef TWO_LAYERS
    return;
#endif

    double p_eq_[MAX_LAYERS] = {0.};

    for (int layer = 0; layer < layers_; layer++)
    {
        p_eq_[layer] = x0_eq_[layer] * partP_[i];
    }

    const double p_A = xA_[i] * partP_[i];

    if (screenflag_ && screen)
    {
        fprintf(screen, "p_eq_I: %f, p_eq_II: %f, p_eq_III: %f, p_A: %f \n", p_eq_[0], p_eq_[1],p_eq_[2],p_A);
    }

    if (layers_ == 3)
    {
        const double A0 = Aterm[i][0];
        const double A1 = Aterm[i][1];
        const double B1 = Bterm[i][1];
        const double A1plusB1         = A1 + B1;
        const double A2plusB2         = Aterm[i][2] + Bterm[i][2];
        const double B0plusMass       = Bterm[i][0] + Massterm[i];
        const double B0plusB1plusMass = Bterm[i][0] + B1 + Massterm[i];

        // including reaction resistance and diffusion coeff terms
        const double W = A2plusB2 * (A0 * (A1 + B0plusB1plusMass) + A1plusB1 * B0plusMass)
                       + A1       * (A0 * (     B0plusB1plusMass) +       B1 * B0plusMass);

        // hematite to magnetite
        dY[i][2] = ((A0 * (A1 + B0plusB1plusMass) + A1plusB1 * B0plusMass) * (p_A - p_eq_[2])
                  - (A0 * (     B0plusB1plusMass) +       B1 * B0plusMass) * (p_A - p_eq_[1])
                  - (                               A1       * B0plusMass) * (p_A - p_eq_[0])) / W;

        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][2] < 0.0)
            dY[i][2] = 0.0;

        // magnetite to wustite
        dY[i][1] = (((A2plusB2 + B1) * (A0 + B0plusMass) + A0 * B0plusMass) * (p_A - p_eq_[1])
                  - (            B1  * (A0 + B0plusMass) + A0 * B0plusMass) * (p_A - p_eq_[2])
                  - ( A2plusB2       *       B0plusMass)                    * (p_A - p_eq_[0])) / W;

        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][1] < 0.0)
            dY[i][1] = 0.0;

        // wustite to iron
        // if magnetite is not reducing, wustite does also not reduce
        if (dY[i][1] == 0.0)
        {
            dY[i][0] = 0.0;
        }
        else
        {
            dY[i][0] = ((A2plusB2 * (A1 + B0plusB1plusMass) + A1 * B0plusB1plusMass) * (p_A - p_eq_[0])
                    -   (                                     A1 * B0plusMass)       * (p_A - p_eq_[2])
                    -   (A2plusB2 *       B0plusMass)                                * (p_A - p_eq_[1])) / W;
        }

        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][0] < 0.0)
            dY[i][0] = 0.0;

    }
    else if (layers_ == 2)
    {
        const double A0 = Aterm[i][0];
        const double A1plusB1         = Aterm[i][1] + Bterm[i][1];
        const double B0plusMass       = Bterm[i][0] + Massterm[i];
        const double A0plusB0plusMass = A0 + B0plusMass;

        const double W = A1plusB1 * A0plusB0plusMass + A0 * B0plusMass;

        // hematite to magnetite
        dY[i][2] = 0.0;

        // magnetite to wustite
        dY[i][1]     = (A0plusB0plusMass * (p_A - p_eq_[1])  -  B0plusMass * (p_A - p_eq_[0])) / W;

        if (dY[i][1] < 0.0)
            dY[i][1] = 0.0;

        // wustite to iron
        if (dY[i][1] == 0.0)
        {
            dY[i][0] = 0.0;
        }
        else
        {
            dY[i][0] = ((A1plusB1 + B0plusMass) * (p_A - p_eq_[0])  -  B0plusMass * (p_A - p_eq_[1])) / W;
        }

        if (dY[i][0] < 0.0)
            dY[i][0] = 0.0;

    }
    else if (layers_ == 1)
    {
        // rate of chemical reaction for 1 active layer
        const double W = Aterm[i][0] + Bterm[i][0] + Massterm[i];

        // hematite to magnetite
        dY[i][2]   =   0.0;

        // magnetite to wustite
        dY[i][1]   =   0.0;

        // wustite to iron
        dY[i][0]   =   (p_A - p_eq_[0]) / W;

        if (dY[i][0] < 0.0)
            dY[i][0] = 0.0;
    }

    double dmA_sum = 0.0;
    for (int j = 0 ; j < layers_; j++)
    {
        // mass flow rate for reactant gas species
        // dmA is a positive value
        dmA_[j] =   dY[i][j] * (1.0 / (Runiv * T_[i])) * molMass_A_
                   * (MY_4PI * ((radius_[i] * radius_[i]) / (cg_ * cg_)))
                   * TimeStep * nevery;
        dmA_sum += dmA_[j];
        // fix property added so values are outputted to file
        dmA_f_[i][j] = dmA_[j];
    }

    // limit mass change - can't remove more than present in cell
    // limit it with species mass per volume x voidfraction x cell volume / particles in cell x relaxation factor
    // if reacted mass is very small, do not limit it because of numerical stability (division by dmA_sum below)
    if (limit_reactant_consumption_ && dmA_sum > 1e-12)
    {
        if (screenflag_ && screen) fprintf(screen,"checking reactant limitation\n");

        double dAmax = p_A / (Runiv * T_[i]) * molMass_A_ * reactantPerParticle_[i] * maxReactantConsumptionFrac_;

        if (dmA_sum > dAmax)
        {
            for (int j = 0 ; j < layers_; j++)
            {
                dmA_[j] = dAmax * dmA_[j] / dmA_sum;
                dmA_f_[i][j] = dmA_[j];
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, const double *dmA_,const double *v_reac, const double* v_prod)
{
    if (screenflag_ && screen)
        fprintf(screen,"run update atom props \n");

    // based on material change: update relative radii, average density and mass of atom i

    // initialize radius, mass change of layer and sum of particle
    double rad[MAX_LAYERS+1] = {0.};
    double dmL_[MAX_LAYERS+1] = {0.};     // mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass_p_new = 0.0;
    double Cp = 0.0;
    double layer_Cp[MAX_LAYERS+1] = {0.};
    double Tpart = Tpart_[i];
    double H = 0.0;
    double layer_H[4] = {0.};

    if (variableCp_)
    {
        layer_Cp[0] = spec_heat(a_coeff_nasa_Fe,Tpart);
        layer_Cp[1] = spec_heat(a_coeff_nasa_FeO,Tpart);
        layer_Cp[2] = spec_heat(a_coeff_nasa_Fe3O4,Tpart);
        layer_Cp[3] = spec_heat(a_coeff_nasa_Fe2O3,Tpart);
    }
    
    layer_H[0] = conv_enthalpy(a_coeff_nasa_Fe, i);
    layer_H[1] = conv_enthalpy(a_coeff_nasa_FeO,i);
    layer_H[2] = conv_enthalpy(a_coeff_nasa_Fe3O4,i);
    layer_H[3] = conv_enthalpy(a_coeff_nasa_Fe2O3,i);

    // Mass Change of Layers
    // dmL is a positive value, therefore it will be subtracted from the total mass
    // Fe2O3 (iniital)
    dmL_[layers_] = dmA_[layers_-1] * v_reac[layers_-1] * (layerMolMasses_[layers_] / molMass_A_);

    // Initial FeO & Fe3O4 layers
    // changes with active layers
    for (int layer = 1; layer < layers_; layer++)
        dmL_[layer] = -dmA_[layer]   * v_prod[layer]   * (layerMolMasses_[layer] / molMass_A_)
                     + dmA_[layer-1] * v_reac[layer-1] * (layerMolMasses_[layer] / molMass_A_);

    // Fe
    dmL_[0] = -dmA_[0] * v_prod[0] * (layerMolMasses_[0] / molMass_A_);

    // slow decay of 4FeO -> Fe + Fe3O4 at low temperatures could be incorporated at this point

    // New layer masses
    for (int j = 0; j <= layers_; j++)
    {
        massLayer_[i][j] -= dmL_[j]*scale_reduction_rate;
        if (massLayer_[i][j] < 0.0) massLayer_[i][j] = 0.0;
        Cp += massLayer_[i][j] * layer_Cp[j] / layerMolMasses_[j];
        H += massLayer_[i][j] * layer_H[j] / layerMolMasses_[j];
    }
    for (int j = 0; j <= layers_; j++)
    {
        // calculate total mass of particle
        // since there is a minimum radius for layers, there is always a
        // non-zero contribution (at least from the innermost layer)
        sum_mass_p_new += massLayer_[i][j];
    }
    if (variableCp_)
    {
        Cp /= sum_mass_p_new;
        fix_thermal_capacity_->vector_atom[i] = Cp;
    }

    // if (screen) fprintf(screen,"total mass of particle = %f \n", sum_mass_p_new);

    // Total mass of particle with coarse-graining
    pmass_[i] = sum_mass_p_new*cg_*cg_*cg_;

    fix_internal_energy_->vector_atom[i] = H*cg_*cg_*cg_;

    // if (screen) fprintf(screen, "pmass = %f \n",pmass_[i]);

    // Core layer radius (initial Fe2O3)
    rad[layers_] = cbrt((0.75*massLayer_[i][layers_])/(rhoeff_[i][layers_]*M_PI));

    // Outer layer radii (Fe3O4, FeO)
    for (int layer = layers_ - 1; layer > 0; layer--)
        rad[layer]   =   cbrt((0.75*massLayer_[i][layer]/(rhoeff_[i][layer]*M_PI))+rad[layer+1]*rad[layer+1]*rad[layer+1]);

    if (fix_polydisp_)
    {
        for (int layer = layers_; layer > 0; layer--)
        {
            rad[layer] /= cbrt(effvolfactors_[i]);
        }
    }

    // Iron Layer / Ore radius = constant.
    // NOTE: keeping the particle radius constant may introduce a mismatch between
    //       the mass calculated from the chemical reaction and the mass calculated
    //       from the (relative) radii (as in calcMassLayer())!
    rad[0] = radius_[i]/cg_;

    // Determine new relative radii after reduction
    for (int j = 1; j <= layers_; j++)
    {
        relRadii_[i][j] = rad[j]/rad[0];
    }

    relRadii_[i][1] = std::min(0.9999, relRadii_[i][1]);
    relRadii_[i][2] = std::min(0.9998, relRadii_[i][2]);
    relRadii_[i][3] = std::min(0.9997, relRadii_[i][3]);

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

void FixChemShrinkCore::update_gas_properties(int i, const double *dmA_)
{
    /*double kch2_ = 0.0;
    kch2_ = xA_[i] + xC_[i];*/
    double dmA = 0.0;
    // based on material change: update gas-phase source terms for mass and heat
#ifdef PSEUDO_THREE_LAYERS
    for (int j = (T_[i]<SWITCH_LOW_HIGH_TEMPERATURE)?1:0; j < MAX_LAYERS; j++)
#else
    for (int j = 0; j < MAX_LAYERS; j++)
#endif
    {
        dmA = dmA_[j]*cg_*cg_*cg_;
        // Reactant gas mass change
        changeOfA_[i]   -=  dmA;
        // Limit maximum reactant gas
        if (changeOfA_[i] > 0.0) changeOfA_[i] = 0.0;
        // Product gas mass change
        changeOfC_[i]   +=  dmA*molMass_C_/molMass_A_;
    }

    // Limit product gas to the total amount of carbon or hydrogen content
    //changeOfC_[i]   = std::max(kch2_, changeOfC_[i]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::FractionalReduction(int i)
{
    // calculate the fractional reduction as defined in Tang et al. (2012)
    // "Simulation study on performance of Z-path Moving-fluidized Bed for Gaseous Reduction
    // of Iron Ore Fines" ISIJ International, Vol. 52, No. 7, pp. 1241 - 1249

    /*
    const double f_WF = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    const double f_MW = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    const double f_HM = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    */

    // this formulation assumes that molar concentrations n_i of the layers satisfy n_i / nu_i = n_j / nu_j, e.g. n_Hem / 3 = n_Mag / 2
    double f_WF = 1.0 - relRadii_[i][1]*relRadii_[i][1]*relRadii_[i][1];
    double f_MW = 1.0 - relRadii_[i][2]*relRadii_[i][2]*relRadii_[i][2];
    double f_HM = 1.0 - relRadii_[i][3]*relRadii_[i][3]*relRadii_[i][3];

    if (f_HM > 1.0 - SMALL) f_HM = 1.0 - SMALL;
    if (f_MW > 1.0 - SMALL) f_MW = 1.0 - SMALL;
    if (f_WF > 1.0 - SMALL) f_WF = 1.0 - SMALL;

    fracRed_[i][0] = f_WF;
    fracRed_[i][1] = f_MW;
    fracRed_[i][2] = f_HM;

}

/* ---------------------------------------------------------------------- */

/* Heat of Reaction Calcualtion Depending on JANAF thermochemical tables */
void FixChemShrinkCore::heat_of_reaction(int i, const double *dmA_, const double *v_reac, const double *v_prod)
{
    double HR[MAX_LAYERS] = {0.};
    // reaction enthalpy
    double delta_h[MAX_LAYERS] = {0.};
    // conventional enthalpy
    double conv_h[6] = {0.};
    double conv_h54 = 0.0;

    {
#ifdef TWO_LAYERS
        conv_h[0] = conv_enthalpy(a_coeff_nasa_Fe,i);
        conv_h[1] = conv_enthalpy(a_coeff_nasa_Fe3O4,i);
        conv_h[2] = conv_enthalpy(a_coeff_nasa_Fe2O3,i);
#else
        conv_h[0] = conv_enthalpy(a_coeff_nasa_Fe, i);
        conv_h[1] = conv_enthalpy(a_coeff_nasa_FeO,i);
        conv_h[2] = conv_enthalpy(a_coeff_nasa_Fe3O4,i);
        conv_h[3] = conv_enthalpy(a_coeff_nasa_Fe2O3,i);
#endif
    }

    if (strcmp(speciesA, "CO") == 0)
    {
        conv_h[4] = conv_enthalpy(a_coeff_nasa_CO,i);
        conv_h[5] = conv_enthalpy(a_coeff_nasa_CO2,i);
    }
    else if (strcmp(speciesA,"H2")==0)
    {
        conv_h[4] = conv_enthalpy(a_coeff_nasa_H2,i);
        conv_h[5] = conv_enthalpy(a_coeff_nasa_H2O,i);
    }
    conv_h54 = conv_h[5] - conv_h[4];

    // enthalpy changes due to iron oxides
    for (int j = 0; j < layers_; j++)
    {
        delta_h[j] = v_prod[j] * conv_h[j] - v_reac[j] * conv_h[j+1];
    }
    // enthalpy changes due to reduction agent
    for (int j = 0; j < layers_; j++)
    {
#ifdef PSEUDO_THREE_LAYERS
        if (T_[i]<SWITCH_LOW_HIGH_TEMPERATURE && j == 1) continue;
#endif
        delta_h[j] += conv_h54;
    }

    if (screenflag_ && screen) {
        fprintf(screen, "delta_h w %s for reaction 1 is %f \n", speciesA, delta_h[0]);
        fprintf(screen, "delta_h w %s for reaction 2 is %f \n", speciesA, delta_h[1]);
#ifndef TWO_LAYERS
        fprintf(screen, "delta_h w %s for reaction 3 is %f \n", speciesA, delta_h[2]);
#endif
    }

    for (int k = 0; k < layers_; k++)
    {
        HR[k] = delta_h[k]*dmA_[k]/molMass_A_*cg_*cg_*cg_;
    }

    if (screenflag_ && screen) {
        fprintf(screen, "heatFlux of reaction w %s for reaction 1 is %f \n", speciesA, HR[0]);
        fprintf(screen, "heatFlux of reaction w %s for reaction 2 is %f \n", speciesA, HR[1]);
#ifndef TWO_LAYERS
        fprintf(screen, "heatFlux of reaction w %s for reaction 3 is %f \n", speciesA, HR[2]);
#endif
    }

    // add per-particle reactionHeat flux
    for (int k = 0; k < layers_; k++)
        reactionHeat_[i] += HR[k];
}

/* ---------------------------------------------------------------------- */

/* Calculate conventional enthalpies of species */

double FixChemShrinkCore::conv_enthalpy (const double *a, int i)
{
    double value = 0.;

    if (T_[i] < a[0]) { // Temperature smaller than lower bound
        const double Tbound_low = a[0];
        const double Tbound_low_sq = Tbound_low*Tbound_low;
        const double Tbound_low_cb = Tbound_low_sq*Tbound_low;
        value =   a[10]*Tbound_low
                + a[11]*Tbound_low_sq*0.5
                + a[12]*Tbound_low_cb/3.0
                + a[13]*Tbound_low_sq*Tbound_low_sq*0.25
                + a[14]*Tbound_low_sq*Tbound_low_cb*0.20
                + a[15];
    } else if (T_[i] < a[2]) {
        const double Ti = T_[i];
        const double Ti_sq = Ti*Ti;
        const double Ti_cb = Ti_sq*Ti;
        value =   a[10]*Ti
                + a[11]*Ti_sq*0.5
                + a[12]*Ti_cb/3.0
                + a[13]*Ti_sq*Ti_sq*0.25
                + a[14]*Ti_sq*Ti_cb*0.20
                + a[15];
    } else if (T_[i] < a[1]) {
        const double Ti = T_[i];
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

double FixChemShrinkCore::spec_heat (const double *a, double Ti)
{
    double value = 0.;

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

// Equilibrium constant for low temperature reactions

double FixChemShrinkCore::K_eq_low(int layer, int i)
{
    double Keq_low = 0.0;
#ifdef TWO_LAYERS
    // 0 = Fe3O4 (magnetite) -> Fe (iron)
    // 1 = Fe2O3 (hematite)  -> Fe3O4 (magnetite)
    if (strcmp(speciesA, "CO") == 0)
    {
        if (layer == 1) {
            Keq_low = exp(3968.37/T_[i]+3.94);
        } else if (layer == 0) {
            Keq_low = 0.97949;// =pow(10.0,-0.009);
        }
    }
    else if(strcmp(speciesA,"H2")==0)
    {
        if (layer == 1)
            Keq_low = exp(-362.6/T_[i] + 10.334);
        else if (layer == 0)
            Keq_low = pow(10.0,(-1742.0/T_[i]+1.568));
    }
#else
    if (strcmp(speciesA, "CO") == 0)
    {
        switch(layer) {
        case LAYER_HEMATITE:
            Keq_low = exp(3968.37/T_[i]+3.94); break;
        case LAYER_MAGNETITE:
            Keq_low = 0.97949; break; // = pow(10.0,-0.009);
        case LAYER_WUSTITE:
            Keq_low = 0.97949; break; // = pow(10.0,-0.009);
        default:
            error->fix_error(FLERR, this, "invalid layer in equilibrium constant method");
            break;
        }
    }
    else if (strcmp(speciesA,"H2")==0)
    {
        switch(layer) {
        case LAYER_HEMATITE:
            Keq_low = exp(-362.6/T_[i] + 10.334); break;
        case LAYER_MAGNETITE:
            Keq_low = pow(10.0,(-1742.0/T_[i]+1.568)); break;
        case LAYER_WUSTITE:
            Keq_low = pow(10.0,(-1742.0/T_[i]+1.568)); break;
        default:
            error->fix_error(FLERR, this, "invalid layer in equilibrium constant method");
            break;
        }
    }
#endif
    else
    {
        error->fix_error(FLERR, this, "Undefined Reaction \n");
    }

    if (screenflag_ && screen)
        fprintf(screen,"Keq_low : %f \n",Keq_low);

    return Keq_low;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction_low(int i, double *dmA_, const double *x0_eq_)
{
    double p_eq_[MAX_LAYERS] = {0.};

    for (int layer = 0; layer < layers_; layer++)
    {
        p_eq_[layer] = x0_eq_[layer] * partP_[i];
    }

    const double p_A = xA_[i] * partP_[i];

    if (screenflag_ && screen)
    {
        fprintf(screen, "p_eq_I: %f, p_eq_II: %f, p_A: %f \n", p_eq_[0], p_eq_[1], p_A);
    }

#ifdef TWO_LAYERS
    if (layers_ == 2)
    {
        const double A0 = Aterm[i][0];
        const double A1plusB1         = Aterm[i][1] + Bterm[i][1];
        const double B0plusMass       = Bterm[i][0] + Massterm[i];
        const double A0plusB0plusMass = A0 + B0plusMass;

        const double W = A1plusB1 * A0plusB0plusMass + A0 * B0plusMass;

        // hematite to magnetite
        dY[i][1]   =   (A0plusB0plusMass * (p_A - p_eq_[1]) - B0plusMass * (p_A - p_eq_[0])) / W;

        if (dY[i][1] < 0.0)
            dY[i][1] = 0.0;

        // magnetite to iron
        if (dY[i][1] == 0.0)
            dY[i][0] = 0.0;
        else
            dY[i][0] = ((A1plusB1 + B0plusMass) * (p_A - p_eq_[0]) - B0plusMass * (p_A - p_eq_[1])) / W;

        if (dY[i][0] < 0.0)
            dY[i][0] = 0.0;

    }
    else if (layers_ == 1)
    {
        // rate of chemical reaction for 1 active layer
        const double W = Aterm[i][0] + Bterm[i][0] + Massterm[i];

        // hematite to magnetite
        dY[i][1] = 0.0;

        // magnetite to iron
        dY[i][0]   =   (p_A - p_eq_[0]) / W;

        if (dY[i][0] < 0.0)
            dY[i][0] = 0.0;
    }
#else
#ifdef PSEUDO_THREE_LAYERS
    if (layers_ == 3)
    {
        const double A0 = Aterm[i][1];
        const double A1plusB1         = Aterm[i][2] + Bterm[i][2];
        const double B0plusMass       = Bterm[i][1] + Massterm[i];
        const double A0plusB0plusMass = A0 + B0plusMass;

        const double W = A1plusB1 * A0plusB0plusMass + A0 * B0plusMass;

        // hematite to magnetite
        dY[i][2]   =   (A0plusB0plusMass * (p_A - p_eq_[2]) - B0plusMass * (p_A - p_eq_[1])) / W;

        if (dY[i][2] < 0.0)
            dY[i][2] = 0.0;

        // magnetite to iron
        if (dY[i][2] == 0.0)
            dY[i][1] = 0.0;
        else
            dY[i][1] = ((A1plusB1 + B0plusMass) * (p_A - p_eq_[1]) - B0plusMass * (p_A - p_eq_[2])) / W;

        if (dY[i][1] < 0.0)
            dY[i][1] = 0.0;
    }
    else if (layers_ == 2)
    {
        // rate of chemical reaction for 1 active layer
        const double W = Aterm[i][1] + Bterm[i][1] + Massterm[i];

        // hematite to magnetite
        dY[i][2] = 0.0;

        //magnetite to iron
        dY[i][1]   =   (p_A - p_eq_[1]) / W;

        if (dY[i][1] < 0.0)
            dY[i][1] = 0.0;
    }
    else if (layers_ == 1)
    {
        // should never get here (see post_force)
        error->fix_error(FLERR, this, "Trying to reduce wuestite layer at low T");
    }
#else

#endif
    dY[i][0] = dY[i][1];
#endif

    double dmA_sum = 0.0;
    for (int j = 0 ; j < layers_; j++)
    {
        // mass flow rate for reactant gas species
        // dmA is a positive value
        dmA_[j] =   dY[i][j]*(1.0/(Runiv*T_[i]))*molMass_A_*(MY_4PI*((radius_[i]*radius_[i])/(cg_*cg_)))*TimeStep*nevery;
        dmA_sum += dmA_[j];
        // fix property added so values are outputted to file
        dmA_f_[i][j] = dmA_[j];
    }

    // limit mass change - can't remove more than present in cell
    // limit it with species mass per volume x voidfraction x cell volume / particles in cell x relaxation factor
    // if reacted mass is very small, do not limit it because of numerical stability (division by dmA_sum below)
    if (limit_reactant_consumption_ && dmA_sum > 1e-12)
    {
        if (screenflag_ && screen) fprintf(screen,"checking reactant limitation\n");

        double dAmax = p_A / (Runiv * T_[i]) * molMass_A_ * reactantPerParticle_[i] * maxReactantConsumptionFrac_;

        if (dmA_sum > dAmax)
        {
            for (int j = 0 ; j < layers_; j++)
            {
                dmA_[j] = dAmax * dmA_[j] / dmA_sum;
                dmA_f_[i][j] = dmA_[j];
            }
        }
    }

#ifndef TWO_LAYERS
#ifdef PSEUDO_THREE_LAYERS
    dmA_f_[i][0] = dmA_[0] = dmA_[1];
#else

#endif
#endif
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::FractionalReduction_low(int i)
{
    double f_WF = 1.0 - relRadii_[i][1]*relRadii_[i][1]*relRadii_[i][1];
    double f_MW = 1.0 - relRadii_[i][2]*relRadii_[i][2]*relRadii_[i][2];

    if (f_MW > 1.0 - SMALL) f_MW = 1.0 - SMALL;
    if (f_WF > 1.0 - SMALL) f_WF = 1.0 - SMALL;

    fracRed_[i][0] = f_WF;
    fracRed_[i][1] = f_MW;
    fracRed_[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getXi_low(int i, double *x0_eq_)
{
    const double kch2_ = xA_[i] + xC_[i];

    for (int j = 0; j < layers_; j++)
    {
        x0_eq_[j]  =   kch2_/(1.0+K_eq_low(j,i));
    }

#ifdef TWO_LAYERS
    x0_eq_[2] = 0.;
#endif

    if (screenflag_ && screen)
        fprintf(screen, "x0_eq_0: %f, x0_eq_1: %f, x0_eq_2: %f \n", x0_eq_[0], x0_eq_[1], x0_eq_[2]);
}

/* ---------------------------------------------------------------------- */

// 0 = magnetite interface, 1 = hematite interface
void FixChemShrinkCore::getA_low(int i)
{
    if (strcmp(speciesA, "CO") == 0)
    {
        for (int j = 0; j < layers_ ; j++)
        {
            Aterm[i][j] = (k0_low_CO[j] * exp(-Ea_low_CO[j] / (Runiv * T_[i])))
                        * cbrt(square(1.0 - fracRed_[i][j]))
                        * (1.0 + 1.0 / K_eq_low(j,i));
            Aterm[i][j] = 1.0 / Aterm[i][j];
        }
    }
    else if(strcmp(speciesA,"H2") == 0)
    {
        for (int j = 0; j < layers_ ; j++)
        {
            Aterm[i][j] = (k0_low_H2[j] * exp(-Ea_low_H2[j] / (Runiv * T_[i])))
                        * cbrt(square(1.0 - fracRed_[i][j]))
                        * (1.0 + 1.0 / K_eq_low(j,i));
            Aterm[i][j] = 1.0 / Aterm[i][j];
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::init_defaults()
{
    molMass_A_ = molMass_C_ = 0.0;
    rhoeff_ = NULL;
    porosity_ = pore_diameter_ = NULL;
    tortuosity_ = 0.0;
    relRadii_ = massLayer_ = NULL;
    k0_ = Ea_ = NULL;
    xA_ = xC_ = NULL;
    scale_reduction_rate = 1.;
    layerDensities_ = NULL;

    // particle properties total
    radius_ = pmass_ = pdensity_ = NULL;

    // initialise fix handles
    changeOfA_ = changeOfC_ = T_ = Tpart_ = molecularDiffusion_ = nuf_ = Rep_ = partP_ = Massterm = reactionHeat_ = NULL;
    Aterm = Bterm = effDiffBinary = effDiffKnud = fracRed_ = NULL;

    dY = dmA_f_ = NULL;

    TimeStep = 0.0;

    fix_changeOfA_ = NULL;
    fix_changeOfC_ = NULL;
    fix_tgas_ = NULL;       // [K]
    fix_tpart_ = NULL;       // [K]
    fix_reactionHeat_= NULL;
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
    fix_layerDens_ = NULL;      //  [kg/m^3]
    fix_k0_ = NULL;             //  [m/s]
    fix_Ea_ = NULL;             //  [J/mol] - [kg*m^2/s^2*mol]
    fix_porosity_ = NULL;       //  [%]
    fix_rhoeff_ = NULL;
    fix_thermal_capacity_ = NULL;
    fix_internal_energy_ = NULL;
    fix_tortuosity_ = NULL;
    fix_pore_diameter_ = NULL;  // [m]
    fix_dY_ = NULL;
    fix_dmA_ = NULL;

    fix_reactantPerParticle_ = NULL;
    reactantPerParticle_ = NULL;
    limit_reactant_consumption_ = false;
    maxReactantConsumptionFrac_ = 0.5;

    massA = massC = NULL;
    diffA = moleFracA = moleFracC = NULL;
    speciesA = speciesC = NULL;

    variableCp_ = false;
}

void FixChemShrinkCore::update_fix(int narg, char **arg)
{
    setup(0);
    int nlocal = atom->nlocal;
    int *mask  = atom->mask;

    for (int i = 0; i < nlocal; ++i)
    {
        if (mask[i] & groupbit)
        {
            active_layers(i);
            double m = 0.0;
            double Cp = 0.0;
            double layer_Cp[MAX_LAYERS+1] = {0.};
            double Tpart = Tpart_[i];

            for (int layer = 0 ; layer <= layers_; layer++)
            {
                m += massLayer_[i][layer];
            }

            if (variableCp_)
            {
                layer_Cp[0] = spec_heat(a_coeff_nasa_Fe,Tpart);
                layer_Cp[1] = spec_heat(a_coeff_nasa_FeO,Tpart);
                layer_Cp[2] = spec_heat(a_coeff_nasa_Fe3O4,Tpart);
                layer_Cp[3] = spec_heat(a_coeff_nasa_Fe2O3,Tpart);
                for (int layer = 0 ; layer <= layers_; layer++)
                {
                    Cp += massLayer_[i][layer] * layer_Cp[layer] / layerMolMasses_[layer];
                }
                Cp /= m;
                fix_thermal_capacity_->vector_atom[i] = Cp;
            }

            pmass_[i] = m*cg_*cg_*cg_;
            pdensity_[i] = 0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);
            if (fix_polydisp_)
            {
                pdensity_[i] /= effvolfactors_[i];
            }
        }
    }
}
