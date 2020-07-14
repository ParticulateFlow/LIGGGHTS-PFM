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



/* ---------------------------------------------------------------------- */

FixChemShrinkCoreSingle::FixChemShrinkCoreSingle(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    layers_(1),
    minMolarFrac_(1e-3),
    rmin_(1e-5),      //  [m]
    k0_(-1.0),
    T0_(-1.0),
    Tmin_(0.0),
    nPreFactor_(0),
    created_fix_porosity_(false)
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
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"Tmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            Tmin_ = atof(arg[iarg_+1]);
            if (Tmin_ <= 0.)
                error -> fix_error(FLERR, this, "Tmin is not (well-)defined");
            iarg_ +=2;
        }
        else if (strcmp(arg[iarg_],"nPreFactor") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            nPreFactor_ = atoi(arg[iarg_+1]);
            if (nPreFactor_ < 0)
                error -> fix_error(FLERR, this, "nPreFactor is not (well-)defined");
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
        else
        {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    if (k0_ < 0.0 || T0_ < 0.0)
    {
        error->fix_error(FLERR,this,"specify positive values for k0 and T0");
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
        if (fix_porosity_ && created_fix_porosity_) { modify->delete_fix(fix_porosity_->id); porosity_ = NULL; }
#ifdef NO_CFD_COUPLING
        if (fix_changeOfA_)     { modify->delete_fix(fix_changeOfA_->id); changeOfA_ = NULL; }
        if (fix_changeOfC_)     { modify->delete_fix(fix_changeOfC_->id); changeOfC_ = NULL; }
        if (fix_reactionHeat_)  { modify->delete_fix(fix_reactionHeat_->id); reactionHeat_ = NULL; }
#endif
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
    reactionHeat_   =   fix_reactionHeat_   ->  vector_atom;
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
    //
    tortuosity_     =   fix_tortuosity_     ->  compute_scalar();
    pore_diameter_  =   fix_pore_diameter_  ->  compute_scalar();
    //
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
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, style));
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
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::setup(int)
{
    updatePtrs();

    // set initial values for rhoeff, and use them to calculate mass of layers
    for (int i = 0; i < atom->nlocal; ++i)
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
    double Keq_ = 1e6;
    return Keq_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::getXi(int i, double x0_eq_)
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
// Use binary diffusion for mixture, and knudsen diffusion to determine the effective diffusion term
void FixChemShrinkCoreSingle::getB(int i)
{
    double fracRedThird_ = 0.;
    double diffEff_ = 0.;


    // calculate fractional reduction to the power of 1/3 for simpler use
    fracRedThird_ = cbrt(1.0-fracRed_[i][0]);

    // Calculate the effective molecular diffusion
    effDiffBinary[i] = molecularDiffusion_[i]*(porosity_[0]/tortuosity_) + SMALL;

    // Calculate the knudsen diffusion
    effDiffKnud[i]  =  (pore_diameter_/6.0)*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_))*(porosity_[0]/tortuosity_) + SMALL;

    // total effective diffusivity
    diffEff_ =   effDiffKnud[i]*effDiffBinary[i]/(effDiffBinary[i]+effDiffKnud[i]) + SMALL;


    // calculation of diffusion term
    if (molecularDiffusion_[i] < SMALL)
    {
        Bterm[i]= 0.0;
    }
    else
    {
        // diffusion resistance through Fe
        Bterm[i]   =   ((1-fracRedThird_)/fracRedThird_)*((radius_[i]/cg_)/(diffEff_));
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

void FixChemShrinkCoreSingle::reaction(int i, double dmA_, const double x0_eq_)
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
    // fix property added so values are otuputted to file
    dmA_f_[i] = dmA_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::update_atom_properties(int i, const double dmA_)
{
    if (screenflag_ && screen)
        fprintf(screen,"run update atom props \n");

    // based on material change: update relative radii, average density and mass of atom i

    // initialize radius, mass change of layer and sum of particle
    double rad[2] = {0.};
    double dmL_[2] = {0.};     // mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass_p_new = 0.0;

    // Mass change of inner layer
    dmL_[1] = dmA_ * nu_B_ * (molMass_B_ / molMass_A_);

    // New layer mass
    massLayer_[i][1] -= dmL_[1]*scale_reduction_rate;
    if (massLayer_[i][1] < 0.0)
        massLayer_[i][1] = 0.0;

    // NOTE: This keeps the particle radius and layer porosities and densities constant,
    //       hence there is a slight mass conservation error.
    //       Alternatively, one could provide details on the chemical reaction for the product layer.

    dmL_[0] = dmL_[1] * rhoeff_[i][0] * porosity_[0] / (rhoeff_[i][1] * porosity_[1]);
    massLayer_[i][0] += dmL_[0]*scale_reduction_rate;

///
    for (int j = 0; j <= 1; j++)
    {
        // calculate total mass of particle
        // since there is a minimum radius for layers, there is always a
        // non-zero contribution (at least from the innermost layer)
        sum_mass_p_new += massLayer_[i][j];
    }

    // if (screen) fprintf(screen,"total mass of particle = %f \n", sum_mass_p_new);

    // Total mass of particle with coarse-graining
    pmass_[i] = sum_mass_p_new*cg_*cg_*cg_;

    // if (screen) fprintf(screen, "pmass = %f \n",pmass_[i]);
    rad[1] = cbrt((0.75*massLayer_[i][1])/(rhoeff_[i][1]*M_PI));
    rad[0] = radius_[i]/cg_;

    if (fix_polydisp_)
    {
        rad[1] /= cbrt(effvolfactors_[i]);
    }

    // Determine new relative radii after reduction

    relRadii_[i][1] = rad[1]/rad[0];
    relRadii_[i][1] = std::min(0.9999, relRadii_[i][1]);

    // total particle effective density
    // pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);

    if (screenflag_ && screen)
        fprintf(screen, "radius_: %f, pmass_: %f, pdensity_: %f\n ", radius_[0], pmass_[0], pdensity_[0]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::update_gas_properties(int i, const double dmA_)
{
    /*double kch2_ = 0.0;
    kch2_ = xA_[i] + xC_[i];*/

    // based on material change: update gas-phase source terms for mass and heat

    // Reactant gas mass change
    changeOfA_[i]   -=  dmA_;
    // Limit maximum reactant gas
    if (changeOfA_[i] > 0.0) changeOfA_[i] = 0.0;
    // Product gas mass change
    changeOfC_[i]   +=  dmA_*nu_C_*molMass_C_/molMass_A_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::FractionalReduction(int i)
{
    const double f = 1.0 - relRadii_[i][1]*relRadii_[i][1]*relRadii_[i][1];

    fracRed_[i][0] = f;
}

/* ---------------------------------------------------------------------- */

/* Heat of Reaction Calcualtion Depending on JANAF thermochemical tables */
void FixChemShrinkCoreSingle::heat_of_reaction(int i, const double dmA_)
{
 
}


/* ---------------------------------------------------------------------- */

void FixChemShrinkCoreSingle::init_defaults()
{
    molMass_A_ = molMass_C_ = 0.0;
    rhoeff_ = NULL;
    porosity_ = NULL;
    pore_diameter_ = tortuosity_ = 0.0;
    relRadii_ = massLayer_ = NULL;
    xA_ = xC_ = NULL;
    scale_reduction_rate = 1.;
    layerDensities_ = NULL;

    // particle properties total
    radius_ = pmass_ = pdensity_ = NULL;

    // initialise fix handles
    changeOfA_ = changeOfC_ = T_ = molecularDiffusion_ = nuf_ = Rep_ = partP_ = Massterm = reactionHeat_ = NULL;
    Aterm = Bterm = effDiffBinary = effDiffKnud = NULL;
    fracRed_ = NULL;

    dY = dmA_f_ = NULL;

    TimeStep = 0.0;

    fix_changeOfA_ = NULL;
    fix_changeOfC_ = NULL;
    fix_tgas_ = NULL;       // [K]
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
    fix_porosity_ = NULL;       //  [%]
    fix_rhoeff_ = NULL;
    fix_tortuosity_ = NULL;
    fix_pore_diameter_ = NULL;  // [m]
    fix_dY_ = NULL;
    fix_dmA_ = NULL;

    massA = massC = NULL;
    diffA = moleFracA = moleFracC = NULL;
    speciesA = speciesC = NULL;
}

void FixChemShrinkCoreSingle::update_fix(int narg, char **arg)
{
    setup(0);
}
