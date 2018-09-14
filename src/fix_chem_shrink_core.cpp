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
#include "fix_chem_shrink_core.h"
#include "fix_property_atom.h"
#include "pair_gran.h"
#include "compute_pair_gran_local.h"
#include "fix_property_global.h"
#include "properties.h"
#include "property_registry.h"
#include "global_properties.h"
#include "force.h"
#include "group.h"
#include "vector_liggghts.h"
#include "math_const.h"
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define SMALL   1e-10
#define Runiv   8.3144

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    nmaxlayers_(3),
    layers_(nmaxlayers_),
    minMolarFrac_(1e-3),
    rmin_(1e-5)       //  [m]
{
    if ((strncmp(style, "chem/shrink/core", 15) == 0) && ((!atom->radius_flag) || (!atom->rmass_flag)))
        error->all(FLERR, "Fix chem/shrink/core needs per particle radius and mass");

    // set defaults
    init_defaults();
    screenflag_ =   0;
    massA = massC = NULL;
    diffA = moleFracA = moleFracC = NULL;
    speciesA = speciesC = NULL;
    molarConc_  =   NULL;
    dY_previous3 = dY_previous2 = false;
    cg_ = 0.0;
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
        else if (strcmp(arg[iarg_],"screen") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) screenflag_ = 1;
            else if (strcmp(arg[iarg_+1],"no") == 0) screenflag_ = 0;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            iarg_ += 2; // iarg = 13
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nevery") == 0)
        {
            nevery = atoi(arg[iarg_+1]);
            if (nevery <= 0) error->fix_error(FLERR,this,"nevery must be larger than 0");
            iarg_+=2;
            hasargs = true; //iarg = 15
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
    cg_ = force->cg();
    peratom_flag = 1;
    peratom_freq = 1;
    global_freq = 1;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_create()
{
    if (fix_Aterm==NULL) {
        char *fixname = new char [strlen("Aterm_")+strlen(id)+1];
        strcpy (fixname,"Aterm_");
        strcat(fixname,id);

        const char* fixarg[11];
        // register property/atom for chemical reaction resistance
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_Bterm==NULL) {
        char *fixname = new char[strlen("Bterm_")+strlen(id)+1];
        strcpy (fixname,"Bterm_");
        strcat(fixname,id);

        const char* fixarg[11];
        // register property/atom for diffusion resistance term
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_Massterm==NULL) {
        char *fixname = new char[strlen("Massterm_")+strlen(id)+1];
        strcpy (fixname,"Massterm_");
        strcat(fixname,id);

        const char* fixarg[9];
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_effDiffBinary==NULL) {
        char *fixname = new char[strlen("effDiffBinary_")+strlen(id)+1];
        strcpy (fixname,"effDiffBinary_");
        strcat(fixname,id);

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_effDiffKnud==NULL)    {
        char *fixname = new char[strlen("effDiffKnud_")+strlen(id)+1];
        strcpy (fixname,"effDiffKnud_");
        strcat(fixname,id);

        const char* fixarg[11];
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_dY_==NULL)     {
        char* fixname = new char[strlen("dY_")+strlen(id)+1];
        strcpy (fixname,"dY_");
        strcat(fixname,id);

        const char* fixarg[12];
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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

    if (fix_dmA_==NULL)     {
        char* fixname = new char[strlen("dmA_")+strlen(id)+1];
        strcpy (fixname,"dmA_");
        strcat(fixname,id);

        const char* fixarg[12];
        fixarg[0]=fixname;            // fixid
        fixarg[1]="all";
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
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
{
    char *fixname;

    if (unfixflag)
    {
        if (fix_changeOfA_)     { modify  ->  delete_fix(massA); delete [] changeOfA_;}
        if (fix_changeOfC_)     { modify  ->  delete_fix(massC); delete [] changeOfC_;}
        if (fix_tgas_)          { modify  ->  delete_fix("partTemp"); delete [] T_;}
        if (fix_reactionHeat_)  { modify  ->  delete_fix("reactionHeat"); delete [] reactionHeat_;}
        if (fix_diffcoeff_)     { modify  ->  delete_fix(diffA); delete [] molecularDiffusion_; }
        if (fix_nuField_)       { modify  ->  delete_fix("partNu"); delete [] nuf_; }
        if (fix_partRe_)        { modify  ->  delete_fix("partRe"); delete [] Rep_; }

        if (fix_moleFractionA_) { modify  ->  delete_fix(moleFracA); delete [] xA_; }
        if (fix_moleFractionC_) { modify  ->  delete_fix(moleFracC); delete [] xC_; }
        if (fix_fracRed)        { modify  ->  delete_fix("fracRed"); delete [] fracRed_;}

        fixname = new char[strlen("Aterm_")+strlen(id)+1];
        strcpy (fixname,"Aterm_");
        strcat(fixname,id);
        if (fix_Aterm)          { modify  ->  delete_fix(fixname); delete [] Aterm; }
        delete []fixname;
        fixname = new char[strlen("Bterm_")+strlen(id)+1];
        strcpy (fixname,"Bterm_");
        strcat(fixname,id);
        if (fix_Bterm)          { modify  ->  delete_fix(fixname); delete [] Bterm; }
        delete []fixname;
        fixname = new char[strlen("Massterm_")+strlen(id)+1];
        strcpy (fixname,"Massterm_");
        strcat(fixname,id);
        if (fix_Massterm)       { modify  ->  delete_fix(fixname); delete [] Massterm; }
        delete []fixname;
        fixname = new char[strlen("effDiffBinary_")+strlen(id)+1];
        strcpy (fixname,"effDiffBinary_");
        strcat(fixname,id);
        if (fix_effDiffBinary)  { modify  ->  delete_fix(fixname); delete [] effDiffBinary;}
        delete []fixname;
        fixname = new char[strlen("effDiffKnud_")+strlen(id)+1];
        strcpy (fixname,"effDiffKnud_");
        strcat(fixname,id);
        if (fix_effDiffKnud)    { modify  ->  delete_fix(fixname); delete [] effDiffKnud;}
        delete [] fixname;

        if (fix_partPressure_)  { modify  ->  delete_fix("partP"); delete [] partP_;}
        if (fix_layerRelRad_)   { modify  ->  delete_fix("relRadii"); delete [] relRadii_; }
        if (fix_layerMass_)     { modify  ->  delete_fix("massLayer"); delete [] massLayer_; }
        if (fix_porosity_)      { modify  ->  delete_fix("porosity_"); delete [] porosity_;}
        if (fix_rhoeff_)        { modify  ->  delete_fix("rhoeff"); delete [] rhoeff_;}

        if(fix_totalMole_)      { modify  ->  delete_fix("partMolarConc"); delete [] molarConc_;}

        fixname = new char[strlen("dY_")+strlen(id)+1];
        strcpy (fixname,"dY_");
        strcat(fixname,id);
        if(fix_dY_)      { modify  ->  delete_fix(fixname); delete [] dY;}
        delete []fixname;

        fixname = new char[strlen("dmA_")+strlen(id)+1];
        strcpy (fixname,"dmA_");
        strcat(fixname,id);
        if(fix_dmA_)      { modify  ->  delete_fix(fixname); delete [] dmA_f_;}
        delete []fixname;
    }
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
    delete []massA;
    delete []massC;
    delete []diffA;
    delete []moleFracA;
    delete []moleFracC;

    delete []speciesA;
    delete []speciesC;
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCore::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::updatePtrs()
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
    layerDensities_ =   fix_layerDens_           ->  get_values();
    layerMolMasses_ =   fix_layerMolMass_        ->  get_values();
    k0_             =   fix_k0_             ->  get_values();
    Ea_             =   fix_Ea_             ->  get_values();
    porosity_       =   fix_porosity_       ->  array_atom;
    rhoeff_         =   fix_rhoeff_         ->  array_atom;
    //
    tortuosity_     =   fix_tortuosity_     ->  compute_scalar();
    pore_diameter_  =   fix_pore_diameter_  ->  compute_scalar();
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

    molarConc_  =   fix_totalMole_  ->  vector_atom;
    dY          =   fix_dY_         ->  array_atom;
    dmA_f_    =   fix_dmA_ -> array_atom;
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

    // look up pre-exponential factor k0
    // differs for every ore id
    int ntype = atom -> ntypes;
    char* fixname = new char[strlen("k0_")+strlen(id)+1];
    strcpy (fixname,"k0_");
    strcat(fixname,id);
    fix_k0_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,style));
    delete []fixname;

    // look up activation energies Ea
    // differs for every ore id
    fixname = new char [strlen("Ea_")+strlen(id)+1];
    strcpy(fixname, "Ea_");
    strcat(fixname, id);
    fix_Ea_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname, "property/global", "vector", ntype, 0, style));
    delete[]fixname;

    // Layer Molar Mass
    fixname = new char [strlen("molMass_")+strlen(group->names[igroup])+1];
    strcpy(fixname,"molMass_");
    strcat(fixname,group->names[igroup]);
    fix_layerMolMass_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,style));
    delete []fixname;

    // Layer Density
    fixname = new char [strlen("density_")+strlen(group->names[igroup])+1];
    strcpy(fixname,"density_");
    strcat(fixname,group->names[igroup]);
    fix_layerDens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,style));
    delete []fixname;

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
    fix_fracRed         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("fracRed", "property/atom", "vector", 0, 0, style));
    fix_partPressure_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partP", "property/atom", "scalar", 0, 0, style));
    fix_layerRelRad_    =   static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii", "property/atom", "vector", 0, 0, style));
    fix_layerMass_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer","property/atom","vector",0,0,style));
    fix_porosity_       =   static_cast<FixPropertyAtom*>(modify->find_fix_property("porosity_", "property/atom", "vector", 0, 0, style));
    fix_rhoeff_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("rhoeff", "property/atom", "vector", 0, 0, style));
    // references for global properties - valid for every particle equally
    fix_tortuosity_     =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("tortuosity_", "property/global", "scalar", 0, 0, style));
    fix_pore_diameter_  =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("pore_diameter_", "property/global", "scalar", 0, 0,style));
    fix_totalMole_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partMolarConc","property/atom","scalar",0,0,style));

    fixname = new char [strlen("Aterm_")+strlen(id)+1];
    strcpy (fixname,"Aterm_");
    strcat(fixname,id);
    fix_Aterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 0, 0, style));
    delete []fixname;


    fixname = new char[strlen("Bterm_")+strlen(id)+1];
    strcpy (fixname,"Bterm_");
    strcat(fixname,id);
    fix_Bterm           =   static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 3, 0, style));
    delete [] fixname;

    fixname = new char[strlen("Massterm_")+strlen(id)+1];
    strcpy (fixname,"Massterm_");
    strcat(fixname,id);
    fix_Massterm        =   static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "scalar", 0, 0, style));
    delete [] fixname;

    fixname = new char[strlen("effDiffBinary_")+strlen(id)+1];
    strcpy (fixname,"effDiffBinary_");
    strcat(fixname,id);
    fix_effDiffBinary = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 3, 0, style));
    delete [] fixname;

    fixname = new char[strlen("effDiffKnud_")+strlen(id)+1];
    strcpy (fixname,"effDiffKnud_");
    strcat(fixname,id);
    fix_effDiffKnud = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 3, 0, style));
    delete [] fixname;

    fixname = new char[strlen("dY_")+strlen(id)+1];
    strcpy (fixname,"dY_");
    strcat(fixname,id);
    fix_dY_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 0, 0, style));
    delete [] fixname;

    fixname = new char[strlen("dmA_")+strlen(id)+1];
    strcpy (fixname,"dmA_");
    strcat(fixname,id);
    fix_dmA_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixname, "property/atom", "vector", 0, 0, style));
    delete [] fixname;

    updatePtrs();

    // get initial values for rhoeff, and use them to calculate mass of layers
    for (int i = 0; i < atom->nlocal; ++i)
    {
        rhoeff_[i][layers_] = pdensity_[i];
        for (int layer=0; layer < layers_; layer++)
        {
            rhoeff_[i][layer] = (1.0 - porosity_[i][layer])*layerDensities_[layer];
        }
        calcMassLayer(i);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    updatePtrs();
    int i = 0;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    double x0_eq_[nmaxlayers_] = {0.};          // molar fraction of reactant gas
    double dmA_[nmaxlayers_] = {0.};            // mass flow rate of reactant gas species for each layer at w->fe, m->w & h->m interfaces

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
                active_layers(i);
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
                }
        }
    }
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCore::active_layers(int i)
{
    for(int j  = layers_; j > 0; j--)
    {
        if (relRadii_[i][j]*(radius_[i]/cg_) < rmin_)
        {
            layers_ -= 1;
            calcMassLayer(i);
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
    double rad[nmaxlayers_+1] = {0.};
    for (int layer = 0; layer <= layers_ ; layer++)
        rad[layer] = (radius_[i]/cg_)*relRadii_[i][layer];

    massLayer_[i][layers_]   =   MY_4PI3*rad[layers_]*rad[layers_]*rad[layers_]*rhoeff_[i][layers_];
    for (int layer = 0 ; layer < layers_; layer++)
    {
        massLayer_[i][layer]   =   MY_4PI3*(rad[layer]*rad[layer]*rad[layer]-rad[layer+1]*rad[layer+1]*rad[layer+1])*rhoeff_[i][layer];
    }
}

/* ---------------------------------------------------------------------- */
double FixChemShrinkCore::K_eq(int layer, int i)
{
    // 0 = w->fe , 1 = m->w, 2 = h->m;
    double Keq_ = 0.;

    if (strcmp(speciesA, "CO") == 0)
    {
        if (layer == 2)
            Keq_ = exp(3968.37/T_[i]+3.94);
        else if (layer == 1)
            Keq_ = pow(10.0,(-1834.0/T_[i]+2.17));
            //Keq_ = exp(-3585.64/T_[i]+4.58);
        else if (layer == 0)
            Keq_ = pow(10.0,(914.0/T_[i]-1.097));
            //Keq_ = exp(2744.63/T_[i]-2.946);
     }
     else if(strcmp(speciesA,"H2")==0)
     {
        if (layer == 2)
            Keq_   =   exp(-362.6/T_[i] + 10.334);
        else if (layer == 1)
            Keq_    =   pow(10.0,(-3577.0/T_[i]+3.74));
            //Keq_   =   exp(-7916.6/T_[i] + 8.46);
        else if (layer == 0)
            //Keq_    =   pow(10.0,(-827.0/T_[i]+0.468)); // --> With this equilibrium constant R1 is occuring with all temperatures for H2
            Keq_    =   pow(10.0,(-856.66/T_[i]+0.4387)); // Equilibrium constant from Turkdogan "Physical Chemistry of High Temperature Technology"
            //Keq_   =   exp(-1586.9/T_[i] + 0.9317);
     }
     else
     {
         printf("Error : Undefined Reaction \n");
     }

     if(screenflag_ && screen)
         fprintf(screen,"Keq_ : %f \n",Keq_);

    return Keq_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getXi(int i, double *x0_eq_)
{
    double kch2_ = 0.0;
    kch2_ = xA_[i] + xC_[i];

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
    updatePtrs();
    for (int j = 0; j < layers_ ; j++)
    {
        Aterm[i][j]   =   (k0_[j]*exp(-Ea_[j]/(Runiv*T_[i])))*cbrt((1.0-fracRed_[i][j])*(1.0-fracRed_[i][j]))*(1+1/K_eq(j,i));
        Aterm[i][j]   =   1.0/Aterm[i][j];
    }

}

/* ---------------------------------------------------------------------- */

// Calculate B [s/m] - the diffusion resistance term
// Use binary diffusion for mixture, and knudsen diffusion to determine the effective diffusion term
// 0 : diffusion through iron layer, 1 : diffusion through wüstite, 2 : diffusion through magnetite.
// there is no diffusion through the hematite layer
void FixChemShrinkCore::getB(int i)
{
    updatePtrs();
    double fracRedThird_[nmaxlayers_] = {0.};
    double diffEff_[nmaxlayers_] = {0.};

    for (int layer = 0; layer < layers_; layer++)
    {
        // calculate fractional reduction to the power of 1/3 for simpler use
        fracRedThird_[layer] = cbrt(1.0-fracRed_[i][layer]);
        //fracRedThird_[layer] = std::max(fracRedThird_[layer], SMALL);

        // Calculate the effective molecular diffusion
        effDiffBinary[i][layer] = molecularDiffusion_[i]*(porosity_[i][layer]/tortuosity_) + SMALL;

        // Calculate the knudsen diffusion
        effDiffKnud[i][layer]  =  (pore_diameter_/6.0)*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_))*(porosity_[i][layer]/tortuosity_) + SMALL;

        // total effective diffusivity
        diffEff_[layer] =   effDiffKnud[i][layer]*effDiffBinary[i][layer]/(effDiffBinary[i][layer]+effDiffKnud[i][layer]) + SMALL;
    }

    if (screenflag_ && screen)
    {
        fprintf(screen, "diffEff_[0]: %f, diffEff_[1]: %f, diffEff_[2]: %f \n fracRedThird_[0]: %f, fracRedThird_[1]: %f, fracRedThird_[2] : %f \n"
                        "fracRed_[0]: %f, fracRed_[1]: %f, fracRed_[2]: %f \n"
                ,diffEff_[0], diffEff_[1], diffEff_[2], fracRedThird_[0], fracRedThird_[1] , fracRedThird_[2],fracRed_[0][0], fracRed_[0][1], fracRed_[0][2]);
    }

    // calculation of diffusion term
    if (molecularDiffusion_[i] < SMALL)
    {
        for (int layer = 0; layer < layers_; layer++)
            Bterm[i][layer] = 0.0;
    }
    else
    {
        Bterm[i][0]   =   ((1-fracRedThird_[0])/fracRedThird_[0])*((radius_[i]/cg_)/(diffEff_[0]));
        for (int layer = 1; layer < layers_; layer++)
        {
            Bterm[i][layer] = (fracRedThird_[layer-1]-fracRedThird_[layer])/(fracRedThird_[layer-1]*fracRedThird_[layer])*((radius_[i]/cg_)/(diffEff_[layer]));
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getMassT(int i)
{
    updatePtrs();
    // initialize sherwood & schmidt numbers for every particle
    double Sc_[atom->nlocal] = {0.};
    double Sh_[atom->nlocal] = {0.};

    // if molecular diffusion is around 0, overwrite to avoid numerical errors.
    if (molecularDiffusion_[i] < SMALL)
        molecularDiffusion_[i] += SMALL;

    Sc_[i]  =   nuf_[i]/molecularDiffusion_[i];
    Sh_[i]  =   2.0+0.6*sqrt(Rep_[i])*cbrt(Sc_[i]);

    Massterm[i] = Sh_[i]*molecularDiffusion_[i]/(2.0*(radius_[i]/cg_)+SMALL);
    Massterm[i] = 1.0/(Massterm[i]);

    if (screenflag_ && screen)
        fprintf(screen, "Schmidt number: %f \n",Sc_[0]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction(int i, double *dmA_, double *x0_eq_)
{
    updatePtrs();
    double W = 0.;
    double p_eq_[nmaxlayers_] = {0.};
    double p_A = 0.;

    for (int layer = 0; layer < layers_; layer++)
    {
        p_eq_[layer] = x0_eq_[layer]*partP_[i];
    }
    p_A = xA_[i]*partP_[i];
    if (screenflag_ && screen)
            fprintf(screen, "p_eq_I: %f, p_eq_II: %f, p_eq_III: %f, p_A: %f \n", p_eq_[0], p_eq_[1],p_eq_[2],p_A);

    if (layers_ == nmaxlayers_)
    {
        // including reaction resistance and diffusion coeff terms
        W = (Aterm[i][2]+Bterm[i][2])*(Aterm[i][0]*(Aterm[i][1]+Bterm[i][1]+Bterm[i][0]+Massterm[i])+(Aterm[i][1]+Bterm[i][1])*(Bterm[i][0]+Massterm[i]))
                +Aterm[i][1]*(Aterm[i][0]*(Bterm[i][1]+Bterm[i][0]+Massterm[i])+Bterm[i][1]*(Bterm[i][0]+Massterm[i]));
        // hematite to magnetite
        dY[i][2]   =   ((Aterm[i][0]*(Aterm[i][1]+Bterm[i][1]+Bterm[i][0]+Massterm[i])+(Bterm[i][0]+Massterm[i])*(Aterm[i][1]+Bterm[i][1]))*(p_A-p_eq_[2])
                -   (Aterm[i][0]*(Bterm[i][1]+Bterm[i][0]+Massterm[i])+Bterm[i][1]*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[1])
                -   (Aterm[i][1]*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[0]))/W;

        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][2] < 0.0) dY[i][2] = 0.0;

        // magnetite to wustite
        dY[i][1]   =   (((Aterm[i][2]+Bterm[i][2]+Bterm[i][1])*(Aterm[i][0]+Bterm[i][0]+Massterm[i])+Aterm[i][0]*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[1])
                -   (Bterm[i][1]*(Aterm[i][0]+Bterm[i][0]+Massterm[i])+Aterm[i][0]*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[2])
                -   ((Aterm[i][2]+Bterm[i][2])*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[0]))/W;
        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][1] < 0.0) dY[i][1] = 0.0;

        // wustite to iron
        // if magnetite is not reducing, wustite does also not reduce
        if ((dY[i][1] == 0.0))
            dY[i][0] = 0.0;
        else
            dY[i][0]   =   (((Aterm[i][2]+Bterm[i][2])*(Aterm[i][1]+Bterm[i][1]+Bterm[i][0]+Massterm[i])+Aterm[i][1]*(Bterm[i][1]+Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[0])
                -   (Aterm[i][1]*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[2])
                -   ((Aterm[i][2]+Bterm[i][2])*(Bterm[i][0]+Massterm[i]))*(p_A-p_eq_[1]))/W;
        // reaction doesn't happen if chemical reaction rate is negative
        if (dY[i][0] < 0.0) dY[i][0] = 0.0;
        dY_previous3 = true;
    }
    else if (layers_ == 2)
    {
        W = (Aterm[i][1]+Bterm[i][1])*(Aterm[i][0]+Bterm[i][0]+Massterm[i])+Aterm[i][0]*(Bterm[i][0]+Massterm[i]);
        // hematite to magnetite
        dY[i][2]   =   0.0;

        // magnetite to wustite
        if (dY_previous3 == true)
        {
            if (dY[i][1] == 0.0)
                dY[i][1] = 0.0;
            else
                dY[i][1]   =   ((Aterm[i][0]+Bterm[i][0]+Massterm[i])*(p_A-p_eq_[1])-(Bterm[i][0]+Massterm[i])*(p_A-p_eq_[0]))/W;
        } else
            dY[i][1]   =   ((Aterm[i][0]+Bterm[i][0]+Massterm[i])*(p_A - p_eq_[1])-(Bterm[i][0]+Massterm[i])*(p_A-p_eq_[0]))/W;

        if (dY[i][1] < 0.0) dY[i][1] = 0.0;

        // wustite to iron
        if ((dY[i][1] == 0.0))
            dY[i][0] = 0.0;
        else
            dY[i][0]   =   ((Aterm[i][1]+Bterm[i][1]+Bterm[i][0]+Massterm[i])*(p_A - p_eq_[0])-(Bterm[i][0]+Massterm[i])*(p_A-p_eq_[1]))/W;
        if (dY[i][0] < 0.0) dY[i][0] = 0.0;

        dY_previous2 = true;
    }
    else if (layers_ == 1)
    {
        // rate of chemical reaction for 1 active layer
        W = Aterm[i][0]+Bterm[i][0]+Massterm[i];

        // hematite to magnetite
        dY[i][2]   =   0.0;

        // magnetite to wustite
        dY[i][1]   =   0.0;

        // wustite to iron
        if (dY_previous2 == true)
        {
            if (dY[i][0] == 0.0)
                dY[i][0] = 0.0;
            else
                dY[i][0] = (p_A - p_eq_[0])/W;
        }
        else
            dY[i][0]   =   (p_A - p_eq_[0])/W;

        if (dY[i][0] < 0.0) dY[i][0] = 0.0;
    }

    for (int j = 0 ; j < layers_; j++)
    {
        // mass flow rate for reactant gas species
        // dmA is a positive value
        dmA_[j] =   dY[i][j]*(1.0/(Runiv*T_[i]))*molMass_A_*(MY_4PI*((radius_[i]*radius_[i])/(cg_*cg_)))*TimeStep*nevery;
        // fix property added so values are otuputted to file
        dmA_f_[i][j] = dmA_[j];
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, double *dmA_)
{
    updatePtrs();
    if (screenflag_ && screen)
        fprintf(screen,"run update atom props \n");
    // based on material change: update relative radii, average density and mass of atom i
    // stoichiometric coefficients of reactions
    double v_reac_[3] = {1, 1, 3};
    double v_prod_[3] = {1, 3, 2};
    // initialize radius, mass change of layer and sum of particle
    double rad[nmaxlayers_+1] = {0.};
    double dmL_[nmaxlayers_+1] = {0.};     // mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass_p_new = 0.0;
    double scale_reduction_rate = 10.0;

    /* Mass Change of Layers */
    /* dmL is a positive value, therefore it will be subtracted from the total mass */
    // Fe2O3 (iniital)
    dmL_[layers_] = dmA_[layers_-1]*v_reac_[layers_-1]*(layerMolMasses_[layers_]/molMass_A_);
    // Initial FeO & Fe3O4 layers
    // changes with active layers
    for (int layer = 1; layer < layers_; layer++)
        dmL_[layer] = -dmA_[layer]*v_prod_[layer]*(layerMolMasses_[layer]/molMass_A_) + dmA_[layer-1]*v_reac_[layer-1]*(layerMolMasses_[layer]/molMass_A_);
    // Fe
    dmL_[0] = -dmA_[0]*v_prod_[0]*(layerMolMasses_[0]/molMass_A_);

    for (int layer = 0; layer < layers_; layer++) {
        dmL_[layer] *= scale_reduction_rate;
    }

    // New layer masses
    for (int j = 0; j <= layers_; j++)
    {
        massLayer_[i][j] -= dmL_[j];
        // Limit minimum mass layer to 1e-20
        massLayer_[i][j] = std::max(massLayer_[i][j], 1e-20);

        // calculate total mass of particle
        sum_mass_p_new    +=  massLayer_[i][j];
    }

    // Core layer radius (initial Fe2O3)
    rad[layers_] = cbrt((0.75*massLayer_[i][layers_])/(rhoeff_[i][layers_]*M_PI));     

    // Outer layer radii (Fe3O4, FeO)
    for (int layer = layers_ - 1; layer > 0; layer--)
        rad[layer]   =   cbrt((0.75*massLayer_[i][layer]/(rhoeff_[i][layer]*M_PI))+rad[layer+1]*rad[layer+1]*rad[layer+1]);

    // Iron Layer / Ore radius = constant.
    rad[0] = radius_[i]/cg_;

    // Determine new relative radii after reduction
    for (int j = 1; j <= layers_; j++)
    {
        relRadii_[i][j] = rad[j]/rad[0];
    }

    relRadii_[i][1] = std::min(0.9999999, relRadii_[i][1]);
    relRadii_[i][2] = std::min(0.9999998, relRadii_[i][2]);
    relRadii_[i][3] = std::min(0.9999997, relRadii_[i][3]);

    // Total mass of particle with coarse-graining
    pmass_[i]   =   sum_mass_p_new*cg_*cg_*cg_;

    // total particle effective density
    // pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);

    if(screenflag_ && screen)
        fprintf(screen, "radius_: %f, pmass_: %f, pdensity_: %f\n ", radius_[0], pmass_[0], pdensity_[0]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_gas_properties(int i, double *dmA_)
{
    /*double kch2_ = 0.0;
    kch2_ = xA_[i] + xC_[i];*/

    // based on material change: update gas-phase source terms for mass and heat
    for (int j = 0; j < nmaxlayers_;j++)
    {
        // Reactant gas mass change
        changeOfA_[i]   -=  dmA_[j];
        // Limit maximum reactant gas
        if (changeOfA_[i] > 0.0) changeOfA_[i] = 0.0;
        // Product gas mass change
        changeOfC_[i]   +=  dmA_[j]*molMass_C_/molMass_A_;
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
    /*double f_HM, f_MW, f_WF;

    f_WF = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    f_MW = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    f_HM = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));

    fracRed_[i][0] = f_WF;
    fracRed_[i][1] = f_MW;
    fracRed_[i][2] = f_HM; */
    updatePtrs();

    double f_HM = 0.;
    double f_MW = 0.;
    double f_WF = 0.;

    f_WF = 1.0 - relRadii_[i][1]*relRadii_[i][1]*relRadii_[i][1];
    f_MW = 1.0 - relRadii_[i][2]*relRadii_[i][2]*relRadii_[i][2];
    f_HM = 1.0 - relRadii_[i][3]*relRadii_[i][3]*relRadii_[i][3];

    fracRed_[i][0] = f_WF;
    fracRed_[i][1] = f_MW;
    fracRed_[i][2] = f_HM;
}

/* ---------------------------------------------------------------------- */
void FixChemShrinkCore::init_defaults()
{
    molMass_A_  =   molMass_C_  =  0.0;
    rhoeff_ = NULL;
    porosity_ = NULL;
    pore_diameter_ = tortuosity_ = 0.0;
    relRadii_ = massLayer_ = NULL;
    layerDensities_ = layerMolMasses_ = k0_ = Ea_ = NULL;
    xA_ = xC_ = NULL;

    // particle properties total
    radius_ = pmass_ = pdensity_ = NULL;

    // initialise fix handles
    changeOfA_ = changeOfC_ = T_ =  molecularDiffusion_ = nuf_ = Rep_ =  partP_ = Massterm = reactionHeat_ = NULL;
    Aterm = Bterm = effDiffBinary = effDiffKnud = fracRed_ =   NULL;

    dY = dmA_f_ = NULL;

    TimeStep = 0.0;

    fix_changeOfA_  = NULL;
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
    fix_layerMass_ = NULL;           //  [kg]
    fix_layerDens_ = NULL;           //  [kg/m^3]
    fix_layerMolMass_ = NULL;        //  [kg/mole]
    fix_k0_ = NULL;             //  [m/s]
    fix_Ea_ = NULL;             //  [J/mol] - [kg*m^2/s^2*mol]
    fix_porosity_ = NULL;       //  [%]
    fix_rhoeff_ = NULL;
    fix_tortuosity_ = NULL;
    fix_pore_diameter_ = NULL;   //  [m]
    fix_dY_ = NULL;
    fix_totalMole_      =   NULL;
    fix_dmA_ = NULL;
}
