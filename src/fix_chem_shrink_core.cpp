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

// mass equivalents used to calculate effective densities
#define q_Fe2O3_Fe3O4   0.966603
#define q_Fe3O4_FeO     0.9110
#define q_FeO_Fe        0.699426

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    nmaxlayers_(3),
    layers_(nmaxlayers_),
    rmin_(0.001),       //  [m]
    fix_changeOfA_(0),
    fix_changeOfC_(0),
    fix_rhogas_(0),
    fix_tgas_(0),
    fix_reactionHeat_(0),
    fix_diffcoeff_(0),
    fix_nuField_(0),
    fix_partRe_(0),
    fix_molefraction_(0),
    // for debug purposes, should be deleted afterwards
    fix_fracRed(0),
    //
    fix_layerRelRad_(0),
    fix_layerMass_(0),      //  [kg]
    fix_dens_(0),           //  [kg/m^3]
    fix_molMass_(0),        //  [kg/mole]
    fix_k0_(0),             //  [m/s]
    fix_Ea_(0),             //  [J/mol] - [kg*m^2/s^2*mol]
    //fix_porosity_(0),       //  [%]
    fix_tortuosity_(0),
    fix_pore_diameter_(0)  //  [m]
{
    if ((strncmp(style, "chem/shrink/core", 16) == 0) && ((!atom->radius_flag) || (!atom->rmass_flag)))
        error->all(FLERR, "Fix chem/shrink/core needs per particle radius and mass");

    // defaults    screenflag_ =   0;
    molMass_A_  =   molMass_C_  = 0;
    speciesA    =   speciesC    = NULL;
    pmass_      =   NULL;
    fc_         =   NULL;
    comm_established = false;
    iarg_ = 3;
    rhoeff_Fe2O3 = rhoeff_Fe3O4 = rhoeff_FeO = rhoeff_Fe = 0.;

    if (narg < 11)
        error->all(FLERR, "not enough arguments");

    bool hasargs = true;

    while (iarg_ < narg && hasargs)
    {
        if (strcmp(arg[iarg_], "speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error->fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            speciesA = new char[strlen(arg[iarg_ + 1])];
            strcpy(speciesA, arg[iarg_ + 1]);
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_], "molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_A_ = atof(arg[iarg_ + 1]);
            if (molMass_A_ < 0)
                error->fix_error(FLERR, this, "molar mass of A is not defined");
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_], "speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            speciesC = new char[strlen(arg[iarg_ + 1])];
            strcpy(speciesC, arg[iarg_ + 1]);
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_], "molMassC") == 0)
        {
            if (iarg_ + 2 > narg)
                error->fix_error(FLERR, this, "Wrong number of arguments");
            molMass_C_ = atof(arg[iarg_ + 1]);
            if (molMass_C_ < 0)
                error->fix_error(FLERR, this, "molar mass of C is not defined");
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_],"screen") == 0)
        {
            if (iarg_+2 > narg) error->all(FLERR,"Illegal fix/chem/shrink command");
            if (strcmp(arg[iarg_+1],"yes") == 0) screenflag_ = 1;
            else if (strcmp(arg[iarg_+1],"no") == 0) screenflag_ = 0;
            else error->all(FLERR,"Illegal fix/chem/shrink command");
            hasargs = true;
            iarg_ += 2;
        }
        else if (strcmp(arg[iarg_],"nevery") == 0)
        {
            nevery = atoi(arg[iarg_+1]);
            if (nevery <= 0) error->fix_error(FLERR,this,"");
            iarg_+=2;
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

    // define diffusant species
    diffA = new char [strlen(speciesA)+strlen("_diffCoeff")+1];
    strcpy(diffA,speciesA);
    strcat(diffA,"_diffCoeff");

    // reacting species bulk molar fraction
    moleFrac = new char [strlen("X_")+strlen(speciesA)+1];
    strcpy(moleFrac,"X_");
    strcat(moleFrac,speciesA);

    // nevery = 1;
    time_depend = 1;
    force_reneighbor = 1;
    next_reneighbor = update -> ntimestep + nevery;
    restart_global = 0;
    vector_flag =   1;
    size_vector =   4;
    extvector   =   0;
    local_flag  =   1;

    ts_create_  =   update->ntimestep;
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
    delete massA;
    delete massC;
    delete diffA;
    delete moleFrac;

    delete speciesA;
    delete speciesC;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_changeOfA_)     modify  ->  delete_fix(massA);
        if (fix_changeOfC_)     modify  ->  delete_fix(massC);
        if (fix_rhogas_)        modify  ->  delete_fix("partRho");
        if (fix_tgas_)          modify  ->  delete_fix("partTemp");
        if (fix_reactionHeat_)  modify  ->  delete_fix("reactionHeat");
        if (fix_diffcoeff_)     modify  ->  delete_fix(diffA);
        if (fix_nuField_)       modify  ->  delete_fix("partNu");
        if (fix_partRe_)        modify  ->  delete_fix("partRe");
        if (fix_molefraction_)  modify  ->  delete_fix(moleFrac);
        // fracRed is for debug purposes to see fractional reduction
        if (fix_fracRed)        modify  ->  delete_fix("fracRed");
        if (fix_layerRelRad_) modify->  delete_fix("relRadii");
        if (fix_layerMass_)   modify -> delete_fix("massLayer");
    }
}

/* ---------------------------------------------------------------------- */

int FixChemShrinkCore::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_create()
{
    if (comm -> me == 0 && screen)
        fprintf(screen,"do post_create \n");

    if (!fix_layerRelRad_)
    {
        const char* fixarg[12];
        fixarg[0]="layerRelRad";        // fixid
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="relRadii";           // propertyid
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="1.0";
        fixarg[9]="0.6";
        fixarg[10]="0.5";
        fixarg[11]="0.3";
        modify->add_fix(12,const_cast<char**>(fixarg));
        fix_layerRelRad_ =  static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii","property/atom","vector",atom->ntypes,0,style));
    }

    if (!fix_dens_)
    {
        char* fixname = new char [50];
        strcpy(fixname,"density_");
        strcat(fixname,group->names[igroup]);
        if (screenflag_ && screen)
            fprintf(screen, "density_ = %s \n", fixname);

        const char* fixarg[9];
        fixarg[0]="density";                // fixid
        fixarg[1]=group->names[igroup];
        fixarg[2]="property/global";
        fixarg[3]=fixname;           // propertyid
        fixarg[4]="vector";
        fixarg[5]="7870.0";
        fixarg[6]="5740.0";
        fixarg[7]="5170.0";
        fixarg[8]="5240.0";
        modify->add_fix(9,const_cast<char**>(fixarg),style);
        fix_dens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",atom->ntypes,0,style));
        delete []fixname;
    }

    // check for *fix_layerMass_
    if (!fix_layerMass_)
    {
        const char* fixarg[12];
        fixarg[0]="massLayer";
        // fixarg[1]=group->names[igroup];
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
        modify->add_fix(12,const_cast<char**>(fixarg));
        fix_layerMass_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer","property/atom","vector",atom->ntypes,0,style));
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::updatePtrs()
{
    changeOfA_      =   fix_changeOfA_  ->  vector_atom;
    changeOfC_      =   fix_changeOfC_  ->  vector_atom;
    rhogas_         =   fix_rhogas_     ->  vector_atom;
    T_              =   fix_tgas_       ->  vector_atom;
    reactionHeat_   =   fix_reactionHeat_-> vector_atom;
    // handle for diffusion coefficient
    molecularDiffusion_     =   fix_diffcoeff_  ->  vector_atom;
    // handles for viscosity field and reynolds number
    nuf_            =   fix_nuField_    ->  vector_atom;
    Rep_            =   fix_partRe_     ->  vector_atom;
    // handle for molar fraction
    X0_             =   fix_molefraction_   -> vector_atom;

    // handle for relative radius
    relRadii_       =   fix_layerRelRad_    -> array_atom;
    // mass handle
    massLayer_      = fix_layerMass_    ->  array_atom;
    // density handle
    layerDensities_ =   fix_dens_           -> get_values();
    // molar mass handle
    layerMolMasses_ = fix_molMass_   -> get_values();
    // chemical prop
    // frequency factor
    k0_             =  fix_k0_       -> get_values();
    // activation energy
    Ea_             =  fix_Ea_       -> get_values();
    // particle porosity properties
    //porosity_       =   fix_porosity_       -> get_values();
    tortuosity_     =   fix_tortuosity_     -> get_values();
    pore_diameter_  =   fix_pore_diameter_  -> get_values();

    pmass_          =   atom    -> rmass;
    pdensity_       =   atom    -> density;
    TimeStep        =   update  -> dt;
    radius_         =   atom    -> radius;

    //
    fracRed_        =   fix_fracRed ->array_atom;
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

    // find coupling fix & get coupling interval value
    fc_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    couple = fc_ -> couple_nevery_ + 1;

    int ntype = atom -> ntypes;
    int c = strlen(id) + 1;
    char* fixname = new char[c];

    // look up pre-exponential factor k0
    strcpy (fixname,"k0_");
    strcat(fixname,id);
    fix_k0_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname;

    // look up activation energies Ea
    fixname = new char [c];
    strcpy(fixname, "Ea_");
    strcat(fixname, id);
    fix_Ea_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname, "property/global", "vector", ntype, 0, "FixChemShrinkCore"));
    delete[]fixname;

    // Molar Mass
    fixname = new char [50];
    strcpy(fixname,"molMass_");
    strcat(fixname,group->names[igroup]);
    if (screenflag_ && screen)
        fprintf(screen, "molMass_ = %s \n", fixname);
    fix_molMass_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname;

    // references
    fix_changeOfA_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massA, "property/atom", "scalar", 0, 0, id));
    fix_changeOfC_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massC, "property/atom", "scalar", 0, 0, id));
    fix_rhogas_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partRho", "property/atom", "scalar", 0, 0, id));
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, style));
    fix_reactionHeat_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("reactionHeat", "property/atom", "scalar", 0, 0, id));
    fix_diffcoeff_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffA, "property/atom", "scalar", 0, 0, id));
    fix_nuField_        =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partNu", "property/atom", "scalar", 0, 0, style));
    fix_partRe_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partRe", "property/atom", "scalar", 0, 0, style));
    fix_molefraction_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property(moleFrac, "property/atom", "scalar", 0, 0, id));
    //
    fix_fracRed         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("fracRed", "property/atom", "vector", 0, 0, style));

    // lookup porosity,tortuosity and pore diameter values
    // fix_porosity_       =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("porosity_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));
    fix_tortuosity_     =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("tortuosity_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));
    fix_pore_diameter_  =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("pore_diameter_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));

    updatePtrs();

    for (int i = 0; i < atom->nlocal; i++)
    {
        calcMassLayer(i);

        rhoeff_Fe2O3 = atom->density[i];
        rhoeff_Fe3O4 = rhoeff_Fe2O3*q_Fe2O3_Fe3O4;
        rhoeff_FeO  =   rhoeff_Fe3O4*q_Fe3O4_FeO;
        rhoeff_Fe   =   rhoeff_FeO*q_FeO_Fe;
    }
    if (screenflag_ && screen)
    {
        fprintf(screen, "rhoeff_Fe2O3: %f \n", rhoeff_Fe2O3);
        fprintf(screen, "rhoeff_Fe3O4: %f \n", rhoeff_Fe3O4);
        fprintf(screen, "rhoeff_FeO: %f \n", rhoeff_FeO);
        fprintf(screen, "rhoeff_Fe: %f \n", rhoeff_Fe);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    updatePtrs();
    int i;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    double a_[nmaxlayers_];              // reaction resistance value for each layer
    double x0_eq_[nmaxlayers_];          // molar fraction of reactant gas
    double b_[nmaxlayers_];              // diffusion resistance value
    double dmA_[nmaxlayers_];            // mass flow rate of reactant gas species for each layer
    double masst_[nlocal];
    ts = update->ntimestep;

    // do chem/shrink/core calculations if communication between CFDEM and LIGGGHTS already happened
    // need initial values from CFDEM side
    if (!comm_established)
    {
        if (ts > ts_create_ + couple)
        {
            comm_established = true;
        }
    }

    if (comm_established)
    {
        for (i = 0; i < nlocal; i++)
        {
            if (mask[i] & groupbit)
            {
                active_layers(i);
                if (layers_ > 0)
                {
                    getXi(i,x0_eq_);
                    for (int j = 0; j < nmaxlayers_; j++)
                    {
                        a_[j] = 0.0;
                        dmA_[j] = 0.0;
                    }
                    getA(i,a_);
                    getB(i,b_);
                    // getMassT(i,masst_);
                    // reaction(i, a_, dmA_,x0_eq_, b_, masst_);
                    reaction2(i,a_,dmA_,x0_eq_, b_);
                    // reaction3(i,a_,dmA_,x0_eq_);
                    update_atom_properties(i,dmA_);
                    update_gas_properties(i,dmA_);
                    FractionalReduction(i);
                }
            }
        }
    }

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

    //NP tags and maps
    if (atom->molecular == 0) {
    int *tag = atom->tag;
    for (i = 0; i < atom->nlocal; i++) tag[i] = 0;
    atom->tag_extend();
    }

    if (atom->tag_enable) {
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }
}

/* ----------------- compute particle surface area ------------------------ */

double FixChemShrinkCore::partSurfArea(double radius)
{
     double A_p =   4*M_PI*radius*radius;
     return (A_p);
}

/* ---------------------------------------------------------------------- */

// returns number of active layers
// calculate the new masses of layers, (important if layer number has changed)
int FixChemShrinkCore::active_layers(int i)
{
    for(int j  = 1; j <= layers_; j++)
    {
        if (relRadii_[i][j]*radius_[i] < rmin_)
        {
            layers_ -= 1;
            calcMassLayer(i);
            break;
        }
    }
    if (screenflag_ && screen)
        fprintf(screen,"active layers: %i \n", layers_);

    return layers_;
}

/* ---------------------------------------------------------------------- */

// Determine the mass for every layer -- the latest layer represents the core of the particle.
void FixChemShrinkCore::calcMassLayer(int i)
{
    double rad = radius_[i];
    massLayer_[i][layers_]   =   1.33333*M_PI*relRadii_[i][layers_]*relRadii_[i][layers_]*relRadii_[i][layers_]*layerDensities_[layers_];
    massLayer_[i][layers_]   *= rad*rad*rad;
    for (int j=0; j < layers_; j++)
    {
        massLayer_[i][j]   =   1.33333*M_PI*((relRadii_[i][j]*relRadii_[i][j]*relRadii_[i][j])-(relRadii_[i][j+1]*relRadii_[i][j+1]*relRadii_[i][j+1]))*layerDensities_[j];
        massLayer_[i][j]   *= rad*rad*rad;
    }
}

/* ---------------------------------------------------------------------- */

// Calculate the equilibrium molar fractions
// And the bulk molar fraction
void FixChemShrinkCore::getXi(int i, double *x0_eq_)
{
    // calculate equilibrium concentration from K_eq(layer, T)
    double xp_eq_[nmaxlayers_];
    double kc = 0.4;

    for (int j = 0; j < layers_; j++)
    {
        xp_eq_[j]  =   kc*K_eq(j,T_[i])/(1+K_eq(j,T_[i]));
        x0_eq_[j]  =   kc/(1+K_eq(j,T_[i]));

        // for debug
        // ==========================================
        if (screenflag_ && screen)
        {
            fprintf(screen,"x0_eq : %f \n", x0_eq_[j]);
            fprintf(screen,"xp_eq : %f \n", xp_eq_[j]);
        }
        // ==========================================
    }

    // X0_[i] = std::min(X0_[i],kc);

    // for debug
    // ==========================================
    if (screenflag_ && screen)
        fprintf(screen,"x0_: %f \n", X0_[i]);
    // ==========================================
}

/* ---------------------------------------------------------------------- */

// calculate A - the chemical reaction resistance term - starting from wustite(0).
// Equation available in literature. (Valipour, Natsui, Nietrost...)
// A_[j] [s/m]
void FixChemShrinkCore::getA(int i, double *a_)
{
    for (int j = 0; j < layers_ ; j++)
    {
        // a_[j]   =   (1/(k0_[j]*exp(-Ea_[j]/(Runiv*T_[i]))))*(1/(1+1/K_eq(j,T_[i])))*(1/pow((1-fracRed_[i][j]),0.6666));
        a_[j]   =   (k0_[j]*exp(-Ea_[j]/(Runiv*T_[i])))*pow((1-fracRed_[i][j]),0.6666);
        a_[j]   =   1/a_[j];
    }

    if (screenflag_ && screen)
    {
        fprintf(screen, "reaction resistance term 1 : %f \n",a_[0]);
        fprintf(screen, "reaction resistance term 2 : %f \n",a_[1]);
        fprintf(screen, "reaction resistance term 3 : %f \n",a_[2]);
    }
}

/* ---------------------------------------------------------------------- */

// Calculate B - the diffusion resistance term - 0 = wustite
// Use binary diffusion for mixture, and knudsen diffusion to determine
// effective diffusion coefficient.
void FixChemShrinkCore::getB(int i, double *b_)
{
    if (comm -> me == 0 && screen)
        fprintf(screen,"DO GET B FUNCTION!! \n");

    double BinaryDiff_= 0.0;
    double KnudsenDiff_= 0.0;
    double diffTotal_ = 0.0;
    double fracRedThird_[nmaxlayers_];
    diffEff_ = new double [nmaxlayers_];

    for (int j = 0; j < layers_; j++)
    {
        // calculate fractional reduction to the power of 1/3 for simpler use
        fracRedThird_[j] = pow((1-fracRed_[i][j]),THIRD);
    }

    // Calculate the effective molecular diffusion
    BinaryDiff_ = molecularDiffusion_[i];

    // Due to resistance equation define effective molecular diffusion as 1 over.
    // Make sure that it is not divided by zero
    if (molecularDiffusion_[i] != 0.0) // or (BinaryDiff_[i] != 0)
        BinaryDiff_   =   1.0/BinaryDiff_;

    // Calculate the knudsen diffusion
    KnudsenDiff_  =   pore_diameter_[i]/6.*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_));

    // Due t resistance equation define the Knudsen diffusion as 1 over
    KnudsenDiff_   =   1.0/KnudsenDiff_;

    diffTotal_   = BinaryDiff_ + KnudsenDiff_;

    // debug options -- print out values for effective molecular diffusion to screen
    if (screenflag_ && screen)
    {
        fprintf(screen,"dCoeff value is : %f \n", molecularDiffusion_[i]);
        fprintf(screen,"binary diff 0: %f \n",BinaryDiff_);
        fprintf(screen,"knudsen diff: %f \n",KnudsenDiff_);
        fprintf(screen,"layerPorosity 0: %f \n",layerPorosity_(0));
        fprintf(screen,"layerPorosity 1: %f \n",layerPorosity_(1));
        fprintf(screen,"layerPorosity 2: %f \n",layerPorosity_(2));
    }

    // total effective diffusivity
    // Eq. : 1/D_i,j = 1/D_eff_binary + 1/D_eff_knudsen
    // diffEff_[i] = 1/D_i,j

    for (int j = 0; j < layers_ ; j++)
    {
        diffEff_[j] = (layerPorosity_(j)/tortuosity_[i])*(1/diffTotal_);
    }

    if (screenflag_ && screen)
    {
        fprintf(screen,"eff. diff 0: %f \n",diffEff_[0]);
        fprintf(screen,"eff. diff 1: %f \n",diffEff_[1]);
        fprintf(screen,"eff. diff 2: %f \n",diffEff_[2]);
    }

    // calculation of Diffustion Term from Valipour 2009
    // from wustite to iron
    b_[0]   =   (1-fracRedThird_[0])/fracRedThird_[0]*(radius_[i]/diffEff_[0]);

    for (int j = 1; j < layers_; j++)
    {
        // diffEff is multiplied here, in normal equation b=... radius/effectiveDiffusion
        // b_[j]   =   (relRadii_[i][j] - relRadii_[i][j+1])/(relRadii_[i][j]*relRadii_[i][j+1])*(radius_[i]*diffEff_[i]);
        b_[j] = (fracRedThird_[j-1]-fracRedThird_[j])/(fracRedThird_[j-1]*fracRedThird_[j])*(radius_[i]/diffEff_[j]);
    }

    if (screenflag_ && screen)
    {
        fprintf(screen,"diffusion constant [0]: %f \n",b_[0]);
        fprintf(screen,"diffusion constant [1]: %f \n",b_[1]);
        fprintf(screen,"diffusion constant [2]: %f \n",b_[2]);
    }

    // effecitve binary diffusion coefficient
    // BinaryDiff_[i] = molecularDiffusion_[i]*porosity_[i]/tortuosity_[i];

    // if (screenflag_ && screen)
    // {
    //     fprintf(screen,"eff. binary diff: %f \n",BinaryDiff_[i]);
    //    fprintf(screen,"dCoeff value is : %f \n", molecularDiffusion_[i]);
    // }

    //if (molecularDiffusion_[i] != 0.0) // or (BinaryDiff_[i] != 0)
    //{
    //    BinaryDiff_[i]   =   1.0/BinaryDiff_[i];
    //}

    // Knudsen diff equation is either
    // D_{i,k} = dp/3*sqrt((8*R*T_[i])/(M_PI*molMass_A_))
    // or simplified as
    // D_{i,k} = 4850*dp*sqrt(T_[i]/molMass_A_) dp supposed to be in cm
    // we use si units so convert from meter to cm!

    // KnudsenDiff_[i]   =  48.51*pore_diameter_[i]*sqrt(T_[i]/molMass_A_)*porosity_[i]/tortuosity_[i]; // [m^2/s]
    // KnudsenDiff_[i]  =   pore_diameter_[i]/6.*sqrt((8*Runiv*T_[i])/(M_PI*molMass_A_))*porosity_[j]/tortuosity_[i];
    // KnudsenDiff_[i]   =   1.0/KnudsenDiff_[i];
    // if (screenflag_ && screen)
    // {
    //     fprintf(screen,"eff. knudsen diff: %f \n",KnudsenDiff_[i]);
    // }

    // total effective diffusivity
    // Eq. : 1/D_i,j = 1/D_eff_binary + 1/D_eff_knudsen
    // diffEff_[i] = 1/D_i,j
    // diffEff_[i] = BinaryDiff_[i] + KnudsenDiff_[i];

    // if (screenflag_ && screen)
    //    fprintf(screen,"eff. diff: %f \n",diffEff_[i]);

    // calculation of Diffustion Term from Valipour 2009
    // from wustite to iron
    // diffEff_[i] = 1/D_i,j thus instead of dividing radius by diffEff we multiply it.
    // b_[0]   =   (1-relRadii_[i][1])/(relRadii_[i][1])*(radius_[i]*diffEff_[i]);
    //b_[0]   =   (1-fracRedThird_[0])/fracRedThird_[0]*radius_[i]*diffEff_[i];

    /*for (int j = 1; j < layers_; j++)
    {
        // diffEff is multiplied here, in normal equation b=... radius/effectiveDiffusion
        // b_[j]   =   (relRadii_[i][j] - relRadii_[i][j+1])/(relRadii_[i][j]*relRadii_[i][j+1])*(radius_[i]*diffEff_[i]);
        b_[j] = (fracRedThird_[j-1]-fracRedThird_[j])/(fracRedThird_[j-1]*fracRedThird_[j])*radius_[i]*diffEff_[i];
    }

    if (screenflag_ && screen)
    {
        fprintf(screen,"diffusion constant [0]: %f \n",b_[0]);
        fprintf(screen,"diffusion constant [1]: %f \n",b_[1]);
        fprintf(screen,"diffusion constant [2]: %f \n",b_[2]);
    } */

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getMassT(int i, double *masst_)
{
  /*  // initialize sherwood & schmidt numbers for every particle
    double Sc_[atom->nlocal];
    double Sh_[atom->nlocal];

    Sc_[i]  =   nuf_[i]/diffEff_[i];
    Sh_[i]  =   2+0.6*pow(Rep_[i],0.5)*pow(Sc_[i],THIRD);

    masst_[i] = Sh_[i]*diffEff_[i]/(2*radius_[i]);
    masst_[i] = 1/masst_[i];

    if (screenflag_ && screen)
    {
        fprintf(screen, "Schmidt number: %f \n",Sc_[i]);
        fprintf(screen, "Sherwood number: %f \n",Sh_[i]);
        fprintf(screen, "masst: %f \n",masst_[i]);
    } */
}
/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction(int i, double *a_, double *dmA_, double *x0_eq_, double *b_, double *masst_)
{
    updatePtrs();
    double W;
    double dY[nmaxlayers_];

    if (layers_ == nmaxlayers_)
    {
        // including reaction resistance and diffusion coeff terms
        W = (a_[2]+b_[2])*(a_[0]*(a_[1]+b_[1]+b_[0]+masst_[i])+(a_[1]+b_[1])*(b_[0]+masst_[i]))
                +a_[1]*(a_[0]*(b_[1]+b_[0]+masst_[i])+b_[1]*(b_[0]+masst_[i]));
        // hematite to magnetite
        dY[2]   =   ((a_[0]*(a_[1]+b_[1]+b_[0]+masst_[i])+(b_[0]+masst_[i])*(a_[1]+b_[1]))*(X0_[i]-x0_eq_[2])
                -   (a_[0]*(b_[1]+b_[0]+masst_[i])+b_[1]*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[1])
                -   (a_[1]*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[0]))*1/W;
        // magnetite to wustite
        dY[1]   =   (((a_[2]+b_[2]+b_[1])*(a_[0]+b_[0]+masst_[i])+a_[0]*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[1])
                -   (b_[1]*(a_[0]+b_[0]+masst_[i])+a_[0]*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   (((a_[2]+b_[2])*(a_[1]+b_[1]+b_[0]+masst_[i])+a_[1]*(b_[1]+b_[0]+masst_[i]))*(X0_[i]-x0_eq_[0])
                -   (a_[1]*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*(b_[0]+masst_[i]))*(X0_[i]-x0_eq_[1]))*1/W;
    }
     else if (layers_ == 2)
     {
        W = (a_[1]+b_[1])*(a_[0]+b_[0]+masst_[i])+a_[0]*(b_[0]+masst_[i]);

        // hematite to magnetite
        dY[2]   =   0.0;
        // magnetite to wustite
        dY[1]   =   ((a_[0]+b_[0]+masst_[i])*(X0_[i]-x0_eq_[1])-(b_[0]+masst_[i])*(X0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   ((a_[1]+b_[1]+b_[0]+masst_[i])*(X0_[i]-x0_eq_[0])-(b_[0]+masst_[i])*(X0_[i]-x0_eq_[1]))*1/W;
     }
     else if (layers_ == 1)
     {

        // rate of chemical reaction for 1 active layer
        W = a_[0]+b_[0]+masst_[i];

        // hematite to magnetite
        dY[2]   =   0.0;
        // magnetite to wustite
        dY[1]   =   0.0;
        // wustite to iron
        dY[0]   =   (X0_[i] - x0_eq_[0])*1/W;
     }

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        // mass flow rate for reactant gas species
        dmA_[j] =   dY[j]*rhogas_[i]*partSurfArea(radius_[i])*TimeStep*nevery;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction2(int i, double *a_, double *dmA_, double *x0_eq_, double *b_) // double masst
{
    updatePtrs();
    double W;
    double dY[nmaxlayers_];

    if (layers_ == nmaxlayers_)
    {
        // including reaction resistance and diffusion coeff terms
        W = (a_[2]+b_[2])*(a_[0]*(a_[1]+b_[1]+b_[0])+(a_[1]+b_[1])*b_[0])
                +a_[1]*(a_[0]*(b_[1]+b_[0])+b_[1]*b_[0]);
        // hematite to magnetite
        dY[2]   =   ((a_[0]*(a_[1]+b_[1]+b_[0])+b_[0]*(a_[1]+b_[1]))*(X0_[i]-x0_eq_[2])
                -   (a_[0]*(b_[1]+b_[0])+b_[1]*b_[0])*(X0_[i]-x0_eq_[1])
                -   (a_[1]*b_[0])*(X0_[i]-x0_eq_[0]))*1/W;
        // magnetite to wustite
        dY[1]   =   (((a_[2]+b_[2]+b_[1])*(a_[0]+b_[0])+a_[0]*b_[0])*(X0_[i]-x0_eq_[1])
                -   (b_[1]*(a_[0]+b_[0])+a_[0]*b_[0])*(X0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*b_[0])*(X0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   (((a_[2]+b_[2])*(a_[1]+b_[1]+b_[0])+a_[1]*(b_[1]+b_[0]))*(X0_[i]-x0_eq_[0])
                -   (a_[1]*b_[0])*(X0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*b_[0])*(X0_[i]-x0_eq_[1]))*1/W;
    }
     else if (layers_ == 2)
     {
        W = (a_[1]+b_[1])*(a_[0]+b_[0])+a_[0]*b_[0];

        // hematite to magnetite
        dY[2]   =   0.0;
        // magnetite to wustite
        dY[1]   =   ((a_[0]+b_[0])*(X0_[i]-x0_eq_[1])-b_[0]*(X0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   ((a_[1]+b_[1]+b_[0])*(X0_[i]-x0_eq_[0])-b_[0]*(X0_[i]-x0_eq_[1]))*1/W;
     }
     else if (layers_ == 1)
     {

        // rate of chemical reaction for 1 active layer
        W = a_[0]+b_[0];

        // hematite to magnetite
        dY[2]   =   0.0;
        // magnetite to wustite
        dY[1]   =   0.0;
        // wustite to iron
        dY[0]   =   (X0_[i] - x0_eq_[0])*1/W;
     }

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        // mass flow rate for reactant gas species
        // scaled up Timestep to match 0.01;
        dmA_[j] =   dY[j]*rhogas_[i]*partSurfArea(radius_[i])*TimeStep*nevery;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction3(int i, double *a_, double *dmA_, double *x0_eq_)
{
    updatePtrs();
    double W;
    double dY[nmaxlayers_];

    if (layers_ == nmaxlayers_)
    {
        // if only with reaction resistance term
        // rate of chemical reaction for only chemical reaction resistance.
        W = a_[2]*a_[1]*a_[0];
        // hematite to magnetite
        dY[2]   =   ((a_[0]*a_[1])*(X0_[i]-x0_eq_[2]))*1/W;
        // magnetite to wustite
        dY[1]   =   ((a_[2]*a_[0])*(X0_[i]-x0_eq_[1]))*1/W;
        // wustite to iron
        dY[0]   =   ((a_[2]*a_[1])*(X0_[i]-x0_eq_[0]))*1/W;
    }
     else if (layers_ == 2)
     {
         // rate of chemical reaction for 2 active layers
         W = (a_[1]*a_[0]);

         // hematite to magnetite
         dY[2]   =   0.0;
         // magnetite to wustite
         dY[1]   =   (a_[0]*(X0_[i]-x0_eq_[1]))*1/W;
         // wustite to iron
         dY[0]   =   (a_[1]*(X0_[i] - x0_eq_[0]))*1/W;
     }
     else if (layers_ == 1)
     {
         // rate of chemical reaction for 1 active layer
         W = a_[0];

         // hematite to magnetite
         dY[2]   =   0.0;
         // magnetite to wustite
         dY[1]   =   0.0;
         // wustite to iron
         dY[0]   =   (X0_[i] - x0_eq_[0])*1/W;
     }

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        dmA_[j] =   dY[j]*rhogas_[i]*partSurfArea(radius_[i])*TimeStep*nevery;     }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, double *dmA_)
{
    // based on material change: update relative radii, average density and mass of atom i        //  mass of each layer
    double dmB_[nmaxlayers_+1];               //  mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass_new = 0;                  // just for control case --- delete afterwards
    double vBr_[3] = {1, 1, 3};
    double vBp_[3] = {1, 3, 2};
    double radLayer_[nmaxlayers_+1];

    // double diff_radius;
    // Calculate mass flow rates (dmB[j])
    // Keep in mind the stoichiometric coefficients of reaction prod - reactants
    // (stochiometric ratio is always positive since reactant and product stoichiometry does not have different signs)
    dmB_[layers_] = dmA_[layers_-1]*vBr_[layers_-1]*(layerMolMasses_[layers_]/molMass_A_);
    for (int j = layers_ -1; j>= 1; j--)
    {
        dmB_[j] = - dmA_[j]*vBp_[j]*(layerMolMasses_[j]/molMass_A_) + dmA_[j-1]*vBr_[j-1]*(layerMolMasses_[j]/molMass_A_);
    }
    dmB_[0] = - dmA_[0]*vBp_[0]*(layerMolMasses_[0]/molMass_A_);

    if (layers_ == 2)
        dmB_[layers_+1] = 0.0;
    else if (layers_ == 1)
    {
        dmB_[layers_+2] = 0.0;
        dmB_[layers_+1] = 0.0;
    }

    // print out mass transfer rates for particle layers
    if (screenflag_ && screen)
    {
        fprintf(screen,"dmB[0]: %f \n",dmB_[0]);
        fprintf(screen,"dmB[1]: %f \n",dmB_[1]);
        fprintf(screen,"dmB[2]: %f \n",dmB_[2]);
        fprintf(screen,"dmB[3]: %f \n",dmB_[3]);
    }

    // based on chem. reactions, get new layer masses
    for (int j = 0; j <= layers_; j++)
    {
        // TL: massLayer_[j] should never be negative anyways; just make sure not more is removed than present
        massLayer_[i][j] -= dmB_[j];
        if (massLayer_[i][j] < SMALL)
                 massLayer_[i][j] = SMALL;

        sum_mass_new    +=  massLayer_[i][j];
    }

    // from the new layer masses, get new total mass, layer radii and total density
    pmass_[i]   =   sum_mass_new;

    radLayer_[layers_]   =   pow((0.75*massLayer_[i][layers_]/(layerDensities_[layers_]*M_PI)),THIRD);
    for (int j = layers_-1; j >= 0 ; j--)
    {
        radLayer_[j]   =   pow((0.75*massLayer_[i][j]/(M_PI*layerDensities_[j])+radLayer_[j+1]*radLayer_[j+1]*radLayer_[j+1]),THIRD);
    }
    radius_[i] = radLayer_[0];

    //detemine the new relative layer radii
    for (int j = 0; j <= layers_; j++)
    {
        relRadii_[i][j] = radLayer_[j]/radius_[i];
    }

    // calculate particle total density
    // depending on the new values
    // update total partilce density
    pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);

}

 // WARNING 2: if Wustite layer thickness < rmin AND T < 567 C, direct conversion of Magnetite to Fe is possible
 //            this means that Wu-Fe radius is shrinking together with the Ma-Wu radius
 //            in this case, the equilibrium value and reaction parameters at the Ma layer need to be adapted within
 //            IMPLEMENTATION ONLY IF NECESSARY */

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_gas_properties(int i, double *dmA_)
{
   // based on material change: update gas-phase source terms for mass and heat
    for (int j = 0; j < nmaxlayers_;j++)
    {
        changeOfA_[i]   -=  dmA_[j];
        changeOfC_[i]   +=  dmA_[j]*molMass_C_/molMass_A_;
    }
}

/* ---------------------------------------------------------------------- */

double FixChemShrinkCore::K_eq(int layer, double T)
{
    if (comm -> me == 0 && screen)
        fprintf(screen,"CALCULATE K_EQUILIBIRUM \n");
    // ATTENTION: K is usually given for partial pressures, might be necessary to convert to molar concentrations
    // for the reactions under consideration, this should not be the case (double-check !!!)
    // 0 = wustite, 1 = mangetite, 2 = hematite;
    double Keq_ = 0.;
    if(strcmp(speciesA,"CO")==0)
     {
         if (layer == 0)
             Keq_   =   pow(10,(917/T-1.097));
         else if (layer == 1)
             Keq_   =   pow(10,(-1834/T+2.17));
         else if (layer == 2)
             Keq_   =   exp(3968.37/T+3.94);
     }
     else if(strcmp(speciesA,"H2")==0)
     {
         if (layer == 0)
             Keq_   =   pow(10,(-827/T+0.468));
         else if (layer == 1)
             Keq_   =   pow(10,(-3577/T+3.74));
         else if (layer == 2)
             Keq_   =   exp(-362.6/T+10.344);
     }
    else
     {
         printf("Error : Undefined Reaction \n");
     }
    return Keq_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::FractionalReduction(int i)
{
    double f_[nmaxlayers_];
    // calculate the fractional reduction as defined in Tang et al. (2012)
    // "Simulation study on performance of Z-path Moving-fluidized Bed for Gaseous Reduction
    // of Iron Ore Fines" ISIJ International, Vol. 52, No. 7, pp. 1241 - 1249
    f_[0] = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    f_[1] = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));
    f_[2] = 1 - ((2*massLayer_[i][3]/layerMolMasses_[3])/(2*massLayer_[i][3]/layerMolMasses_[3]+3*massLayer_[i][2]/layerMolMasses_[2]+massLayer_[i][1]/layerMolMasses_[1]+massLayer_[i][0]/layerMolMasses_[0]));

    for (int k=0;k<nmaxlayers_;k++)
    {
        fracRed_[i][k] = f_[k];
    }
}


/* ---------------------------------------------------------------------- */
double FixChemShrinkCore::layerPorosity_(int layer)
{
    double eps_ = 0.;
    if (layer == 2)
        eps_ = 1 - rhoeff_Fe3O4/layerDensities_[2];
    if (layer == 1)
        eps_ = 1 - rhoeff_FeO/layerDensities_[1];
    if (layer == 0)
        eps_ = 1 - rhoeff_Fe/layerDensities_[0];
    return eps_;
}

