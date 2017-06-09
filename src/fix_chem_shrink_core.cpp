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


using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL 1e-10
#define LARGE 1e+10

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    nmaxlayers_(3),
    layers_(nmaxlayers_),
    rmin_(0.001),
    fix_k0_(0), // [m/s]
    fix_molMass_(0),
    fix_Ea_(0), // [J/mol]
    fix_dens_(0),
    fix_layerRelRad_(0),
    fix_layerMass_(0),
    fix_porosity_(0),
    fix_tortuosity_(0),
    fix_pore_diameter_(0)
{
    if (strncmp(style, "chem/shrink/core", 11) == 0 && (!atom->radius_flag) || (!atom->rmass_flag))
        error->all(FLERR, "Fix chem/shrink needs particle radius and mass");

    // defaults
    fix_concA_          =   NULL;
    fix_concC_          =   NULL;
    fix_changeOfA_      =   NULL;
    fix_changeOfC_      =   NULL;
    fix_rhogas_         =   NULL;
    fix_tgas_           =   NULL;
    fix_totalmole_      =   NULL;
    fix_reactionHeat_   =   NULL;
    fix_pressure_       =   NULL;
    fix_diffcoeff_      =   NULL;

    pgas_ = 0;
    molMass_A_ = molMass_C_ = 0;
    Runiv       =   8.3144;  // [J/molK]

    iarg_ = 3;
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
            if (strlen(speciesA) < 1)
                error->fix_error(FLERR, this, "speciesA is not defined");
            molMass_A_ = atof(arg[iarg_ + 1]);
            if (molMass_A_ < 1)
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
            if (strlen(speciesC) < 1)
                error->fix_error(FLERR, this, "speciesC not defined");
            molMass_C_ = atof(arg[iarg_ + 1]);
            if (molMass_C_ < 1)
                error->fix_error(FLERR, this, "molar mass of C is not defined");
            hasargs = true;
            iarg_ += 2;
        }
    }

    // define changed species mass A
    int x = 16;
    char cha[30];
    massA = new char[x];
    strcpy(cha, "Modified_");
    strcat(cha, speciesA);
    strcpy(massA, cha);

    // define changed species mass C
    massC = new char[x];
    strcpy(cha, "Modified_");
    strcat(cha, speciesC);
    strcpy(massC, cha);

    // define diffusant species
    diffA = new char[x];
    strcpy(cha, speciesA);
    strcat(cha, "_diffCoeff");
    strcpy(diffA, cha);

    radLayer_ = new double*[atom->nlocal];
    massLayer_ = new double*[atom->nlocal];
    for (int i = 0; i < atom->nlocal; i++)
    {
        radLayer_[i]    =   new double [nmaxlayers_+1];
        massLayer_[i]    =   new double [nmaxlayers_+1];
    }

    nevery = 1;
    time_depend = 1;
    force_reneighbor = 1;
    next_reneighbor = update -> ntimestep + nevery;
    restart_global = 0;
    vector_flag =   1;
    size_vector =   4;
    extvector   =   0;
    local_flag  =   1;
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
    delete massA;
    delete massC;
    delete diffA;

    delete speciesA;
    delete speciesC;

    for (int i=0; i < atom -> nlocal; i++)
    {
        delete[] radLayer_[i];
        delete[] massLayer_[i];
    }
    delete [] radLayer_;
    delete [] massLayer_;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_concA_)         modify  ->  delete_fix(speciesA);
        if (fix_concC_)         modify  ->  delete_fix(speciesC);
        if (fix_rhogas_)        modify  ->  delete_fix("partRho");
        if (fix_tgas_)          modify  ->  delete_fix("partTemp");
        if (fix_totalmole_)     modify  ->  delete_fix("partN");
        if (fix_reactionHeat_)  modify  ->  delete_fix("reactionHeat");
        if (fix_changeOfA_)     modify  ->  delete_fix(massA);
        if (fix_changeOfC_)     modify  ->  delete_fix(massC);
        if (fix_pressure_)      modify  ->  delete_fix("partP");
        if (fix_diffcoeff_)     modify  ->  delete_fix(diffA);

        // fixes from shrink/core property/globals will not be deleted
        if (fix_layerRelRad_) modify->  delete_fix("relRadii");
        // if (fix_layerMass_) modify -> delete_fix("massLayer");
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

    // TL:
    // look for fix_layerMass; if not there, created it and initialize its values according to the initial layer radii
    // make sure that massLayer_ is connected to that fix so that local changes affect it and other instances (second reaction, e.g. with H2)
    // use the proper values
    // fix_layerRelRad_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii","property/atom","vector",ntype,0,this->style,false));
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
        char* fixname = new char [16];
        strcpy(fixname,"density_");
        strcat(fixname,group->names[igroup]);
        if (screen)
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

    updatePtrsCreate();

    for (int i = 0 ; i < atom->nlocal; i++)
    {
        // layerRad(i);
        calcMassLayer(i);

        /* // check for *fix_layerMass_
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
            char arg8[30];
            sprintf(arg8,"%f",massLayer_[i][0]);
            fixarg[8]=arg8;
            char arg9[30];
            sprintf(arg9,"%f",massLayer_[i][1]);
            fixarg[9]=arg9;
            char arg10[30];
            sprintf(arg10,"%f",massLayer_[i][2]);
            fixarg[10]=arg10;
            char arg11[30];
            sprintf(arg11,"%f",massLayer_[i][3]);
            fixarg[11]=arg11;
            modify->add_fix(12,const_cast<char**>(fixarg),style);
            fix_layerMass_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("massLayer","property/atom","vector",atom->ntypes,0,style));
        }*/
    }
}

/* ---------------------------------------------------------------------- */
void FixChemShrinkCore::updatePtrsCreate()
{
    relRadii_       =   fix_layerRelRad_    -> array_atom;
    layerDensities_ =   fix_dens_           -> get_values();
    radius_     =   atom ->  radius;
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::updatePtrs()
{
    changeOfA_      =   fix_changeOfA_  ->  vector_atom;
    changeOfC_      =   fix_changeOfC_  ->  vector_atom;
    rhogas_         =   fix_rhogas_     ->  vector_atom;
    concA_          =   fix_concA_      ->  vector_atom;
    concC_          =   fix_concC_      ->  vector_atom;
    T_              =   fix_tgas_       ->  vector_atom;
    N_              =   fix_totalmole_  ->  vector_atom;
    reactionHeat_   =   fix_reactionHeat_-> vector_atom;
    pgas_           =   fix_pressure_    -> vector_atom;
    // diffusion coefficient
    dCoeff_         =   fix_diffcoeff_  ->  vector_atom;

    // material prop
    // radius handle & density handle
    // relRadii_       =   fix_layerRelRad_    -> array_atom;
    updatePtrsCreate();
    // mass handle
    // massLayer_      = fix_layerMass_    ->  array_atom;
    // molar mass handle
    layerMolMasses_ = fix_molMass_   -> get_values();
    // chemical prop
    // frequency factor
    k0_             =  fix_k0_       -> get_values();
    // activation energy
    Ea_             =  fix_Ea_       -> get_values();
    // particle porosity properties
    porosity_       =   fix_porosity_       -> get_values();
    tortuosity_     =   fix_tortuosity_     -> get_values();
    pore_diameter_  =   fix_pore_diameter_  -> get_values();

    pmass_          =   atom ->  rmass;
    pdensity_       =   atom -> density;
    TimeStep        =   update -> dt;
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

    int ntype = atom -> ntypes;
    int c = strlen(id) + 3;
    int x = 16;
    char* fixname = new char[c];

    // look up reaction properties
    // they are defined by the user via FixPropertyGlobal (as vector with 3 entries)
    // the fixes' names have to be chosen such that they can be identified by the reaction fix they correspond to
    // example:
    // fix OreReductionCO all chem/shrink/core ...
    // fix k0_OreReductionCO all property/global ...

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

    // lookup porosity,tortuosity and pore diameter values
    fix_porosity_   =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("porosity_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));
    fix_tortuosity_   =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("tortuosity_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));
    fix_pore_diameter_   =   static_cast<FixPropertyGlobal*>(modify->find_fix_property("pore_diameter_", "property/global", "scalar", 0, 0, "FixChemShrinkCore"));

    // Lookup material properties
    // they are defined by the user via FixPropertyGlobal (as vector with 4 entries)
    // the fixes' names have to be chosen such that they can be identified by the reaction fix they correspond to and the group it acts on
    // example:
    // fix OreReductionCO all chem/shrink/core ...
    // fix molMass_all all property/global ...
    // or
    // group Ore ...
    // fix OreReductionCO Ore chem/shrink/core ...
    // fix molMass_Ore all property/global ...

    // density and molar Mass
    /*fixname = new char [x];
    strcpy(fixname,"density_");
    strcat(fixname,group->names[igroup]);
    if (screen)
        fprintf(screen, "density_ = %s \n", fixname);
    fix_dens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname; */

    fixname = new char[x];
    strcpy(fixname,"molMass_");
    strcat(fixname,group->names[igroup]);
    if (screen)
        fprintf(screen, "molMass_ = %s \n", fixname);
    fix_molMass_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname;

    // look up *fix_layerRelRad_;
    // fix_layerRelRad_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii","property/atom","vector",ntype,0,"FixChemShrinkCore"));

    // references
    fix_concA_          =   static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesA, "property/atom", "scalar", 0, 0, id));
    fix_concC_          =   static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesC, "property/atom", "scalar", 0, 0, id));
    fix_changeOfA_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massA, "property/atom", "scalar", 0, 0, id));
    fix_changeOfC_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(massC, "property/atom", "scalar", 0, 0, id));
    fix_rhogas_         =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partRho", "property/atom", "scalar", 0, 0, id));
    fix_tgas_           =   static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, style));
    fix_totalmole_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partN","property/atom","scalar",0,0,style));
    fix_reactionHeat_   =   static_cast<FixPropertyAtom*>(modify->find_fix_property("reactionHeat", "property/atom", "scalar", 0, 0, id));
    fix_pressure_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partP","property/atom","scalar",0,0,style));
    fix_diffcoeff_      =   static_cast<FixPropertyAtom*>(modify->find_fix_property(diffA, "property/atom", "scalar", 0, 0, id));

    updatePtrs();

    if (comm -> me == 0 && screen)
    {
        fprintf(screen, "k0_[0] = %f \n", k0_[0]);
        fprintf(screen, "k0_[1] = %f \n", k0_[1]);
        fprintf(screen, "k0_[2] = %f \n", k0_[2]);
        fprintf(screen, "Ea_[0] = %f \n", Ea_[0]);
        fprintf(screen, "Ea_[1] = %f \n", Ea_[1]);
        fprintf(screen, "Ea_[2] = %f \n", Ea_[2]);
        fprintf(screen, "layerdens_[0] = %f \n",layerDensities_[0]);
        fprintf(screen, "layerdens_[1] = %f \n",layerDensities_[1]);
        fprintf(screen, "layerdens_[2] = %f \n",layerDensities_[2]);
        fprintf(screen, "layerdens_[3] = %f \n",layerDensities_[3]);
        fprintf(screen, "molMass_[0] = %f \n",layerMolMasses_[0]);
        fprintf(screen, "molMass_[1] = %f \n",layerMolMasses_[1]);
        fprintf(screen, "molMass_[2] = %f \n",layerMolMasses_[2]);
        fprintf(screen, "molMass_[3] = %f \n",layerMolMasses_[3]);
    }

    for (int i = 0; i < atom->nlocal; i++)
    {
        // calcMassLayer(i);
        layerRad(i);

        if (screen)
        {
            fprintf(screen, "porosity: %f \n", porosity_[i]);
            fprintf(screen, "tortuosity: %f \n", tortuosity_[i]);
            fprintf(screen, "pore diameter: %f \n", pore_diameter_[i]);
        }
 /*       if (screen)
        {
            fprintf(screen, "relRad_[0][0] = %f \n",relRadii_[i][0]);
            fprintf(screen, "relRad_[0][1] = %f \n",relRadii_[i][1]);
            fprintf(screen, "relRad_[0][2] = %f \n",relRadii_[i][2]);
            fprintf(screen, "relRad_[0][3] = %f \n",relRadii_[i][3]);
        }*/
    }

    if (comm -> me == 0 && screen)
        fprintf(screen, "execute init() \n");
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    updatePtrs();
    int i;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;

    double a_[nmaxlayers_];              // reaction resistance value for each layer
    double b_[nmaxlayers_];              // diffusion resistance value
    double x0_[nlocal];
    double x0_eq_[nmaxlayers_];          // molar fraction of reactant gas
    double dmA_[nmaxlayers_];            // mass flow rate of reactant gas species for each layer
    // double masst;

    for (i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (dCoeff_[i] > 1)
                dCoeff_[i] = 1;

            // debug check
/*            fprintf(screen, "T_ in chem core: %f \n", T_[i]);
            fprintf(screen ,"pgas_ in chem core = %f \n",pgas_[i]); */
            if (screen)
                fprintf(screen,"diffusion coefficient from CFD: %f \n",dCoeff_[i]);

            if(T_[i] > 0.0)
            {         
                active_layers(i);
                if (layers_ > 0)
                {
                    getXi(i,x0_,x0_eq_);
                    if (x0_[i] > 0.0)
                    {
                        for (int j = 0; j < nmaxlayers_; j++)
                        {
                            a_[j] = 0.0;
                            dmA_[j] = 0.0;
                        }
                        getA(i,a_, x0_eq_);
                        getB(i,b_);
                        // getMassT(i,masst);
                        reaction(i,a_,dmA_,x0_,x0_eq_, b_); //masst,
                        // reaction2(i,a_,dmA_,x0_,x0_eq_);
                        update_atom_properties(i,dmA_);
                        update_gas_properties(i,dmA_);
                    }
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

/* ------------------- compute layer radius with relRadii --------------- */

void FixChemShrinkCore::layerRad(int i)
{
   /* if (screen)
        fprintf(screen,"CALCULATE LAYER RADIUS \n"); */

    for (int j = 0 ; j <= nmaxlayers_; j++)
    {
        radLayer_[i][j] =   relRadii_[i][j]*radius_[i];
    }
}

/* ---------------------------------------------------------------------- */

// returns number of active layers
int FixChemShrinkCore::active_layers(int i)
{
    for(int j  = 1; j <= layers_; j++)
    {
        if (radLayer_[i][j] < rmin_)
        {
            layers_ -= 1;
            calcMassLayer(i);
            break;
        }
    }
    if (comm -> me == 0 && screen)
        fprintf(screen,"active layers: %i \n", layers_);
    return layers_;
}

/* ---------------------------------------------------------------------- */

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

void FixChemShrinkCore::getXi(int i, double *x0_, double *x0_eq_)
{
    /*if (screen)
        fprintf(screen,"DO GET X0 FUNCTION!!! \n"); */

    // calculate equilibrium concentration from K_eq(layer, T)
    double xp_eq_[nmaxlayers_];
    double xp_[atom->nlocal];

    for (int j = 0; j < layers_; j++)
    {
        xp_eq_[j]  =   K_eq(j,T_[i])/(1+K_eq(j,T_[i]));
        x0_eq_[j]  =   1-xp_eq_[j];

/*        // for debug
        // ==========================================
        if (screen)
        {
            fprintf(screen,"x0_eq : %f \n", x0_eq_[j]);
            fprintf(screen,"xp_eq : %f \n", xp_eq_[j]);
        }
        // ==========================================*/
    }

    // calculate bulk mole fraction
    // value calculated from N volScalarField
    if (N_[i] == 0.0)
    {
        x0_[i]  =   0.0;
        xp_[i]  =   0.0;
    }
    else
    {
       x0_[i]  =   concA_[i]*rhogas_[i]/(N_[i]*molMass_A_);
       if (x0_[i] > 1.0)
            x0_[i] = 1.0;
        xp_[i] = 1 - x0_[i];
    }

 /*   if (screen)
    {
        fprintf(screen,"check rhogas %f \n", rhogas_[i]);
        fprintf(screen,"check N %f \n", N_[i]);
        fprintf(screen,"molMass_A_ %f \n", molMass_A_);
        fprintf(screen,"check mass frac %f \n", concA_[i]);
    }*/

/*    // for debug
    // ==========================================
    if (screen)
    {
        fprintf(screen,"x0_: %f \n", x0_[i]);
        fprintf(screen,"xp_: %f \n", xp_[i]);
    }*/
    // ==========================================

    // directly guessing x_eq like in the thesis of Nietrost seems to be spurious for multi-component mixtures (H2, H2O, CO, CO2)
    // Negri et al. (1991), Takahashi et al. (1986) and probably many others use a different approach:
    // as driving term, they use (n_A - n_C / K)

    // ATTENTION: molar concentrations n_i [mol / m^3] need to be converted into mass fractions y_i
    //            y_i * rho = n_i * m_{i,mol}
}

/* ---------------------------------------------------------------------- */

// calculate chemical reaction resistance term through literature
void FixChemShrinkCore::getA(int i, double *a_, double *x0_eq_)
{
    /*if (screen)
    fprintf(screen,"DO GET A FUNCTION!! \n");*/

    for (int j = 0; j < layers_ ; j++)
    {
        a_[j]   =   (1/(k0_[j]*exp(-Ea_[j]/(Runiv*T_[i]))))*(1/(1+1/K_eq(j,T_[i])))*(1/(relRadii_[i][j+1]*relRadii_[i][j+1]));
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getB(int i, double *b_)
{
    /*if (screen)
    fprintf(screen,"DO GET B FUNCTION!! \n"); */

    double deBinary_[atom->nlocal];
    double deKnudsen_[atom->nlocal];

    if (dCoeff_[i] < SMALL)
        dCoeff_[i]  =   SMALL;

    // molecular binary diffusion
    deBinary_[i] = dCoeff_[i]*porosity_[i]/tortuosity_[i];

    if (screen)
        fprintf(screen,"eff. binary diff: %6.15f \n",deBinary_[i]);

    // Knudsen diff equation is either
    // Dik = dp/3*sqrt((8*R*T_[i])/(M_PI*molMass_A_))
    // or simplified as
    // Dik = 4850*dp*sqrt(T_[i]/molMass_A_) dp supposed to be in cm
    // we use si units so convert from meter to cm!
    deKnudsen_[i]    =   4850*(pore_diameter_[i]*0.01)*sqrt(T_[i]/molMass_A_)*porosity_[i]/tortuosity_[i];

    if (screen)
        fprintf(screen,"eff. knudsen diff: %f \n",deKnudsen_[i]);

    // total effective diffusivity (diffEff = 1/EffectiveDiffusivity)
    diffEff_[i]     =   1/deBinary_[i] + 1/deKnudsen_[i];
    // diffEff_[i]     =   deKnudsen_[i]*deBinary_[i]/(deKnudsen_[i]+deBinary_[i]);

    if (screen)
        fprintf(screen,"eff. diff: %6.15f \n",diffEff_[i]);

    // calculation of Diffustion Term from Valipour 2009
    // from wustite to iron
    b_[0]   =   (1-relRadii_[i][1])/(relRadii_[i][1])*(radius_[i]*diffEff_[i]);
    for (int j = 1; j < layers_; j++)
    {
        b_[j]   =   (relRadii_[i][j] - relRadii_[i][j+1])/(relRadii_[i][j]*relRadii_[i][j+1])*(radius_[i]*diffEff_[i]);
    }

    if (screen)
    {
        fprintf(screen,"diffusion constant [0]: %f \n",b_[0]);
        fprintf(screen,"diffusion constant [1]: %f \n",b_[1]);
        fprintf(screen,"diffusion constant [2]: %f \n",b_[2]);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction(int i, double *a_, double *dmA_, double *x0_, double *x0_eq_, double *b_) // double masst, double* y0,
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
        dY[2]   =   ((a_[0]*(a_[1]+b_[1]+b_[0])+b_[0]*(a_[1]+b_[1]))*(x0_[i]-x0_eq_[2])
                -   (a_[0]*(b_[1]+b_[0])+b_[1]*b_[0])*(x0_[i]-x0_eq_[1])
                -   (a_[1]*b_[0])*(x0_[i]-x0_eq_[0]))*1/W;
        // magnetite to wustite
        dY[1]   =   (((a_[2]+b_[2]+b_[1])*(a_[0]+b_[0])+a_[0]*b_[0])*(x0_[i]-x0_eq_[1])
                -   (b_[1]*(a_[0]+b_[0])+a_[0]*b_[0])*(x0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*b_[0])*(x0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   (((a_[2]+b_[2])*(a_[1]+b_[1]+b_[0])+a_[1]*(b_[1]+b_[0]))*(x0_[i]-x0_eq_[0])
                -   (a_[1]*b_[0])*(x0_[i]-x0_eq_[2])
                -   ((a_[2]+b_[2])*b_[0])*(x0_[i]-x0_eq_[1]))*1/W;
    }
     else if (layers_ == 2)
     {
        W = (a_[1]+b_[1])*(a_[0]+b_[0])+a_[0]*b_[0];

        // hematite to magnetite
        dY[2]   =   0.0;
        // magnetite to wustite
        dY[1]   =   ((a_[0]+b_[0])*(x0_[i]-x0_eq_[1])-b_[0]*(x0_[i]-x0_eq_[0]))*1/W;
        // wustite to iron
        dY[0]   =   ((a_[1]+b_[1]+b_[0])*(x0_[i]-x0_eq_[0])-b_[0]*(x0_[i]-x0_eq_[1]))*1/W;
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
        dY[0]   =   (x0_[i] - x0_eq_[0])*1/W;
     }

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        // mass flow rate for reactant gas species
        // scaled up Timestep to match 0.01;
        dmA_[j] =   dY[j]*rhogas_[i]*partSurfArea(radius_[i])*TimeStep;
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction2(int i, double *a_, double *dmA_, double *x0_, double *x0_eq_) // double masst, double* y0,
{
   /* updatePtrs();
    double W;
    double dY[nmaxlayers_];

    if (layers_ == nmaxlayers_)
    {
        // if only with reaction resistance term
        // rate of chemical reaction for only chemical reaction resistance.
        W = a_[2]*a_[1]*a_[0];
        // hematite to magnetite
        dY[2]   =   ((a_[0]*a_[1])*(x0_[i]-x0_eq_[2]))*1/W;
        // magnetite to wustite
        dY[1]   =   ((a_[2]*a_[0])*(x0_[i]-x0_eq_[1]))*1/W;
        // wustite to iron
        dY[0]   =   ((a_[2]*a_[1])*(x0_[i]-x0_eq_[0]))*1/W;
    }
     else if (layers_ == 2)
     {
         // rate of chemical reaction for 2 active layers
         W = (a_[1]*a_[0]);

         // hematite to magnetite
         dY[2]   =   0.0;
         // magnetite to wustite
         dY[1]   =   (a_[0]*(x0_[i]-x0_eq_[1]))*1/W;
         // wustite to iron
         dY[0]   =   (a_[1]*(x0_[i] - x0_eq_[0]))*1/W;
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
         dY[0]   =   (x0_[i] - x0_eq_[0])*1/W;
     }

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        dmA_[j] =   dY[j]*rhogas_[i]*partSurfArea(radius_[i])*TimeStep; //*10000;
    } */
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, double *dmA_)
{
    // based on material change: update relative radii, average density and mass of atom i        //  mass of each layer
    double dmB_[nmaxlayers_+1];               //  mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass_new = 0;                  // just for control case --- delete afterwards
    double vBr_[3] = {1, 1, 3};
    double vBp_[3] = {1, 3, 2};

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

/*    // print out mass transfer rates for particle layers
    if (screen)
    {
        fprintf(screen,"dmB[0]: %f \n",dmB_[0]);
        fprintf(screen,"dmB[1]: %f \n",dmB_[1]);
        fprintf(screen,"dmB[2]: %f \n",dmB_[2]);
        fprintf(screen,"dmB[3]: %f \n",dmB_[3]);
    } */
    
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
    
    radLayer_[i][layers_]   =   pow((0.75*massLayer_[i][layers_]/(layerDensities_[layers_]*M_PI)),0.333333);
    for (int j = layers_-1; j >= 0 ; j--)
    {
        radLayer_[i][j]   =   pow((0.75*massLayer_[i][j]/(M_PI*layerDensities_[j])+radLayer_[i][j+1]*radLayer_[i][j+1]*radLayer_[i][j+1]),0.333333);
    }
    radius_[i] = radLayer_[i][0];
    
    //detemine the new relative layer radii
    for (int j = 0; j <= layers_; j++)
    {
        relRadii_[i][j] = radLayer_[i][j]/radius_[i];
    }
    
    // calculate particle total density
    // depending on the new values
    // update total partilce density
    pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);

/*    if (screen)
    {
        fprintf(screen,"post-layerMass[0]: %f \n",massLayer_[i][0]);
        fprintf(screen,"post-layerMass[1]: %f \n",massLayer_[i][1]);
        fprintf(screen,"post-layerMass[2]: %f \n",massLayer_[i][2]);
        fprintf(screen,"post-layerMass[3]: %f \n",massLayer_[i][3]);
        fprintf(screen,"post-new mass: %f \n",pmass_[i]);

        fprintf(screen,"post redox radius of particle (rl): %f \n",radLayer_[i][0]);
        fprintf(screen,"post redox radius of wustite: %f \n",radLayer_[i][1]);
        fprintf(screen,"post redox radius of magnetite: %f \n",radLayer_[i][2]);
        fprintf(screen,"post redox radius of hematite: %f \n",radLayer_[i][3]);
        fprintf(screen,"post redox radius of particle (R): %f \n",radius_[i]);

        fprintf(screen,"post-relrad[0]: %f \n",relRadii_[i][0]);
        fprintf(screen,"post-relrad[1]: %f \n",relRadii_[i][1]);
        fprintf(screen,"post-relrad[2]: %f \n",relRadii_[i][2]);
        fprintf(screen,"post-relrad[3]: %f \n",relRadii_[i][3]);
        fprintf(screen,"post-particle density: %f \n", pdensity_[i]);
    } */
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
 /*       if (screen)
        {
            fprintf(screen, "Available CO for 0: %f \n", changeOfA_[i]);
            fprintf(screen, "Available CO for 1: %f \n", changeOfA_[i]);
            fprintf(screen, "Available CO for 2: %f \n", changeOfA_[i]);
        } */

        changeOfC_[i]   +=  dmA_[j]*molMass_C_/molMass_A_;
    }

/*    if (screen)
    {
        fprintf(screen,"changeOfA_: %f \n",changeOfA_[i]);
        fprintf(screen,"changeOfC_: %f \n",changeOfC_[i]);
    }   */
}

/* ---------------------------------------------------------------------- */

double FixChemShrinkCore::K_eq(int layer, double T)
{
    /*if (screen)
        fprintf(screen,"CALCULATE K_EQUILIBIRUM \n"); */
    // ATTENTION: K is usually given for partial pressures, might be necessary to convert to molar concentrations
    // for the reactions under consideration, this should not be the case (double-check !!!)
    // 0 = wustite, 1 = mangetite, 2 = hematite;
    double Keq_;
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
