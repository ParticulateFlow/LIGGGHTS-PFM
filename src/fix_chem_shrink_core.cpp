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

#define SMALL 0.0001

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    nmaxlayers_(3),
    rmin_(0.01),
    drmin_(0.001),
    fix_k0_(0), // [m/s]
    fix_molMass_(0),
    fix_Ea_(0), // [J/mol]
    fix_dens_(0),
    fix_layerRelRad_(0)
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

    pgas_ = 0;
    molMass_A_ = molMass_C_ = 0;
    Runiv       =   8.3144;  // [J/molK]
    radChange   =   false;

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

    /* // flags for vector output
    vector_flag = 1;  //0/1 per-atom data is stored
    size_vector = 2;
    global_freq = 1;
    extvector   = 1; */

    radLayer_ = new double*[atom->nlocal];
    for (int i = 0; i < atom->nlocal; i++)
    {
        radLayer_[i]    =   new double [nmaxlayers_+1]; // nmaxlayers_ + 1
    }

    nevery = 1;
    time_depend = 1;
    force_reneighbor = 1;
    next_reneighbor = update -> ntimestep + nevery;

    restart_global = 1;

    // if debug
    if (comm-> me == 0 && screen)
    {
        fprintf(screen,"modified species names are: %s \n", massA);
        fprintf(screen,"and: %s \n", massC);
        fprintf(screen,"reactant species: %s \n", speciesA);
        fprintf(screen,"product species: %s \n", speciesC);
    }
}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
    delete massA;
    delete massC;

    delete speciesA;
    delete speciesC;

    for (int i=0; i < atom -> nlocal; i++)
    {
        delete[] radLayer_[i];
    }
    delete [] radLayer_;
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

        // fixes from shrink/core property/globals will not be deleted
        if (fix_layerRelRad_) modify->  delete_fix("relRadii");
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
     if (!fix_layerRelRad_)
    {
        const char* fixarg[12];
        fixarg[0]="layerRelRad";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="relRadii";
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.3";
        fixarg[9]="0.5";
        fixarg[10]="0.6";
        fix_layerRelRad_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
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

    // material prop
    // radius handle
    relRadii_       = fix_layerRelRad_-> array_atom;
    // density handle
    layerDensities_ = fix_dens_      -> get_values();
    // molar mass handle
    layerMolMasses_ = fix_molMass_   -> get_values();

    // chemical prop
    // frequency factor
    k0_             =  fix_k0_       -> get_values();
    // activation energy
    Ea_             =  fix_Ea_       -> get_values();
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
    char *group_name_ = new char [strlen(group->names[igroup])];
    strcpy(group_name_, group->names[igroup]);
    fixname = new char [c];
    strcpy(fixname,"density_");
    strcat(fixname,group_name_);
    if (screen)
        fprintf(screen, "density_ = %s \n", fixname);
    fix_dens_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname;


    /*  group_name_ is not defined again since the group is the same one as the density
        if need be uncomment following;
        delete []group_name_;
        char *group_name_ = new char [strlen(group->names[igroup])];
        strcpy(group_name_, group->names[igroup]); */
    fixname = new char[c];
    strcpy(fixname,"molMass_");
    strcat(fixname,group_name_);
    if (screen)
        fprintf(screen, "molMass_ = %s \n", fixname);
    fix_molMass_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",ntype,0,"FixChemShrinkCore"));
    delete []fixname;
    delete []group_name_;

    // look up *fix_layerRelRad_;
    fix_layerRelRad_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("relRadii","property/atom","vector",ntype,0,"FixChemShrinkCore"));

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

    updatePtrs();

    dens_ = new double [nmaxlayers_+1];
    // change from const double to double to modify and write layerDensities later on
    for (int i = 0; i <= nmaxlayers_; i++)
    {
        dens_[i] = layerDensities_[i];
    }

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
        fprintf(screen, "dens_[0] = %f \n",dens_[0]);
        fprintf(screen, "dens_[1] = %f \n",dens_[1]);
        fprintf(screen, "dens_[2] = %f \n",dens_[2]);
        fprintf(screen, "dens_[3] = %f \n",dens_[3]);
        fprintf(screen, "molMass_[0] = %f \n",layerMolMasses_[0]);
        fprintf(screen, "molMass_[1] = %f \n",layerMolMasses_[1]);
        fprintf(screen, "molMass_[2] = %f \n",layerMolMasses_[2]);
        fprintf(screen, "molMass_[3] = %f \n",layerMolMasses_[3]);
    }

    for (int i = 0; i < atom->nlocal;i++)
    {
        if (screen)
        {
            fprintf(screen, "relRad_[0][0] = %f \n",relRadii_[i][0]);
            fprintf(screen, "relRad_[0][1] = %f \n",relRadii_[i][1]);
            fprintf(screen, "relRad_[0][2] = %f \n",relRadii_[i][2]);
            fprintf(screen, "T_ in chem core: %f \n", T_[i]);
            fprintf(screen ,"pgas_ in chem core = %f \n",pgas_[i]);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    updatePtrs();
    radius_     = atom ->  radius;
    pmass_      = atom ->  rmass;
    pdensity_   = atom -> density;
    TimeStep    = update -> dt;
    current_timestep = update->ntimestep;

    int i;
    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    Q_ = new double [nlocal];

    if (screen)
    {
        fprintf(screen, "check current Timestep = %li \n", current_timestep);
    }

    double a_[nmaxlayers_];                 // reaction resistance value for each layer
    double b_[nmaxlayers_];              // diffusion resistance value
    double x0_[nmaxlayers_];
    double x0_eq_[nmaxlayers_];             // molar fraction of reactant gas
    double dmA_[nmaxlayers_];               // mass flow rate of reactant gas species for each layer
    double diff[nmaxlayers_];
    // double masst;

    for (i = 0; i < nlocal; i++)
    {
        if (screen)
            fprintf(screen, "current Temperature = %f \n", T_[i]);
        if (mask[i] & groupbit)
        {
            if(T_[i] > 0.0)
            {
                getXi(i,x0_,x0_eq_);
                if (x0_[i] > 0.0)
                {
                    // if (!radChange)
                        layerRad(i);
                    getA(i,a_, x0_eq_);
                    // diffcoeff(i,diff);
                    // getB(i,b_,diff);
                    // getDiff(i,diff);
                    // getMassT(i,masst);
                    for(int j = 0; j<nmaxlayers_; j++)
                        dmA_[j] = 0.0;
                    reaction(i,a_,dmA_,x0_,x0_eq_); //diff,masst,y0,
                    update_atom_properties(i,dmA_);
                    update_gas_properties(i,dmA_);
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
int FixChemShrinkCore::active_layers(int i)
{
    /* int layers;
    for(layers=0; layers<nmaxlayers_; layers++)
        if(relRadii_[i][layers] < rmin_)
            break;
    return layers; */
    return 0;
}

/* ------------------- compute layer radius with relRadii --------------- */

void FixChemShrinkCore::layerRad(int i)
{
    if (screen)
        fprintf(screen,"CALCULATE LAYER RADIUS \n");

    for (int j = 0; j <= nmaxlayers_; j++)
    {
        if ( j == nmaxlayers_)
            radLayer_[i][j] =   radius_[i];
        else
            radLayer_[i][j] =   relRadii_[i][j]*radius_[i];
    }

    // print radiuses of layers
    if (screen)
    {
        fprintf(screen,"radius of hematite: %f \n", radLayer_[i][0]);
        fprintf(screen,"radius of magnetite: %f \n", radLayer_[i][1]);
        fprintf(screen,"radius of wustite: %f \n", radLayer_[i][2]);
        fprintf(screen,"particle radius: %f \n", radLayer_[i][3]);
    }
}

/* ---------------------------------------------------------------------- */

// calculate chemical reaction resistance term through literature
void FixChemShrinkCore::getA(int i, double *a_, double *x0_eq_)
{
    if (screen)
    fprintf(screen,"DO GET A FUNCTION!! \n");

    for (int j = 0; j < nmaxlayers_; j++)
    {
        if (radLayer_[i][j] > rmin_)
            a_[j]   =   1/(k0_[j]*exp(-Ea_[j]/(Runiv*T_[i])))*(1-x0_eq_[i])*(1/(relRadii_[i][j]*relRadii_[i][j]));
    }

    if (screen)
    {
        fprintf(screen,"reaction constant [0]: %f \n",a_[0]);
        fprintf(screen,"reaction constant [1]: %f \n",a_[1]);
        fprintf(screen,"reaction constant [2]: %f \n",a_[2]);
        fprintf(screen,"equilibirum constant for hematite %f \n", K_eq(0,T_[i]));
        fprintf(screen,"equilibirum constant for magnetite %f \n", K_eq(1,T_[i]));
        fprintf(screen,"equilibirum constant for wustite %f \n", K_eq(2,T_[i]));
    }
    // min amount of layer thickness must be added
    // if layer J thickness < drmin or layer J radius < rmin, a[J] = LARGE
    // if T < Tcrit1, a = LARGE
}

/* ---------------------------------------------------------------------- */

 // calculate different correlations for effective diffusion
 // such as Fuller-Schetler-Giddings / Chapman-Enskog
 // current diffusion values taken from Natsui et al. The Minerals, Metals & Materials Society and ASM International 2014
void FixChemShrinkCore::diffcoeff(int i, double *diff)
{
    if (screen)
        fprintf(screen,"DO DIFF COEFF CALC \n");

    diff[0] =   exp(4.135-14080/T_[i])/pgas_[i];
    diff[1] =   exp(-1.835-7180/T_[i])/pgas_[i];
    diff[2] =   exp(0.485-8770/T_[i])/pgas_[i];

    if (screen)
    {
        fprintf(screen,"diff coeff [0]: %f \n",diff[0]);
        fprintf(screen,"diff coeff [1]: %f \n",diff[1]);
        fprintf(screen,"diff coeff [2] : %f \n",diff[2]);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getB(int i, double *b_, double *diff)
{
    if (screen)
    fprintf(screen,"DO GET B FUNCTION!! \n");
    double f_[nmaxlayers_];

    for (int j = 0; j < nmaxlayers_ ; j++)
    {
        f_[j]   =   1-relRadii_[i][j]*relRadii_[i][j]*relRadii_[i][j];
    }

    b_[0]   =   ((pow((1-f_[1]),0.3333333)-pow((1-f_[0]),0.3333333))*radius_[i])
                /((pow((1-f_[0]),0.3333333)*pow((1-f_[1]),0.3333333)*diff[0]));

    b_[1]   =   ((pow((1-f_[2]),0.3333333)-pow((1-f_[1]),0.3333333))*radius_[i])
                /((pow((1-f_[1]),0.3333333)*pow((1-f_[2]),0.3333333)*diff[1]));

    b_[2]   =   ((1-pow((1-f_[2]),0.3333333))*radius_[i])
                /((pow((1-f_[2]),0.3333333)*diff[2]));

    fprintf(screen,"diffusion constant [0]: %f \n",b_[0]);
    fprintf(screen,"diffusion constant [1]: %f \n",b_[1]);
    fprintf(screen,"diffusion constant [2]: %f \n",b_[2]);
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::getXi(int i, double *x0_, double *x0_eq_)
{
    if (screen)
        fprintf(screen,"DO GET X0 FUNCTION!!! \n");

    // calculate equilibrium concentration from K_eq(layer, T)
    double xp_eq_[nmaxlayers_];
    double xp_[atom->nlocal];

    for (int j = 0; j < nmaxlayers_; j++)
    {
        xp_eq_[j]  =   K_eq(j,T_[i])/(1+K_eq(j,T_[i]));
        x0_eq_[j]  =   1-xp_eq_[j];

        // for debug
        // ==========================================
        if (screen)
        {
            fprintf(screen,"x0_eq : %f \n", x0_eq_[j]);
            fprintf(screen,"xp_eq : %f \n", xp_eq_[j]);
        }
        // ==========================================
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

    // calculate reaction quotient
    if (x0_[i] > 0.0)
        Q_[i] = xp_[i]/x0_[i];
    else
        Q_[i] = 0.0;

    if (screen)
    {
        fprintf(screen,"check rhogas %f \n", rhogas_[i]);
        fprintf(screen,"check N %f \n", N_[i]);
        fprintf(screen,"molMass_A_ %f \n", molMass_A_);
        fprintf(screen,"check mass frac %f \n", concA_[i]);
    }

    // for debug
    // ==========================================
    if (screen)
    {
        fprintf(screen,"reaction quotient: %f \n", Q_[i]);
        fprintf(screen,"x0_: %f \n", x0_[i]);
        fprintf(screen,"xp_: %f \n", xp_[i]);
    }
    // ==========================================

    // directly guessing x_eq like in the thesis of Nietrost seems to be spurious for multi-component mixtures (H2, H2O, CO, CO2)
    // Negri et al. (1991), Takahashi et al. (1986) and probably many others use a different approach:
    // as driving term, they use (n_A - n_C / K)

    // ATTENTION: molar concentrations n_i [mol / m^3] need to be converted into mass fractions y_i
    //            y_i * rho = n_i * m_{i,mol}
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::reaction(int i, double *a_, double *dmA_, double *x0_, double *x0_eq_) //double *diff, double masst, double* y0,
{
    updatePtrs();
    double W;
    double dY[nmaxlayers_];
    double vA_[nmaxlayers_];

    // rate of chemical reaction for only chemical reaction resistance.
    W = a_[0]*(a_[2]*a_[1]+a_[1])+a_[1]*a_[2];
    dY[0]   =   ((a_[2]*a_[1]+a_[1])*(x0_[i]-x0_eq_[0])-a_[2]*(x0_[i]-x0_eq_[1])-a_[1]*(x0_[i]-x0_eq_[2]))*1/W;
    dY[1]   =   ((a_[0]*a_[2]+a_[2])*(x0_[i]-x0_eq_[1])-(a_[2]+a_[2])*(x0_[i]-x0_eq_[0])-a_[0]*(x0_[i]-x0_eq_[2]))*1/W;
    dY[2]   =   ((a_[0]*a_[1]+a_[1])*(x0_[i]-x0_eq_[2])-a_[1]*(x0_[i]-x0_eq_[0])-a_[0]*(x0_[i]-x0_eq_[1]))*1/W;

    for (int j = 0 ; j < nmaxlayers_; j++)
    {
        // determine direction of reaction
        // this will change the stochiometric coeff
        // only the sign of stoc. coeff v will change depending on direction
        if (Q_[i] < K_eq(j,T_[i]))      // reaction favors products -- forward
            vA_[j] =  1;
        else if (Q_[i] > K_eq(j,T_[i])) // reaction favors reactants -- backward
            vA_[j] = -1;

        // mass flow rate for reactant gas species
        // value for dmA is always positive ?!
        dmA_[j] =   vA_[j]*dY[j]*rhogas_[i]*partSurfArea(radius_[i])*0.01;//Timestep;
    }

    if (screen)
    {
        fprintf(screen,"W = %f \n", W);
        fprintf(screen,"dY_[0]: %lf \n",dY[0]);
        fprintf(screen,"dY_[1]: %lf \n",dY[1]);
        fprintf(screen,"dY_[2]: %lf \n",dY[2]);

        fprintf(screen,"dmA_[0]: %lf \n",dmA_[0]);
        fprintf(screen,"dmA_[1]: %lf \n",dmA_[1]);
        fprintf(screen,"dmA_[2]: %lf \n",dmA_[2]);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, double *dmA_)
{
    // based on material change: update relative radii, average density and mass of atom i
    double layerMass_[nmaxlayers_+1];         //  mass of each layer
    double dmB_[nmaxlayers_+1];               //  mass flow rate between each layer i.e. (btw h->m, m->w, w->Fe) must consider reduction and growth at the same time
    double sum_mass = 0;
    double sum_mass_new = 0;                // just for control case --- delete afterwards
    double vBr_[nmaxlayers_] = {3, 1, 1};
    double vBp_[nmaxlayers_] = {2, 3, 1};

    // double diff_radius;
    // Calculate mass flow rates (dmB[j])
    // Keep in mind the stoichiometric coefficients of reaction prod - reactants
    dmB_[0] =   dmA_[0]*vBr_[0]*(layerMolMasses_[0]/molMass_A_);
    dmB_[1] =   dmA_[0]*(layerMolMasses_[1]/molMass_A_)*vBp_[0] - dmA_[1]*(layerMolMasses_[1]/molMass_A_);
    dmB_[2] =   dmA_[1]*(layerMolMasses_[2]/molMass_A_)*vBp_[1] - dmA_[2]*(layerMolMasses_[2]/molMass_A_);
    dmB_[3] =   dmA_[2]*(layerMolMasses_[3]/molMass_A_); // dmA is positive to increase iron layer negated

    // print out mass transfer rates for particle layers
    if (screen)
    {
        fprintf(screen,"dmB[0]: %f \n",dmB_[0]);
        fprintf(screen,"dmB[1]: %f \n",dmB_[1]);
        fprintf(screen,"dmB[2]: %f \n",dmB_[2]);
        fprintf(screen,"dmB[3]: %f \n",dmB_[3]);
    }

    //if (!radChange)
    //{
        // calculate inital layer masses
        // this should be skipped after first run (if(!radChange))
        for (int j=0; j <= nmaxlayers_; j++)
        {
            if (radLayer_[i][j] > rmin_)
            {
                if (j == 0)
                {
                    layerMass_[j]   =   1.33333*M_PI*radLayer_[i][j]*radLayer_[i][j]*radLayer_[i][j]*dens_[j];
                } else
                {
                    layerMass_[j]   =   1.33333*M_PI*((radLayer_[i][j]*radLayer_[i][j]*radLayer_[i][j])-(radLayer_[i][j-1]*radLayer_[i][j-1]*radLayer_[i][j-1]))*dens_[j];
                }

                // calculatue before reduction pmass and pdensity
                sum_mass        +=  layerMass_[j];
                pmass_[i]       =   sum_mass;
                pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);
            }
        }
    //}

    // print the initial mass out
    if (screen)
    {
        fprintf(screen,"pre-layerMass[0]: %f \n",layerMass_[0]);
        fprintf(screen,"pre-layerMass[1]: %f \n",layerMass_[1]);
        fprintf(screen,"pre-layerMass[2]: %f \n",layerMass_[2]);
        fprintf(screen,"pre-layerMass[3]: %f \n",layerMass_[3]);
        fprintf(screen,"dens0_ = %f \n",dens_[0]);
        fprintf(screen,"dens1_ = %f \n",dens_[1]);
        fprintf(screen,"dens2_ = %f \n",dens_[2]);
        fprintf(screen,"dens3_ = %f \n",dens_[3]);
        fprintf(screen,"pre-particle density: %f \n", pdensity_[i]);
        fprintf(screen,"pre-particle mass: %f \n", pmass_[i]);
    }

    for (int j = 0; j <= nmaxlayers_; j++)
    {
        if (radLayer_[i][j] > rmin_)
        {
            // WARNING 1: do not remove more material of layer J than present --> decrease dm[J] if necessary
            // calculate new masses of each layer
            if (dmB_[j] > layerMass_[j]) //|| -dmB_[j] > layerMass_[j])
                dmB_[j] = SMALL;

            if (layerMass_[j] > 0.0)  // add constraint so that diff layers > drmin
                layerMass_[j] -= dmB_[j];

            // calculatue after reduction pmass and pdensity
            sum_mass_new    +=  layerMass_[j];
            pmass_[i]   =   sum_mass_new;
        }
    }

    if (screen)
    {
        fprintf(screen,"post-layerMass[0]: %f \n",layerMass_[0]);
        fprintf(screen,"post-layerMass[1]: %f \n",layerMass_[1]);
        fprintf(screen,"post-layerMass[2]: %f \n",layerMass_[2]);
        fprintf(screen,"post-layerMass[3]: %f \n",layerMass_[3]);
        fprintf(screen,"post-new mass: %f \n",pmass_[i]);
    }

    // second step: based on predefined molar densities, compute new layer relRadii -- which is used to get new radii
    for (int j = 0; j <= nmaxlayers_; j++)
    {
        if (radLayer_[i][j] > rmin_ && layerMass_[j] > 0.0)
        {
            if (j == 0)
            {
                radLayer_[i][j]   =   pow((0.75*layerMass_[j]/(dens_[j]*M_PI)),0.333333);
            } else
            {
                radLayer_[i][j]   =   pow((0.75*layerMass_[j]/(M_PI*dens_[j])+radLayer_[i][j-1]*radLayer_[i][j-1]*radLayer_[i][j-1]),0.333333);
            }

            if (j == nmaxlayers_)
                radius_[i] = radLayer_[i][j];
        }
    }

    if (screen)
    {
        fprintf(screen,"post redox radius of hematite: %f \n",radLayer_[i][0]);
        fprintf(screen,"post redox radius of magnetite: %f \n",radLayer_[i][1]);
        fprintf(screen,"post redox radius of wustite: %f \n",radLayer_[i][2]);
        fprintf(screen,"post redox radius of particle (rl): %f \n",radLayer_[i][3]);
        fprintf(screen,"post redox radius of particle (R): %f \n",radius_[i]);
    }

    for (int j = 0; j <= nmaxlayers_; j++)
    {
        if (layerMass_[j] > 0.0)
        {
            if (j == 0)
                dens_[j]    =   0.75*layerMass_[j]/(M_PI*radLayer_[i][j]*radLayer_[i][j]*radLayer_[i][j]);
            else
                dens_[j]    =   0.75*layerMass_[j]/(M_PI*(radLayer_[i][j]*radLayer_[i][j]*radLayer_[i][j]-radLayer_[i][j-1]*radLayer_[i][j-1]*radLayer_[i][j-1]));
        }

        // calculate particle total density and mass
        // depending on the new values
        // update total partilce density
        pdensity_[i]    =   0.75*pmass_[i]/(M_PI*radius_[i]*radius_[i]*radius_[i]);
    }

    for (int j = 0; j < nmaxlayers_; j++)
    {
        if (radLayer_[i][j] > rmin_)
        {
            relRadii_[i][j] = radLayer_[i][j]/radius_[i];
        }
    }

    if (screen)
    {
        fprintf(screen,"post-relrad[0]: %f \n",relRadii_[i][0]);
        fprintf(screen,"post-relrad[1]: %f \n",relRadii_[i][1]);
        fprintf(screen,"post-relrad[2]: %f \n",relRadii_[i][2]);
        fprintf(screen,"post-particle density: %f \n", pdensity_[i]);
    }
    radChange = true;
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

    if (screen)
    {
        fprintf(screen,"changeOfA_: %f \n",changeOfA_[i]);
        fprintf(screen,"changeOfC_: %f \n",changeOfC_[i]);
    }
}

//void FixChemShrinkCore::reaction_1(int i, double *a, double *diff, double masst, double* y0, double *dm)
//{
  // obtained from reaction_3 by taking the limits a[1] -> infinity and a[2] -> infinity
  // maybe not needed
//}

//void FixChemShrinkCore::reaction_2(int i, double *a, double *diff, double masst, double* y0, double *dm)
//{
  // obtained from reaction_3 by taking the limit a[2] -> infinity
  // maybe not needed
//}



double FixChemShrinkCore::K_eq(int layer, double T)
{
    if (screen)
        fprintf(screen,"CALCULATE K_EQUILIBIRUM \n");
//     ATTENTION: K is usually given for partial pressures, might be necessary to convert to molar concentrations
//                for the reactions under consideration, this should not be the case (double-check !!!)
    double Keq_;
    if(strcmp(speciesA,"CO")==0)
     {
         if (layer == 0)
             Keq_   =   exp(3968.37/T+3.94);
         else if (layer == 1)
             Keq_   =   pow(10,(-1834/T+2.17));
         else if (layer == 2)
             Keq_   =   pow(10,(917/T-1.097));
     }
     else if(strcmp(speciesA,"H2")==0)
     {
         if (layer == 0)
             Keq_   =   exp(-362.6/T+10.344);
         else if (layer == 1)
             Keq_   =   pow(10,(-3577/T+3.74));
         else if (layer == 2)
             Keq_   =   pow(10,(-827/T+0.468));
     }
    else
     {
         printf("Error : Undefined Reaction \n");
     }
    return Keq_;
}
