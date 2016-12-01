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


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg),
	nmaxlayers_(3),
	rmin_(0.01),
        drmin_(0.001),
        fix_k0_(0),
        fix_molMass_(0),
        fix_Ea_(0),
        fix_dens_(0),
        fix_layerRelRad_(0)
{
        if (strncmp(style, "chem/shrink/core", 11) == 0 && (!atom->radius_flag) || (!atom->rmass_flag))
		error->all(FLERR, "Fix chem/shrink needs particle radius and mass");

	// defaults
	fix_concA_ = NULL;
	fix_concC_ = NULL;
	fix_changeOfA_ = NULL;
	fix_changeOfC_ = NULL;
        fix_rhogas_ = NULL;
        // rmin_ = NULL;
        fix_tgas_       =   NULL;
	// fix_reactionheat_    =   0;
        Runiv = 0.008314458;
        T  = 1073.15;

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

	// changeOfA and changeOfC have to have the same name as in the species.C 
	// therefore new strings of massA and massC are introduced which are defined as 
        // changeOfA and changeOfC as their names in species.c
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

	// flags for vector output
        vector_flag = 1;  //0/1 per-atom data is stored
        size_vector = 2;
        global_freq = 1;
        extvector   = 1;

}	

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
	// delete the strings that are constructed with using "new"	
	delete massA;
	delete massC;
	delete speciesA;
	delete speciesC;

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if (fix_concA_)     modify  ->  delete_fix(speciesA);
        if (fix_concC_)     modify  ->  delete_fix(speciesC);
        if (fix_rhogas_)    modify  ->  delete_fix("partRho");
        if(fix_tgas_)       modify  ->  delete_fix("partTemp");
        // if(fix_reactionHeat_)  modify  ->  delete_fix("reactionHeat");

        if (fix_changeOfA_) modify  ->  delete_fix(massA);
        if (fix_changeOfC_) modify  ->  delete_fix(massC);

        // fixes from shrink/core property/globals will not be deleted
        if (fix_layerRelRad_) modify->  delete_fix("layerradii");
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
        fixarg[3]="layerradii";
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="1.";
        fixarg[9]="1.";
        fixarg[10]="1.";
        fix_layerRelRad_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */
void FixChemShrinkCore::updatePtrs()
{
    changeOfA_      =   fix_changeOfA_  -> vector_atom;
    changeOfC_      =   fix_changeOfC_  -> vector_atom;
    rhogas_         =   fix_rhogas_     -> vector_atom;
    concA_          =   fix_concA_      -> vector_atom;
    concC_          =   fix_concC_      -> vector_atom;
    T_              =   fix_tgas_       -> vector_atom;
    
    relRadii_       = fix_layerRelRad_-> array_atom;
    layerDensities_ = fix_dens_      -> get_values();
    layerMolMasses_ = fix_molMass_   -> get_values();
    
    k0_             =  fix_k0_       -> get_values();
    Ea_             =  fix_Ea_       -> get_values();
    // other chemical parameters

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::init()
{
    if (!atom -> radius_flag || !atom -> density_flag)
        error -> all(FLERR,"Fix chem/shrink/core can only be used with sphere atom style");
    if (!atom->tag_enable || 0 == atom->map_style)
        error->fix_error(FLERR, this, "requires atom tags and an atom map");

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
    fix_layerRelRad_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("layerradii","property/atom","vector",ntype,0,"FixChemShrinkCore"));
    
    // references
    fix_concA_     = static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesA, "property/atom", "scalar", 0, 0, id));
    fix_concC_     = static_cast<FixPropertyAtom*>(modify->find_fix_property(speciesC, "property/atom", "scalar", 0, 0, id));
    fix_changeOfA_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(massA, "property/atom", "scalar", 0, 0, id));
    fix_changeOfC_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(massC, "property/atom", "scalar", 0, 0, id));
    fix_rhogas_    = static_cast<FixPropertyAtom*>(modify->find_fix_property("partRho", "property/atom", "scalar", 0, 0, id));
    fix_tgas_      = static_cast<FixPropertyAtom*>(modify->find_fix_property("partTemp", "property/atom", "scalar", 0, 0, id));
    //fix_reactionheat_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("reactionHeat", "property/atom", "scalar", 0, 0, id));

    updatePtrs();
    if (screen)
    {
        fprintf(screen, "k0_[0] = %f \n", k0_[0]);
        fprintf(screen, "k0_[1] = %f \n", k0_[1]);
        fprintf(screen, "k0_[2] = %f \n", k0_[2]);
        fprintf(screen, "Ea_[0] = %f \n", Ea_[0]);
        fprintf(screen, "Ea_[1] = %f \n", Ea_[1]);
        fprintf(screen, "Ea_[2] = %f \n", Ea_[2]);
        fprintf(screen, "dens_[0] = %f \n",layerDensities_[0]);
        fprintf(screen, "dens_[1] = %f \n",layerDensities_[1]);
        fprintf(screen, "dens_[2] = %f \n",layerDensities_[2]);
        fprintf(screen, "dens_[3] = %f \n",layerDensities_[3]);
        fprintf(screen, "molMass_[0] = %f \n",layerMolMasses_[0]);
        fprintf(screen, "molMass_[1] = %f \n",layerMolMasses_[1]);
        fprintf(screen, "molMass_[2] = %f \n",layerMolMasses_[2]);
        fprintf(screen, "molMass_[3] = %f \n",layerMolMasses_[3]);
        fprintf(screen, "relRad_[0][0] = %f \n",relRadii_[0][0]);
        fprintf(screen, "relRad_[0][1] = %f \n",relRadii_[0][1]);
        fprintf(screen, "relRad_[0][2] = %f \n",relRadii_[0][2]);
        fprintf(screen, "relRad_[0][3] = %f \n",relRadii_[0][3]);
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

    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;

    double a_[nmaxlayers_];
    double dm[nmaxlayers_+1];

    /* double diff[nmaxlayers_];
    double masst;
    double y0[nmaxlayers_]; */

    for (int i = 0; i<nlocal; i++)
    {
        if (mask[i] & groupbit)
	{
            getA(i,a_);
//	    getDiff(i,diff);
//	    getMassT(i,masst);
//	    getY0(i,y0);
	    
           /* for(int j = 0; j<nmaxlayers_; j++)
                dm[j] = 0.0;
            reaction(i,a,diff,masst,y0,dm); */
           update_atom_properties(i,dm);
         /*   update_gas_properties(i,dm); */
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
    int layers;

    for(layers=0; layers<nmaxlayers_; layers++)
        if(relRadii_[i][layers] < rmin_)
            break;
    return layers;
}

/* ---------------------------------------------------------------------- */

// calculated with Valipour's equation
void FixChemShrinkCore::getA(int i, double *a_)

{
    for (int j = 0; j < nmaxlayers_; j++)
    {
        a_[j] = 1/(k0_[j]*exp(-Ea_[j]/(Runiv*T))*((relRadii_[i][j]*relRadii_[i][j])/(radius_[i]*radius_[i])));

        if (screen)
        {
            fprintf(screen,"reaction constant: %f \n",a_[j]);
        }
    }
    // if layer J thickness < drmin or layer J radius < rmin, a[J] = LARGE
    // if T < Tcrit1, a = LARGE
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_atom_properties(int i, double *dm)
{
    // based on material change: update relative radii, average density and mass of atom i

    double A[nmaxlayers_];
    double W;
    double dA[nmaxlayers_];
    double layerMass_[nmaxlayers_+1];

    getA(i,A);
    if (screen)
    {
        fprintf(screen,"other_function a_: %f \n",A[0]);
        fprintf(screen,"other_function a_: %f \n",A[1]);
        fprintf(screen,"other_function a_: %f \n",A[2]);
    }

    // determine concentration flow rate -- use the network resistance calculations
    // concentration flow rate between hematite - magnetite layer
    // ref. Octave Levenspiel pg. 595, IronOre presentation Nietrost.

    W = (A[0]*A[1]+A[1])+A[1]*A[2];
    dA[0] = ((A[2]*A[1]+A[1])*concA_[i]-(A[2])*concA_[i]-(A[1])*concA_[i])*(1/W);//*(partSurfArea(radius_[i])*rhogas_[i]*TimeStep);
    dA[1] = ((A[0]*A[2]+A[2])*concA_[i]-(A[2]+A[2])*concA_[i]-A[0]*concA_[i])*(1/W);//*(partSurfArea(radius_[i])*rhogas_[i]*TimeStep);
    dA[2] = ((A[0]*A[1]+A[1])*concA_[i]-A[0]*concA_[i]-A[1]*concA_[i])*(1/W);//*(partSurfArea(radius_[i])*rhogas_[i]*TimeStep);

    if (screen)
    {
        fprintf(screen,"dA[0]: %f \n",dA[0]);
        fprintf(screen,"dA[1]: %f \n",dA[1]);
        fprintf(screen,"dA[2]: %f \n",dA[2]);
    }

    /*     // mass flow rates from hematite -> magnetite
        dm[0]           =   dA[0]*(layerMolMasses_[0]/molMass_A_);  //  negative sign since this layer is shrinking
        dm[1]           =   dA[0]*(layerMolMasses_[1]/molMass_A_);  //  positive - this layer is gaining

    // mass flow rates from magentite -> wustite
        dm[1]           =   dA[1]*(layerMolMasses_[1]/molMass_A_);  //  negative
        dm[2]           =   dA[1]*(layerMolMasses_[2]/molMass_A_);  //  positive

    // mass flow rates from wustite -> iron
        dm[2]           =   dA[2]*(layerMolMasses_[2]/molMass_A_);  //  negative
        dm[3]           =   dA[2]*(layerMolMasses_[3]/molMass_A_);  //  positive */

    // mass flow rates for every layer
    dm[0]           =   -dA[0]*(layerMolMasses_[0]/molMass_A_);
    dm[1]           =   dA[0]*(layerMolMasses_[1]/molMass_A_) - dA[1]*(layerMolMasses_[1]/molMass_A_);
    dm[2]           =   dA[1]*(layerMolMasses_[2]/molMass_A_) - dA[2]*(layerMolMasses_[2]/molMass_A_);
    dm[3]           =   dA[2]*(layerMolMasses_[3]/molMass_A_);

    if (screen)
    {
        fprintf(screen,"dm[0]: %f \n",dm[0]);
        fprintf(screen,"dm[1]: %f \n",dm[1]);
        fprintf(screen,"dm[2]: %f \n",dm[2]);
        fprintf(screen,"dm[3]: %f \n",dm[3]);
    }

    // calculate new masses of each layer
    for (int j = 0 ; j < nmaxlayers_+1 ; j++)
    {
        // calculate mass of layers for initial radius
        layerMass_[j] = 1.3333*M_PI*relRadii_[i][j]*relRadii_[i][j]*relRadii_[i][j]*layerDensities_[j];

        // first step: get new masses of each layer
        // mass of layers
        layerMass_[j] += dm[j];

        // second step: based on predefined molar densities, compute new layer radii
        // new layer radii
        // relRadii_[i][j] = pow((0.75*layerMass_[j]/(M_PI*layerDensities_[j])),0.33333);
        // relRadii_[i][j] = 0.02;
    }

    if (screen)
    {
        fprintf(screen,"layerMass: %f \n",layerMass_[0]);
        fprintf(screen,"layerMass: %f \n",layerMass_[1]);
        fprintf(screen,"layerMass: %f \n",layerMass_[2]);
        fprintf(screen,"layerMass: %f \n",layerMass_[3]);
    }

 // WARNING 1: do not remove more material of layer J than present --> decrease dm[J] if necessary
  
 // WARNING 2: if Wustite layer thickness < rmin AND T < 567 C, direct conversion of Magnetite to Fe is possible
 //            this means that Wu-Fe radius is shrinking together with the Ma-Wu radius 
 //            in this case, the equilibrium value and reaction parameters at the Ma layer need to be adapted within
 //            IMPLEMENTATION ONLY IF NECESSARY
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::update_gas_properties(int i, double *dm)
{
   // based on material change: update gas-phase source terms for mass and heat
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

void FixChemShrinkCore::reaction(int i, double *a, double *diff, double masst, double* y0, double *dm)
{
//   dm[0] = ... ;
//   dm[1] = ... ;
//   dm[2] = ... ;
//    
//   dm[0] *= update->dt * nevery;
//   dm[1] *= update->dt * nevery;
//   dm[2] *= update->dt * nevery;
}

void FixChemShrinkCore::getY0(int i, double *y0)
{
  // calculate equilibrium concentration from K_eq(layer, T)
  
  // directly guessing x_eq like in the thesis of Nietrost seems to be spurious for multi-component mixtures (H2, H2O, CO, CO2)
  // Negri et al. (1991), Takahashi et al. (1986) and probably many others use a different approach:
  // as driving term, they use (n_A - n_C / K)
  
  // ATTENTION: molar concentrations n_i [mol / m^3] need to be converted into mass fractions y_i
  //            y_i * rho = n_i * m_{i,mol}

}

double FixChemShrinkCore::K_eq(int layer, double T)
{
//     ATTENTION: K is usually given for partial pressures, might be necessary to convert to molar concentrations
//                for the reactions under consideration, this should not be the case (double-check !!!)  
  
//     if(speciesA == 'CO')
//     {
//         K_eq depending on layer and T       
//     }
//     else if(speciesA == 'H2')
//     {
//         K_eq depending on layer and T
//     }
//     else
//     {
//         //error: undefined reaction
//     }
    return 0;
}
