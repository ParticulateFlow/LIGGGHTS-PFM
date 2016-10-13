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
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::FixChemShrinkCore(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg),
    nmaxlayers_(3),
    rmin_(0.01),
    drmin_(0.001)
{
    if (strncmp(style,"chem/shrink/core",14) == 0 && (!atom->radius_flag)||(!atom->rmass_flag))
            error -> all (FLERR,"Fix chem/shrink needs particle radius and mass");

    // defaults
    fix_concA_   =   NULL;
    fix_concC_   =   NULL;
    fix_changeOfA_   =   NULL;
    fix_changeOfC_   =   NULL;
    fix_rhogas_      =   NULL;
    // fix_tpart_       =   NULL;
    // fix_reactionheat_    =   0;

    iarg_ = 3;
    if (narg < 14)
        error -> all (FLERR,"not enough arguments");

    // check and define species A
    if (strcmp(arg[iarg_++],"speciesA") != 0)
        error -> all (FLERR, "missing keyword speciesA");
    speciesA = new char [strlen(arg[iarg_])+1];
    strcpy(speciesA, arg[iarg_++]);
    if (speciesA == 0)
        error -> all (FLERR, "speciesA not defined");

    // check and define molecular mass of A
    if (strcmp(arg[iarg_++],"molMassA") != 0)
        error -> all (FLERR, "keyword molMassA for species A is missing");
    molMass_A_  =   atof(arg[iarg_++]);
    if (molMass_A_ < 1)
        error -> all (FLERR,"molMass species A is not defined");

    // check and define Species C
    if (strcmp(arg[iarg_++],"speciesC") != 0)
        error -> all (FLERR, "missing keyword speciesC");
    speciesC = new char [strlen(arg[iarg_])+1];
    strcpy(speciesC, arg[iarg_++]);
    if (speciesC == 0)
        error -> all (FLERR, "speciesC not defined");

    // check and define molecular mass of C
    if (strcmp(arg[iarg_++],"molMassC") != 0)
        error -> all (FLERR, "keyword molMassC for species C is missing");
    molMass_C_  =   atof(arg[iarg_++]);
    if (molMass_C_ < 1)
        error -> all (FLERR,"molMass species C is not defined");

    

    // flags for vector output
    peratom_flag =  1;  //0/1 per-atom data is stored
    peratom_freq =  1;
    scalar_flag =   1;
    global_freq =   1;
    extscalar   =   1;
    
    
    

}

/* ---------------------------------------------------------------------- */

FixChemShrinkCore::~FixChemShrinkCore()
{
    if(fix_changeOfA_)  delete  []fix_changeOfA_;
    if(fix_changeOfC_)  delete  []fix_changeOfC_;

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_concA_) modify  ->  delete_fix(speciesA);
    if(unfixflag && fix_concC_) modify  ->  delete_fix(speciesC);
    if(unfixflag && fix_rhogas_)    modify  -> delete_fix("partRho");
    // if(unfixflag && fix_tgas_)  modify  ->  delete_fix("partTemp");

    if(unfixflag && fix_changeOfA_) modify  ->  delete_fix("changeOfSpeciessMass_A");
    if(unfixflag && fix_changeOfC_) modify  ->  delete_fix("changeOfSpeciesMass_C");


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
}

/* ---------------------------------------------------------------------- */
void FixChemShrinkCore::updatePtrs()
{
    changeOfA_  =   fix_changeOfA_  -> vector_atom;
    changeOfC_  =   fix_changeOfC_  -> vector_atom;
    rhogas_     =   fix_rhogas_     -> vector_atom;
    concA_      =   fix_concA_      -> vector_atom;
    concC_      =   fix_concC_      -> vector_atom;
    
    relRadii_   =   fix_layerRelRad_-> array_atom;
    layerDensities_ = fix_dens_     -> get_values();
    layerMolMasses_ = fix_molMass   -> get_values();
    
    k0_             = fix_k0_       -> get_values();
    // other chemical parameters

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::init()
{
    if (!atom -> radius_flag || !atom -> density_flag)
        error -> all(FLERR,"Fix chem/shrink can only be used with sphere atom style");

    // references
    
    
    // look up reaction properties
    // they are defined by the user via FixPropertyGlobal
    // the fixes' names have to be chosen such that they can be identified by the reaction fix they correspond to
    int n = strlen(id) + 1 + 3;
    char* fixname = new char[n];
    strcpy (fixname,"k0_");
    strcat(fixname,id);      
    fix_k0_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property(fixname,"property/global","vector",0,0,"FixChemShrinkCore"));
    delete []fixname;
      
    // look up *fix_dens_, *fix_molMass
    // their names are of the form property + material, e.g. "molMass_Ore"
    
    // look up *fix_layerRelRad_;
    // name something like "layerradii"
    
    
    fix_concA_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesA,"property/atom","scalar",0,0,id));
    fix_concC_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesC,"property/atom","scalar",0,0,id));
    fix_changeOfA_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("changeOfSpeciessMass_A","property/atom","scalar",0,0,id));
    fix_changeOfC_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("changeOfSpeciessMass_C","property/atom","scalar",0,0,id));
    fix_tgas_        =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,id));
    fix_rhogas_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,id));
    fix_reactionheat_=   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,id));
    updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixChemShrinkCore::post_force(int)
{
    updatePtrs();
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;

    int nlocal  =   atom->nlocal;
    int *mask   =   atom->mask;
    
    double dm[nmaxlayers_];
    double a[nmaxlayers_];
    double diff[nmaxlayers_];
    double masst;
    double y0[nmaxlayers_];

    for (int i = 0; i<nlocal; i++)
    {
        if (mask[i] & groupbit)
	{
//	    getA(i,a);
//	    getDiff(i,diff);
//	    getMassT(i,masst);
//	    getY0(i,y0);
	    
	    for(int j = 0; j<nmaxlayers_; j++)
                dm[j] = 0.0;
	    reaction(i,a,diff,masst,y0,dm);
	    
	    update_atom_properties(i,dm);
	    update_gas_properties(i,dm);
	}
    }
}



/* ----------------- compute particle surface area ------------------------ */
 double FixChemShrinkCore::partSurfArea(double radius)
{
        double A_p =   4*M_PI*radius*radius;
        return (A_p);
};

    
// returns number of active layers    
int FixChemShrinkCore::active_layers(int i)
{
    int layers;
    
    for(layers=0; layers<nmaxlayers_; layers++)
        if(fix_layerRelRad_->array_atom[i][layers] < rmin_)
	    break;

    return layers;
}

void FixChemShrinkCore::update_atom_properties(int i, double *dm)
{
 // based on material change: update relative radii, average density and mass of atom i
  
 // WARNING 1: do not remove more material of layer J than present --> decrease dm[J] if necessary
  
 // WARNING 2: if Wustite layer thickness < rmin AND T < 567 C, direct conversion of Magnetite to Fe is possible
 //            this means that Wu-Fe radius is shrinking together with the Ma-Wu radius 
 //            in this case, the equilibrium value and reaction parameters at the Ma layer need to be adapted within
 //            IMPLEMENTATION ONLY IF NECESSARY
}

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

void FixChemShrinkCore::getA(int i, double *a)
{
  // if layer J thickness < drmin or layer J radius < rmin, a[J] = LARGE
  // if T < Tcrit1, a = LARGE

}