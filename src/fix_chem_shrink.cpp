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
#include "region.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_chem_shrink.h"
#include "pair_gran.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "properties.h"
#include "property_registry.h"
#include "global_properties.h"
#include "force.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrink::FixChemShrink(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp,narg,arg)
{
    if (strncmp(style,"chem/shrink",14) == 0 && (!atom->radius_flag)||(!atom->rmass_flag))
            error -> all (FLERR,"Fix chem/shrink needs particle radius and mass");

    // defaults
    fix_concA_          =   NULL;
    fix_concC_          =   NULL;
    fix_changeOfA_      =   NULL;
    fix_changeOfC_      =   NULL;
    fix_rhogas_         =   NULL;
    fix_tgas_           =   NULL;
    fix_reactionheat_   =   NULL;

    iarg_ = 3;
    k = 0;
    molMass_A_ = 0;
    molMass_C_ = 0;
    molMass_B_ = 0;
    pmass_ = NULL;
    rmin = 0.;
    rdefault = 0.;
    hertzpct = 0.;
    spcA = 0;
    spcC = 0;
    Yeff_ = NULL;
    int n = 16;
    char cha[30];

    bool hasargs = true;
    rdef = false;
    while (iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            spcA = 1;
            iarg_++;
            speciesA = new char [strlen(arg[iarg_]) + 1];
            strcpy(speciesA,arg[iarg_]);
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassA") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            if (spcA != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesA' before 'molMassA'");
            if (strlen(speciesA) < 1)
                error -> fix_error(FLERR, this, "speciesA not defined");
            iarg_++;
            molMass_A_ = atof(arg[iarg_]);
            if (molMass_A_ < 1)
                error -> fix_error(FLERR, this, "molar mass of A is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            spcC = 1;
            iarg_++;
            speciesC = new char [strlen(arg[iarg_])+1];
            strcpy(speciesC,arg[iarg_]);
            iarg_ ++;
            hasargs = true;
        }
       else if (strcmp(arg[iarg_],"molMassC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            if (spcC != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesC' before 'molMassC'");
            if (strlen(speciesC) < 1)
                error -> fix_error(FLERR, this, "speciesC not defined");
            iarg_++;
            molMass_C_ = atof(arg[iarg_]);
            if (molMass_C_ < 1)
                error -> fix_error(FLERR, this, "molar mass of C is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassB") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_++;
            molMass_B_ = atof(arg[iarg_]);
            if (molMass_B_ < 1)
                error -> fix_error(FLERR, this, "molar mass of B is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"k") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_ ++;
            k = atof(arg[iarg_]);
            if (k == 0)
                error -> fix_error(FLERR, this, "k is not defined");
            iarg_ ++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"rmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            iarg_++;
            rmin = atof(arg[iarg_]);
            if (rmin == 0)
                error -> fix_error(FLERR, this, "rmin is not defined");
            iarg_++;
            hasargs = true;
        }else if (strcmp(arg[iarg_],"rdef") == 0)
        {
            if (screen)
                fprintf(screen, "rmin is not given, default radius value is used! \n");
            rdef = true;
            iarg_++;
            if (strcmp(arg[iarg_],"hertzpct") != 0)
                error -> fix_error(FLERR, this, "please enter a percentage to caclculate the hertz timestep for");
            iarg_++;
            hertzpct = atof(arg[iarg_]);
            if (hertzpct == 0)
                error -> fix_error(FLERR, this, "hertzpct can not be 0");
            iarg_++;
            hasargs = true;
        }
    }

    // define changed species mass A
    massA = new char [n];
    strcpy(cha,"Modified_");
    strcat(cha,speciesA);
    strcpy(massA,cha);

    // define changed species mass C
    massC = new char [n];
    strcpy(cha,"Modified_");
    strcat(cha,speciesC);
    strcpy(massC,cha);

    time_depend = 1;
    force_reneighbor =1;
    nevery = 1;
    next_reneighbor = update ->ntimestep + nevery;

    restart_global = 1;
}


/* ---------------------------------------------------------------------- */

FixChemShrink::~FixChemShrink()
{
    delete massA;
    delete massC;

    delete speciesA;
    delete speciesC;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        if(fix_concA_)         modify  ->  delete_fix(speciesA);
        if(fix_concC_)         modify  ->  delete_fix(speciesC);
        if(fix_rhogas_)        modify  ->  delete_fix("partRho");
        if(fix_tgas_)          modify  ->  delete_fix("partTemp");
        if(fix_reactionheat_)  modify  ->  delete_fix("reactionHeat");
        if(fix_changeOfA_)     modify  ->  delete_fix(massA);
        if(fix_changeOfC_)     modify  ->  delete_fix(massC);
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
    concA_      =   fix_concA_      -> vector_atom;
    concC_      =   fix_concC_      -> vector_atom;
    changeOfA_  =   fix_changeOfA_  -> vector_atom;
    changeOfC_  =   fix_changeOfC_  -> vector_atom;
    rhogas_     =   fix_rhogas_     -> vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::reaction()
{
        updatePtrs();
        int nlocal  =   atom -> nlocal;
        double dr;
        double totalChangeOfC, totalRhogas, totalpmass;
        double aveChangeOfC, aveRhogas, avePmass;

        for (int i = 0 ; i < nlocal; i++)
        {
            if(radius_[i] > rmin)
            {
                // Current time step concentration change of reactant A and product C
                double dA   =   -k*rhogas_[i]*concA_[i]*partSurfArea(radius_[i])*TimeStep;
                double dC   =   -dA*(molMass_C_/molMass_A_);
                // mass removed from particle
                double dB   =    dA*(molMass_B_/molMass_A_);

                // total rate of change for species A
                changeOfA_[i]       +=  dA;

                // total rate of change for species C
                changeOfC_[i]       +=  dC;

                // Mass of single particle
                pmass_[i]           +=  dB;

                // radius removed
                dr   =   -k*concA_[i]*rhogas_[i]*TimeStep*molMass_B_/(molMass_A_*pdensity_[i]);

                // radius particle
                radius_[i]          +=  dr;

                // change of radius of particle -assumption: density of particle is constant
                //radius_[i]           =   pow((0.75*pmass_[i]/(M_PI*pdensity_[i])),0.333333);

                // add up the values that are going to be printed
                totalChangeOfC += changeOfC_[i];
                totalpmass     += pmass_[i];
                totalRhogas    += rhogas_[i];

                // take the mean over every particle
                aveChangeOfC   = totalChangeOfC/nlocal;
                aveRhogas      = totalRhogas/nlocal;
                avePmass       = totalpmass/nlocal;
            }

            if(screen)
            {
                fprintf(screen,"total change of C = %f \n", totalChangeOfC);
                fprintf(screen,"total change of rhogas = %f \n", totalRhogas);
                fprintf(screen,"total change of pmass = %f \n", totalpmass);
                fprintf(screen,"mean CO2 change = %f \n", aveChangeOfC);
                fprintf(screen,"mean rhogas change = %f \n", aveRhogas);
                fprintf(screen,"mean pmass change = %f \n", avePmass);
            }

            //if (comm -> me == 0 && screen)
            //  //fprintf(screen, "rhogas: %f \n", rhogas_[i];
        }
}

/* ----------------- compute particle surface area ------------------------ */

double FixChemShrink::partSurfArea(double radius)
    {
        double A_p =   4*M_PI*radius*radius;
        return (A_p);
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
    fix_concA_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesA,"property/atom","scalar",0,0,style));
    fix_concC_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesC,"property/atom","scalar",0,0,style));
    fix_changeOfA_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massA,"property/atom","scalar",0,0,style));
    fix_changeOfC_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massC,"property/atom","scalar",0,0,style));

    fix_tgas_        =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_=   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));

    if (rdef)
    {
        PairGran *pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
        if(!pair_gran)
            error->fix_error(FLERR,this,"'area_correction' requires using a granular pair style");
        int max_type = pair_gran->get_properties()->max_type();

        Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
        nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();

        force->registry.registerProperty("Yeff_", &MODEL_PARAMS::createYeff);
        force->registry.connect("Yeff_", Yeff_,this->style);
    }
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;
    TimeStep = update -> dt;
    int nlocal = atom -> nlocal;
    int i;

    if (rdef)   default_radius();

    // allocate and initialize deletion list and shrink list
    memory->create(dlist,nlocal,"delete_atoms:dlist");
    for (i = 0; i < nlocal; i++) dlist[i] = 0;

    reaction();

    delete_atoms();

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

/* ---------------------------------------------------------------------- */

void FixChemShrink::delete_atoms()
{
    int nlocal = atom->nlocal;
    int i = 0;

    while (i < nlocal) {
        dlist[i] = radius_[i] < rmin;
        if (dlist[i]) {
        atom -> avec -> copy(nlocal-1,i,1);
        dlist[i] = dlist[nlocal-1];
        nlocal--;
        } else i++;
      }

      atom->nlocal = nlocal;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::default_radius()
{
    int nlocal  = atom->nlocal;
    double **v  = atom -> v;
    double vmag;
    int *mask   = atom -> mask;
    int *type   = atom -> type;

    double m_;
    double denom;
    double numer;

    vmax_ = 0;

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            vmag = sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
            if(vmag > vmax_) vmax_=vmag;
        }
    }

    MPI_Max_Scalar(vmax_,world);
    vmax_ = (2.*vmax_);

    if (screen)
    {
        fprintf(screen,"vmax value is: %f \n", vmax_);
    }

    PairGran *pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
    int max_type = pair_gran->get_properties()->max_type();

    for (int a = 1; a < max_type + 1; a++)
    {
        for (int b = a; b < max_type + 1; b++)
        {
            const double Eeff = Yeff_[a][b];

            for (int i = 0; i < nlocal; i++)
            {
                if (mask[i] & groupbit)
                {
                    if (type[i]!=a || type[i]!=b) continue;
                    // Hertzain time step can be calculated with
                    // 2.87 * radius *( m / (Y*vmax))^0.2
                    // here it is considered that Hertzian ts is 20 percent of normal dt
                    // with this the minimum radius is calculated
                    m_ = 4 * M_PI / (3 * pdensity_[i]);
                    denom = Eeff*Eeff*vmax_;
                    numer = 2*m_*m_;
                    rdefault = (hertzpct / 2.87)  * TimeStep * pow(numer/denom,-0.2);
                }
            }
        }
    }

    MPI_Max_Scalar(rdefault,world);
    rmin = rdefault;
}
