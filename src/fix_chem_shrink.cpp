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
#include "fix_cfd_coupling_chemistry.h"
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
    if (strncmp(style,"chem/shrink",16) == 0 && (!atom->radius_flag)||(!atom->rmass_flag))
            error -> all (FLERR,"Fix chem/shrink needs particle radius and mass");

    // defaults
    //fix_concA_          =   NULL;
    //fix_concC_          =   NULL;
    fix_changeOfA_      =   NULL;
    fix_changeOfC_      =   NULL;
    fix_rhogas          =   NULL;
    fix_tgas            =   NULL;
    fix_reactionheat_   =   NULL;
    fix_totalMole_      =   NULL;
    fix_nuField_        =   NULL;
    fix_partRe_         =   NULL;
    fix_moleFractionA_  =   NULL;
    fix_moleFractionC_  =   NULL;

    iarg_ = 3;
    k = 0;
    molMass_A_ = 0;
    molMass_C_ = 0;
    molMass_B_ = 0;
    nu_A_ = 1;
    nu_B_ = 1;
    nu_C_ = 1;
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

    if (narg < 16)
        error -> all (FLERR,"not enough arguments");
    while (iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"speciesA") == 0)
        {
            if (narg < iarg_ + 2)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesA'");
            spcA = 1;
            speciesA = new char [strlen(arg[iarg_+1])];
            strcpy(speciesA,arg[iarg_+1]);
            if(strlen(speciesA) < 1)
                error -> fix_error(FLERR,this,"speciesA not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassA") == 0)
        {
            if (spcA != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesA' before 'molMassA'");
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_A_ = atof(arg[iarg_+1]);
            if (molMass_A_ < 1)
                error -> fix_error(FLERR, this, "molar mass of A is not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nuA") == 0)
        {
            nu_A_ = atoi(arg[iarg_+1]);
            if (nu_A_ < 1)
                error -> fix_error(FLERR, this, "nuA is not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"speciesC") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "not enough arguments for 'speciesC'");
            spcC = 1;
            speciesC = new char [strlen(arg[iarg_+1])];
            strcpy(speciesC,arg[iarg_+1]);
            if (strlen(speciesC) < 1)
                error -> fix_error(FLERR, this, "speciesC not defined");
            iarg_ +=2;
            hasargs = true;
        }
       else if (strcmp(arg[iarg_],"molMassC") == 0)
        {
            if (spcC != 1)
                error -> fix_error(FLERR, this, "have to define keyword 'speciesC' before 'molMassC'");
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_C_ = atof(arg[iarg_+1]);
            if (molMass_C_ < 1)
                error -> fix_error(FLERR, this, "molMassC > 0 is required");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nuC") == 0)
        {
            nu_C_ = atoi(arg[iarg_+1]);
            if (nu_C_ < 1)
                error -> fix_error(FLERR, this, "nuC is not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"molMassB") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            molMass_B_ = atof(arg[iarg_+1]);
            if (molMass_B_ < 1)
                error -> fix_error(FLERR, this, "molMassB > 0 is required");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"nuB") == 0)
        {
            nu_B_ = atoi(arg[iarg_+1]);
            if (nu_B_ < 1)
                error -> fix_error(FLERR, this, "nuB is not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"k") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            k = atof(arg[iarg_+1]);
            if (k == 0)
                error -> fix_error(FLERR, this, "k is not defined");
            iarg_ +=2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg_],"rmin") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            rmin = atof(arg[iarg_+1]);
            if (rmin == 0)
                error -> fix_error(FLERR, this, "rmin is not defined");
            iarg_+=2;
            hasargs = true;
        }else if (strcmp(arg[iarg_],"rdef") == 0)
        {
            if (screen)
                fprintf(screen, "rmin is not given, default radius value is used! \n");
            rdef = true;
            hertzpct = atof(arg[iarg_+1]);
            if (hertzpct == 0)
                error -> fix_error(FLERR, this, "hertzpct can not be 0");
            iarg_+=3;
            hasargs = true;
        }else if (strcmp(arg[iarg_],"nevery") == 0)
        {
            nevery = atoi(arg[iarg_+1]);
            if (nevery <= 0) error->fix_error(FLERR,this,"");
            iarg_+=2;
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

    // EKI:
    // or we write directly in liggghts script instead of just speciesA as O2, instead X_O2.
    // if so the massA and C should also be changed.
    // reacting species molar fraction
    moleFracA  =   new char[n];
    strcpy(cha,"X_");
    strcat(cha,speciesA);
    strcpy(moleFracA, cha);

    // product species molar fraction
    moleFracC  =   new char[n];
    strcpy(cha,"X_");
    strcat(cha,speciesC);
    strcpy(moleFracC, cha);

    time_depend = 1;
    force_reneighbor =1;
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
        //if(fix_concA_)         modify  ->  delete_fix(speciesA);
        //if(fix_concC_)         modify  ->  delete_fix(speciesC);
        if(fix_moleFractionA_)  modify  ->  delete_fix(moleFracA);
        if(fix_moleFractionC_)  modify  ->  delete_fix(moleFracC);
        if(fix_rhogas)         modify  ->  delete_fix("partRho");
        if(fix_tgas)           modify  ->  delete_fix("partTemp");
        if(fix_reactionheat_)  modify  ->  delete_fix("reactionHeat");
        if(fix_changeOfA_)     modify  ->  delete_fix(massA);
        if(fix_changeOfC_)     modify  ->  delete_fix(massC);
        if(fix_totalMole_)     modify  ->  delete_fix("partN");
        if(fix_nuField_)       modify  ->  delete_fix("partNu");
        if(fix_partRe_)        modify  ->  delete_fix("partRe");
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
    //concA_      =   fix_concA_      -> vector_atom;
    //vector_atom = concA_;
    xA_         =   fix_moleFractionA_  -> vector_atom;

    //concC_      =   fix_concC_      -> vector_atom;
    //vector_atom = concC_;
    xC_         =   fix_moleFractionC_  ->  vector_atom;

    changeOfA_  =   fix_changeOfA_  ->  vector_atom;
    changeOfC_  =   fix_changeOfC_  ->  vector_atom;
    rhogas_     =   fix_rhogas      ->  vector_atom;
    N_          =   fix_totalMole_  ->  vector_atom;
    nuf_        =   fix_nuField_    ->  vector_atom;
    Rep_        =   fix_partRe_     ->  vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::reaction()
{
        updatePtrs();
        int nlocal  =   atom -> nlocal;
        // double dr;

        for (int i = 0 ; i < nlocal; i++)
        {
            if(radius_[i] > rmin)
            {
                if (screen)
                {
                    // fprintf(screen,"check mass frac %f \n",concA_[i]);
                    fprintf(screen,"check TimeStep %f \n",TimeStep);
                    fprintf(screen,"check current timestep %d \n", current_timestep);
                }

                // Current time step concentration change of reactant A and product C
                double dA   =   0.0;
		if(nu_A_ < 2)
		{
                    // dA   =   -k*rhogas_[i]*concA_[i]*partSurfArea(radius_[i])*TimeStep*nevery;
            dA  =   -k*xA_[i]*N_[i]*molMass_A_*partSurfArea(radius_[i])*TimeStep*nevery;
		}
		else
		{
            // dA   =   -k*pow((rhogas_[i]*concA_[i]/molMass_A_),nu_A_)*partSurfArea(radius_[i])*molMass_A_*nu_A_*TimeStep*nevery;
            dA   =   -k*pow((xA_[i]*N_[i]),nu_A_)*partSurfArea(radius_[i])*molMass_A_*nu_A_*TimeStep*nevery;
		}
                double dC   =   -dA*molMass_C_*nu_C_/(molMass_A_*nu_A_);
                // mass removed from particle
                double dB   =    dA*molMass_B_*nu_B_/(molMass_A_*nu_A_);

                // total rate of change for species A
                changeOfA_[i]       +=  dA;

                // total rate of change for species C
                changeOfC_[i]       +=  dC;

                // Mass of single particle
                pmass_[i]           +=  dB;

                // radius removed
                // dr   =   -k*concA_[i]*rhogas_[i]*TimeStep*molMass_B_*3/(molMass_A_*pdensity_[i]);

                // radius particle
                // radius_[i]          +=  dr;

                // change of radius of particle -assumption: density of particle is constant
                radius_[i]           =   pow((0.75*pmass_[i]/(M_PI*pdensity_[i])),0.333333);

                // uncomment if postproc (verification)
                if (screen)
                {
                    fprintf(screen,"Co2 Mass %f \n",changeOfC_[i]);
                    fprintf(screen,"O2 Mass %f \n",changeOfA_[i]);
                    fprintf(screen,"Particle Mass %f \n",pmass_[i]);
                    fprintf(screen,"Gas Density %f \n",rhogas_[i]);
                    fprintf(screen,"total number of moles in chem/shrink: %f \n",N_[i]);
                    // testin communication nufield and Rep
                    fprintf(screen,"nufield from DEM: %f \n", nuf_[i]);
                    fprintf(screen,"particle Reynolds from DEM: %f \n", Rep_[i]);
                }
            }
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
    // fix_concA_           =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesA,"property/atom","scalar",0,0,style));
    // fix_concC_           =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(speciesC,"property/atom","scalar",0,0,style));
    fix_moleFractionA_  =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(moleFracA,"property/atom","scalar",0,0,style));
    fix_moleFractionC_  =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(moleFracC,"property/atom","scalar",0,0,style));


    fix_changeOfA_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massA,"property/atom","scalar",0,0,style));
    fix_changeOfC_       =   static_cast<FixPropertyAtom*>(modify -> find_fix_property(massC,"property/atom","scalar",0,0,style));

    fix_tgas            =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partTemp","property/atom","scalar",0,0,style));
    fix_rhogas          =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRho","property/atom","scalar",0,0,style));
    fix_reactionheat_   =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("reactionHeat","property/atom","scalar",0,0,style));
    fix_totalMole_      =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partN","property/atom","scalar",0,0,style));
    fix_nuField_        =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partNu","property/atom","scalar",0,0,style));
    fix_partRe_         =   static_cast<FixPropertyAtom*>(modify -> find_fix_property("partRe","property/atom","scalar",0,0,style));

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

    updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixChemShrink::post_force(int)
{
    radius_ = atom ->  radius;
    pmass_  = atom ->  rmass;
    pdensity_ = atom -> density;
    TimeStep = update -> dt;
    current_timestep = update->ntimestep;

    int nlocal = atom -> nlocal;
    int i;
    
    // skip if integration not wanted at this timestep
    if (current_timestep % nevery) 
    {
        return;
    }

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
