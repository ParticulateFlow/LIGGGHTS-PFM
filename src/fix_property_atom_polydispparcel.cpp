/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_property_atom_polydispparcel.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

//NP modified C.K.
#include "pair_gran.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "group.h"
#include "timer.h"
#include "neighbor.h"

#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.001

/* ---------------------------------------------------------------------- */

FixPropertyAtomPolydispParcel::FixPropertyAtomPolydispParcel(LAMMPS *lmp, int narg, char **arg, bool parse) :
  FixPropertyAtom(lmp, narg, arg, false)
{
    if(parse) parse_args(narg,arg);
}

void FixPropertyAtomPolydispParcel::parse_args(int narg, char **arg)
{
    // Check args
    if (narg < 5) error->all(FLERR,"Illegal fix property/atom/polydispparcel command, not enough arguments");

    // Read args
    //NP 4 values for base stuff
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    data_style = FIXPROPERTY_ATOM_SCALAR;

    restart_peratom = 1;
    restart_global = 1;
 
    commGhost = 0;
    commGhostRev = 0;

    ndefaultvalues = narg - 4;
    nvalues = 1;

    int ntypes = atom->ntypes;

    defaultvalues = new double[ndefaultvalues];

    // fix handles properties that need to be initialized at particle creation
    create_attribute = 1;

    propertyname = NULL;
 
    for (int j = 0; j < ndefaultvalues; j++)
    {
        defaultvalues[j] = force->numeric(FLERR,arg[4+j]);
    }
    


    size_peratom_cols = 0;

    peratom_flag=1; //NP per-atom data is stored
    peratom_freq=1;
    extvector=0; //NP intensive

    // perform initial allocation of atom-based array
    // register with Atom class
    vector_atom = NULL; array_atom = NULL;
    grow_arrays(atom->nmax); //NP allocation of arrays
    atom->add_callback(0); //NP register that fix handles per-particle properties
    if (restart_peratom) atom->add_callback(1); //NP register that fix handles per-particles restart properties

    // init all arrays since dump may access it on timestep 0
    // or a variable may access it before first run
    //NP do not do this if caller wants to init himself (none values in defaultvalues)
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
    {
        vector_atom[i] = 1.0;
    }

    // check if there is already a fix that tries to register a property with the same name
    //NP need if(modify->fix[ifix]) check in case fix is replaced by fix with same ID
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyAtomPolydispParcel*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->fix_error(FLERR,this,"there is already a fix that registers a variable of the same name");

    // flags for vector output
    //vector_flag = 1;
    size_vector = nvalues;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomPolydispParcel::~FixPropertyAtomPolydispParcel()
{
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyAtomPolydispParcel::check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
    char errmsg[400];

    if(strcmp(varname,variablename) == 0)
    {
        // success
        return static_cast<Fix*>(this);
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

void FixPropertyAtomPolydispParcel::init()
{
    me = comm->me;

    char errmsg[300];
    int ntypes = atom->ntypes;

    if(ndefaultvalues != ntypes)
    {
        sprintf(errmsg,"Fix property/atom/polydispparcel: Length not correct for variable %s, length should be equal to %d (= number of atom types)",
                variablename,ntypes);
        error->fix_error(FLERR,this,errmsg);
    }
    
    int *types = atom->type;
    int nlocal = atom->nlocal;
    int type = 0;
    for (int i = 0; i < nlocal; i++)
    {
        type = types[i];
        vector_atom[i] = defaultvalues[type-1];
    }

}

void FixPropertyAtomPolydispParcel::pre_set_arrays()
{
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixPropertyAtomPolydispParcel::set_arrays(int i)
{
    int type = atom->type[i];
    vector_atom[i] = defaultvalues[type-1];

    // update mass of particle accordingly
    double *mass = atom->rmass;
    mass[i] *= vector_atom[i];
}



void FixPropertyAtomPolydispParcel::set_all(double value)
{
}

/* ----------------------------------------------------------------------
   set parcel factor of atom i to values and update its mass
------------------------------------------------------------------------- */
void FixPropertyAtomPolydispParcel::set_vector(int i, double value)
{
    FixPropertyAtom::set_vector(i, value);
    double *mass = atom->rmass;
    mass[i] *= vector_atom[i];
}