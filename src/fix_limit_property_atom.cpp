/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

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
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_limit_property_atom.h"
#include "atom.h"
#include "memory.h"
#include "error.h"

//NP modified C.K.
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"

#include "mpi_liggghts.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLimitPropertyAtom::FixLimitPropertyAtom(LAMMPS *lmp, int narg, char **arg, bool parse) :
  Fix(lmp, narg, arg)
{
    /*NL*/ //if (screen) fprintf(screen,"HERE parse for id %s\n",id);
    if(parse) parse_args(narg,arg);
}

void FixLimitPropertyAtom::parse_args(int narg, char **arg)
{
    int n = strlen(arg[3]) + 1;
    targetfixname = new char[n];
    strcpy(targetfixname,arg[3]);
    fix_target_ = NULL;

    if (narg < 8 || (narg - 6) % 2)
    {
        error->all(FLERR,"Wrong number of arguments in fix limit/property/atom.");
    }
    nvalues = (narg - 6)/2;
    maxvalues = new double[nvalues];
    minvalues = new double[nvalues];

    if (strcmp(arg[4],"min") == 0)
    {
        for (int j = 0; j < nvalues; j++)
        {
            minvalues[j] = atof(arg[5+j]);
        }
    }
    else error->all(FLERR,"Expected to find minimum limiting values");

    if (strcmp(arg[5+nvalues],"max") == 0)
    {
        for (int j = 0; j < nvalues; j++)
        {
            maxvalues[j] = atof(arg[6+nvalues+j]);
        }
    }
    else error->all(FLERR,"Expected to find maximum limiting values");
}

/* ---------------------------------------------------------------------- */

FixLimitPropertyAtom::~FixLimitPropertyAtom()
{
    delete[] targetfixname;
    delete[] maxvalues;
    delete[] minvalues;
}

/* ---------------------------------------------------------------------- */

void FixLimitPropertyAtom::init()
{
    if (nvalues == 1)
    {
        fix_target_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(targetfixname, "property/atom", "scalar", 0, 0, "fix/limit/property/atom"));
    }
    else
    {
        fix_target_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(targetfixname, "property/atom", "vector", 0, 0, "fix/limit/property/atom"));
    }

    if (fix_target_->num_defaultvalues() != nvalues)
    {
        char errmsg[400];
        sprintf(errmsg,"Number of values to be limited by FixLimitPropertyAtom does not match that stored in %s",targetfixname);
        error->all(FLERR,errmsg);    
    }
}

/* ---------------------------------------------------------------------- */

int FixLimitPropertyAtom::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLimitPropertyAtom::end_of_step()
{
    int nlocal = atom->nlocal;

    if (nvalues == 1)
    {
        double *targetvalue = fix_target_->vector_atom;
        for (int i = 0; i < nlocal; ++i)
        {
            if (targetvalue[i] < minvalues[0])
            {
                targetvalue[i] = minvalues[0];
            }
            else if (targetvalue[i] > maxvalues[0])
            {
                targetvalue[i] = maxvalues[0];
            }
        }
    }
    else
    {
        double **targetvalue = fix_target_->array_atom;
        for (int i = 0; i < nlocal; ++i)
        {
            for (int j = 0; j < nvalues; j++)
            {
                if (targetvalue[i][j] < minvalues[j])
                {
                    targetvalue[i][j] = minvalues[j];
                }
                else if (targetvalue[i][j] > maxvalues[j])
                {
                    targetvalue[i][j] = maxvalues[j];
                }
            }
        }
    }
}
