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

#include <string.h>
#include <math.h>
#include "error.h"
#include "fix_chem_shrink_Arrhenius.h"
#include "fix_property_atom.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixChemShrinkArrhenius::FixChemShrinkArrhenius(LAMMPS *lmp, int narg, char **arg) :
    FixChemShrink(lmp,narg,arg),
    T0(0.0)
{
    // ATTENTION: keyword T needs to come last in input script in fix/chem/shrink/Arrhenius
    while (iarg_ < narg)
    {
        if (strcmp(arg[iarg_],"T") == 0)
        {
            if (iarg_ + 2 > narg)
                error -> fix_error(FLERR, this, "Wrong number of arguments");
            T0 = atof(arg[iarg_+1]);
            if (T0 <= 0.)
                error -> fix_error(FLERR, this, "T is not (well-)defined");
            iarg_ +=2;
        }
        else
        {
            iarg_++;
        }
    }
}

/* ---------------------------------------------------------------------- */

FixChemShrinkArrhenius::~FixChemShrinkArrhenius()
{

}

/* ---------------------------------------------------------------------- */

void FixChemShrinkArrhenius::updatePtrs()
{
    FixChemShrink::updatePtrs();
    tgas_ = fix_tgas->vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixChemShrinkArrhenius::reactionRatConst(int i)
{
    const double t = tgas_[i];

    if(t < 100.0)
    {
        return 0.0;
    }
    else
    {
        return k0*exp(-T0/t);
    }
}

/* ---------------------------------------------------------------------- */

