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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_sph_general_simple.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallSphGeneralSimple::FixWallSphGeneralSimple(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg)
{
    if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments.");

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        if(strcmp(arg[iarg_++],"r0"))
            error->fix_error(FLERR,this,"illegal argument, expecting keyword 'r0'");
        r0 = force->numeric(arg[iarg_++]);
        if(strcmp(arg[iarg_++],"D"))
            error->fix_error(FLERR,this,"illegal argument, expecting keyword 'D'");
        D  = force->numeric(arg[iarg_++]);
    }

    if (r0 <= 0. || D < 0.)
      error->fix_error(FLERR,this,"values for r0 or D are invalid");

    set_r0(r0);
}

FixWallSphGeneralSimple::~FixWallSphGeneralSimple()
{

}

/* ----------------------------------------------------------------------
   compute force on particle
------------------------------------------------------------------------- */

void FixWallSphGeneralSimple::compute_force(int ip,double deltan,double rsq,double meff_wall,
                            double dx,double dy,double dz,double *vwall,
                            double *c_history,double area_ratio)
{
    // do not use deltan here for SPH

    double **f = atom->f;
    double fwall,rinv, r;
    double frac,frac2,frac4; // for penetration force

    if (rsq == 0.) return; // center of the cylinder ... no repulsive force!

    r = sqrt(rsq);
    rinv = 1./r;

    // repulsive penetration force
    if (r <= r0) {

        frac = r0*rinv;
        frac2 = frac*frac;
        frac4 = frac2*frac2;

        fwall = D * (frac4 - frac2) * rinv; // second rinv

        f[ip][0] += fwall * dx;
        f[ip][1] += fwall * dy;
        f[ip][2] += fwall * dz;
    }

}
