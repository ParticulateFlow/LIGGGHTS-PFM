/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2021- Eindhoven University of Technology

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Tim Nijssen (TU/e)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "math_const.h"
#include "fix_cfd_coupling_dissolve.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingDissolve::FixCfdCouplingDissolve(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_convectiveFlux = NULL;
    rmin = -1.0;

    int iarg = 3;

    bool hasargs = true;
    while (iarg < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg],"rmin") == 0)
        {
            if (narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'rmin'");
            iarg++;
            rmin = atof(arg[iarg]);
            iarg++;
            hasargs = true;
        }
    }

    if (rmin < 0.) error->all(FLERR,"Fix couple/cfd/dissolve: rmin must be >= 0. Specify a meaningful value.");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingDissolve::~FixCfdCouplingDissolve()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::pre_delete(bool unfixflag)
{
    if (fix_convectiveFlux) modify->delete_fix("convectiveMassFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingDissolve::setmask()
{
    int mask = 0;
    mask |= PRE_EXCHANGE;
    mask |= POST_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::post_create()
{
    // register convective flux
    if (!fix_convectiveFlux)
    {
        const char* fixarg[11];
        fixarg[0]="convectiveMassFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="convectiveMassFlux";
        fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
        fixarg[5]="no";     //NP restart
        fixarg[6]="yes";    //NP communicate ghost
        fixarg[7]="no";     //NP communicate rev
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::init()
{
    // make sure there is only one fix of this style
    if (modify->n_fixes_style(style) != 1)
        error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if (!fix_coupling)
        error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    // values to come from OF
    fix_coupling->add_pull_property("convectiveMassFlux","scalar-atom");

    fix_convectiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveMassFlux","property/atom","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::post_force(int)
{
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double dt = update->dt;
    double *radius_ = atom->radius;
    double *pmass_ = atom->rmass;
    double *pdensity_ = atom->density;

    double vmin = MY_4PI3*rmin*rmin*rmin;

    // communicate convective flux to ghosts, there might be new data

    if (0 == neighbor->ago)
    {
        fix_convectiveFlux->do_forward_comm();
    }

    double *convectiveFlux = fix_convectiveFlux->vector_atom;


    //NP add convective flux to per-particle mass flux

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (radius_[i] > rmin)
            {
                pmass_[i] += convectiveFlux[i]*dt;

                if (pmass_[i] <= (vmin*pdensity_[i]))
                {
                    // set to minimum radius and mark for delete
                    pmass_[i] = vmin*pdensity_[i];
                    radius_[i] = rmin;
                    atom_tags_delete_.push_back(atom->tag[i]);
                    //fprintf(screen, "Marking atom %i for deletion \n",atom->tag[i]);
                }
                else
                {
                    // set new radius
                    radius_[i] = cbrt(0.75*pmass_[i]/(MY_PI*pdensity_[i]));
                }
            }
        }
    }
}

void FixCfdCouplingDissolve::pre_exchange()
{
    delete_atoms();
}

void FixCfdCouplingDissolve::delete_atoms()
{
    int nparticles_deleted_this_ = 0.;
    int *atom_map_array = atom->get_map_array();
    AtomVec *avec = atom->avec;
    int nlocal = atom->nlocal;

    while (atom_tags_delete_.size() > 0)
    {
        int iPart = atom->map(atom_tags_delete_[0]);

        avec->copy(nlocal-1,iPart,1);

        nparticles_deleted_this_++;
        nlocal--;

        //NP manipulate atom map array
        //NP need to do this since atom map is needed for deletion
        atom_map_array[atom->tag[atom->nlocal-1]] = iPart;

        atom_tags_delete_.erase(atom_tags_delete_.begin());
    }

    atom_tags_delete_.clear();

    //NP tags and maps
    atom->nlocal = nlocal;

    //fprintf(screen,"nparticles_deleted_this_ = %i \n", nparticles_deleted_this_);

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

    if (atom->tag_enable)
    {
        if (atom->map_style)
        {
            atom->nghost = 0;
            atom->map_init();
            atom->map_set();
        }
    }
}
