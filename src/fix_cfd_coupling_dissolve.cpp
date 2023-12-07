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
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_dissolve.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingDissolve::FixCfdCouplingDissolve(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling_ = NULL;
    fix_convectiveFlux_ = NULL;
    rmin_ = -1.0;

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
            rmin_ = atof(arg[iarg]);
            iarg++;
            hasargs = true;
        }
    }

    if (rmin_ < 0.) error->all(FLERR,"Fix couple/cfd/dissolve: rmin must be >= 0. Specify a meaningful value.");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingDissolve::~FixCfdCouplingDissolve()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::pre_delete(bool unfixflag)
{
    if (fix_convectiveFlux_) modify->delete_fix("convectiveMassFlux");
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
    if (!fix_convectiveFlux_)
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
        fix_convectiveFlux_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::init()
{
    // make sure there is only one fix of this style
    if (modify->n_fixes_style(style) != 1)
        error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if (!fix_coupling_)
        error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    // values to come from OF
    fix_coupling_->add_pull_property("convectiveMassFlux","scalar-atom");

    fix_convectiveFlux_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveMassFlux","property/atom","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDissolve::post_force(int)
{
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double dt = update->dt;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *density = atom->density;

    double vmin = MY_4PI3*rmin_*rmin_*rmin_;

    // communicate convective flux to ghosts, there might be new data

    if (0 == neighbor->ago)
    {
        fix_convectiveFlux_->do_forward_comm();
    }

    double *convectiveFlux = fix_convectiveFlux_->vector_atom;


    //NP add convective flux to per-particle mass flux

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if (radius[i] > rmin_)
            {
                rmass[i] += convectiveFlux[i]*dt;

                if (rmass[i] <= (vmin*density[i]))
                {
                    // set to minimum radius and mark for delete
                    rmass[i] = vmin*density[i];
                    radius[i] = rmin_;
                    atom_tags_delete_.push_back(atom->tag[i]);
                    //fprintf(screen, "Marking atom %i for deletion \n",atom->tag[i]);
                }
                else
                {
                    // set new radius
                    radius[i] = cbrt(0.75*rmass[i]/(MY_PI*density[i]));
                }
            }
        }
    }
}

void FixCfdCouplingDissolve::pre_exchange()
{
    // delete atoms

    int nparticles_deleted_this = 0;
    AtomVec *avec = atom->avec;
    int nlocal = atom->nlocal;

    while (atom_tags_delete_.size() > 0)
    {
        int iPart = atom->map(atom_tags_delete_[0]);

        if(iPart >= 0)
        {
            nparticles_deleted_this++;

            avec->copy(nlocal-1,iPart,1);

            // update atom map
            // need to do this since atom map is needed to get local index for deletion
            atom->map_one(atom->tag[nlocal-1], iPart);

            nlocal--;
        }
        else
        {
            // particle may have been removed already by a different deleting command
            error->fix_warning(FLERR, this, "failed to find atom for deletion (possibly already deleted by another deleting command)");
        }

        atom_tags_delete_.erase(atom_tags_delete_.begin());
    }

    atom_tags_delete_.clear();

    MPI_Sum_Scalar(nparticles_deleted_this,world);

    //NP tags and maps
    atom->nlocal = nlocal;

    //fprintf(screen,"nparticles_deleted_this_ = %i \n", nparticles_deleted_this_);

    bigint nblocal = atom->nlocal;
    MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

    if (nparticles_deleted_this)
    {
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
}
