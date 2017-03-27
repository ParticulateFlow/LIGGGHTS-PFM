/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
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
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "particleToInsert.h"
#include <math.h>
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "modify.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) : Pointers(lmp)
{
        groupbit = 0;

        nspheres = ns;

        memory->create(x_ins,nspheres,3,"x_ins");
        radius_ins = new double[nspheres];

        atom_type_vector = new int[nspheres];
        atom_type_vector_flag = false;

        fix_property = 0;
        fix_property_value = 0.;
        fix_property_ivalue = 0;
        property_index = -1;
        property_iindex = -1;
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy(x_ins);
        delete []radius_ins;
        delete []atom_type_vector;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    for(int i = 0; i < nspheres; i++)
    {
        /*NL*/ //if (screen) fprintf(screen,"proc %d tyring to insert particle at pos %f %f %f\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
        //NP do not need subdomain check any longer since have processor-local lists anyway
        //if (domain->is_in_extended_subdomain(x_ins[i]))
        //{
                /*NL*/ //if (screen) fprintf(screen,"   proc %d inserting particle at pos %f %f %f\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
                inserted++;
                if(atom_type_vector_flag)
                    atom->avec->create_atom(atom_type_vector[i],x_ins[i]);
                else
                    atom->avec->create_atom(atom_type,x_ins[i]);
                int m = atom->nlocal - 1;
                atom->mask[m] = 1 | groupbit;
                vectorCopy3D(v_ins,atom->v[m]);
                vectorCopy3D(omega_ins,atom->omega[m]);
                atom->radius[m] = radius_ins[0];
                atom->density[m] = density_ins;
                atom->rmass[m] = mass_ins;

                //pre_set_arrays() called via FixParticleDistribution
                for (int j = 0; j < nfix; j++)
                    if (fix[j]->create_attribute) fix[j]->set_arrays(m);

                // apply fix property setting coming from fix insert
                // this overrides the set_arrays call above
                if(fix_property)
                    fix_property->vector_atom[m] = fix_property_value;
                else if(property_iindex >= 0)
                    atom->ivector[property_iindex][m] = fix_property_ivalue;
                else if(property_index >= 0)
                    atom->dvector[property_index][m] = fix_property_value;
        //}
    }

    return inserted;
}

/* ---------------------------------------------------------------------- */
//NP checks against xnear list
//NP returns # inerted spheres
//NP if >1, increase nbodies

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    if(nspheres > 1)
        error->one(FLERR,"check_near_set_x_v_omega not implemented yet for nspheres>1");

    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        /*NL*///if (screen) printVec3D(screen,"x_ins[0]",x_ins[0]);
        /*NL*///if (screen) fprintf(screen,"radius_ins[0] %f xnear[i][3] %f \n",radius_ins[0],xnear[i][3]);
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

    return 1;
}

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList)
{
    if(nspheres > 1)
        error->one(FLERR,"check_near_set_x_v_omega not implemented yet for nspheres>1");

    vectorCopy3D(x,x_ins[0]);

    if(neighList.hasOverlap(x_ins[0], radius_ins[0], get_atom_type())) {
        return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    neighList.insert(x_ins[0], radius_ins[0], get_atom_type());

    return 1;
}

/* ---------------------------------------------------------------------- */
//NP no check against xnear list
//NP returns # inserted spheres
//NP if >1, increase nbodies

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double rel[3];

    // add insertion position
    // relative position of spheres to each other already stored at this point
    // also take quat into account
    for(int j = 0; j < nspheres; j++)
    {
        // if only one sphere, then x_bound = x_ins
        if(1 == nspheres)
            vectorAdd3D(x_ins[j],x,x_ins[j]);
        else
        {
            vectorSubtract3D(x_ins[j],x_bound_ins,rel);
            MathExtraLiggghts::vec_quat_rotate(rel,quat);
            vectorAdd3D(rel,x,x_ins[j]);
        }
    }

    // set velocity and omega
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    return nspheres;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsert::scale_pti(double r_scale)
{
    double r_scale3 = r_scale*r_scale*r_scale;

    for(int i = 0; i < nspheres; i++)
    {
        radius_ins[i] *= r_scale;
        vectorScalarMult3D(x_ins[i],r_scale);
    }

    volume_ins *= r_scale3;
    mass_ins *= r_scale3;

    r_bound_ins *= r_scale;
}
