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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include "particleToInsert.h"
#include <math.h>
#include <limits.h>
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "math_const.h"
#include "modify.h"
#include "force.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) :
  Pointers(lmp),
  nspheres(ns),
  groupbit(0),
  atom_type(0),
  bond_type(0),
  density_ins(0.0),
  volume_ins(0.0),
  mass_ins(0.0),
  r_bound_ins(0.0),
  atom_type_vector_flag(false),
  fix_property(NULL),
  n_fix_property(0),
  fix_property_nentry(NULL),
  fix_property_value(NULL),
  local_start(-1),
  needs_bonding(false)
{
    memory->create(x_ins,nspheres,3,"x_ins");
    radius_ins = new double[nspheres]();
    atom_type_vector = new int[nspheres]();
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
    memory->destroy(x_ins);
    delete []radius_ins;
    delete []atom_type_vector;

    if (fix_property_value)
    {
        for (int i = 0; i < n_fix_property; i++)
            delete [] fix_property_value[i];
        delete [] fix_property_value;
    }

    delete [] fix_property;
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

    local_start = atom->nlocal;
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
                atom->radius[m] = radius_ins[i];
                atom->density[m] = density_ins;
                atom->rmass[m] = (nspheres==1)? (mass_ins) : (MathConst::MY_4PI3*radius_ins[i]*radius_ins[i]*radius_ins[i]*density_ins);

                //pre_set_arrays() called above
                for (int j = 0; j < nfix; j++)
                   if (fix[j]->create_attribute) fix[j]->set_arrays(m);

                // apply fix property setting coming from fix insert
                // this overrides the set_arrays call above
                if(fix_property)
                {
                    for (int j = 0; j < n_fix_property; j++)
                    {
                        if (fix_property_nentry[j] == 1)
                        {
                            fix_property[j]->vector_atom[m] = fix_property_value[j][0];
                            if(strcmp(fix_property[j]->id,"bond_random_id") == 0)
                            {
                                if (atom->molecule_flag)
                                {
                                    needs_bonding = true;
                                    // use random part as base for dummy molecule ID
                                    // see FixTemplateMultiplespheres::randomize_ptilist
                                    double dmol = (fix_property_value[j][0] - static_cast<double>(update->ntimestep));
                                    if(dmol > 1.0 || dmol < 0.0)
                                        error->one(FLERR, "Internal error (particle to insert: mol id)");
                                    atom->molecule[m] = -static_cast<int>(dmol * INT_MAX);
                                    // actual molecule value needs to be created afterwards via atom->mol_extend()
                                }
                            }
                        }
                        else
                        {
                            for (int k = 0; k < fix_property_nentry[j]; k++)
                                fix_property[j]->array_atom[m][k] = fix_property_value[j][k];
                        }
                    }
                }
        //}
    }

    return inserted;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::create_bonds_implicit()
{
    if(nspheres == 1 || !needs_bonding || local_start < 0)
        return 0;

    needs_bonding = false; // reset in case pti gets reused

    // find bond partners
    int *npartner = new int[nspheres](); // convert to member variable to be set from fix template
    int **partner = new int*[nspheres];  // convert to member variable to be set from fix template

    for(int i = 0; i < nspheres; ++i)
        partner[i] = new int[nspheres-1]();

    int create_bonds = 0;

    for(int i = 0; i < nspheres-1; ++i)
    {
        const double xtmp = x_ins[i][0];
        const double ytmp = x_ins[i][1];
        const double ztmp = x_ins[i][2];

        for(int j = i+1; j < nspheres; ++j)
        {
            const double max_bonding_dist = radius_ins[i] + radius_ins[j] + neighbor->skin;
            const double delx = xtmp - x_ins[j][0];
            const double dely = ytmp - x_ins[j][1];
            const double delz = ztmp - x_ins[j][2];
            const double rsq = delx * delx + dely * dely + delz * delz;

            if(rsq < max_bonding_dist*max_bonding_dist)
            {
                if(npartner[i] == nspheres-1 || npartner[j] == nspheres-1)
                {
                    // should not happen print warning
                    continue;
                }
                partner[i][npartner[i]] = j;
                partner[j][npartner[j]] = i;
                npartner[i]++;
                npartner[j]++;
                create_bonds = 1;
            }
        }
    }

    int ncreate = 0;

    if(create_bonds)
    {
        // create bonds
        int **_bond_type = atom->bond_type;
        int **_bond_atom = atom->bond_atom;
        int *num_bond = atom->num_bond;
        int newton_bond = force->newton_bond;
        int n_bondhist = atom->n_bondhist;
        double ***bond_hist = atom->bond_hist;

        // Note: atoms are created with dummy tag = 0, but
        //       actual tags must be available at this point,
        //       i.e. atom->tag_extend() must have been called
        for(int i = 0; i < nspheres; ++i)
        {
            if (npartner[i] == 0) continue;

            for(int k = 0; k < npartner[i]; ++k)
            {
                const int j = partner[i][k];
                if (!newton_bond || i < j)
                {
                    const int ilocal = local_start + i;

                    if (num_bond[ilocal] == atom->bond_per_atom)
                    {
                        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
                    }

                    _bond_type[ilocal][num_bond[ilocal]] = bond_type;
                    _bond_atom[ilocal][num_bond[ilocal]] = atom->tag[local_start+j];

                    // reset history
                    for (int ih = 0; ih < n_bondhist; ++ih)
                    {
                        bond_hist[ilocal][num_bond[ilocal]][ih] = 0.;
                    }
                    num_bond[ilocal]++;
                }

                if(i < j)
                    ++ncreate;
            }
        }
    }

    for(int i = 0; i < nspheres; ++i)
        delete [] partner[i];
    delete [] partner;
    delete [] npartner;

    return ncreate;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::create_bonds_explicit(int *npartner, int **partner)
{
    if(nspheres == 1 || !needs_bonding || local_start < 0)
        return 0;

    needs_bonding = false; // reset in case pti gets reused

    int create_bonds = 1;

    int ncreate = 0;

    if(create_bonds)
    {
        // create bonds
        int **_bond_type = atom->bond_type;
        int **_bond_atom = atom->bond_atom;
        int *num_bond = atom->num_bond;
        int newton_bond = force->newton_bond;
        int n_bondhist = atom->n_bondhist;
        double ***bond_hist = atom->bond_hist;

        // Note: atoms are created with dummy tag = 0, but
        //       actual tags must be available at this point,
        //       i.e. atom->tag_extend() must have been called
        for(int i = 0; i < nspheres; ++i)
        {
            if (npartner[i] == 0) continue;

            for(int k = 0; k < npartner[i]; ++k)
            {
                const int j = partner[i][k];
                if (!newton_bond || i < j)
                {
                    const int ilocal = local_start + i;

                    if (num_bond[ilocal] == atom->bond_per_atom)
                    {
                        error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
                    }

                    _bond_type[ilocal][num_bond[ilocal]] = bond_type;
                    _bond_atom[ilocal][num_bond[ilocal]] = atom->tag[local_start+j];

                    // reset history
                    for (int ih = 0; ih < n_bondhist; ++ih)
                    {
                        bond_hist[ilocal][num_bond[ilocal]][ih] = 0.;
                    }
                    num_bond[ilocal]++;
                }

                if(i < j)
                    ++ncreate;
            }
        }
    }

    return ncreate;
}

/* ---------------------------------------------------------------------- */
//NP checks against xnear list
//NP returns # inerted spheres
//NP if >1, increase nbodies

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    if(nspheres > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,xnear,nnear);

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

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nspheres; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        for(int i = 0; i < nnear; i++)
        {
           vectorSubtract3D(xins_j_try,xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           // no success in overlap
           if (rsq <= radsum*radsum)
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nspheres; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nspheres; j++)
    {
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
    }

    return nspheres;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList)
{
    if(nspheres > 1)
        return check_near_set_x_v_omega_ms(x,v, omega,quat,neighList);

    vectorCopy3D(x,x_ins[0]);

    if(neighList.hasOverlap(x_ins[0], radius_ins[0])) {
        return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    neighList.insert(x_ins[0], radius_ins[0]);

    return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList)
{
    // x is position where insertion should take place
    // v and omega are the velocity and omega for the newly inserted particles
    double rel[3],xins_j_try[3];
    //double del[3], rsq, radsum;

    // check insertion position, take quat into account
    // relative position of spheres to each other already stored at this point
    // check sphere against all others in xnear
    for(int j = 0; j < nspheres; j++)
    {
        // take orientation into account; x_bound_ins is in the global coordinate system
        // calculate xins_j_try for every sphere and check if would work
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,xins_j_try);

        if(neighList.hasOverlap(xins_j_try, radius_ins[j])) {
            return 0;
        }
    }

    // no overlap with any other - success
    // set x_ins, v_ins and omega_ins
    for(int j = 0; j < nspheres; j++)
    {
        vectorSubtract3D(x_ins[j],x_bound_ins,rel);
        MathExtraLiggghts::vec_quat_rotate(rel,quat);
        vectorAdd3D(rel,x,x_ins[j]);
    }
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear for future checks
    for(int j = 0; j < nspheres; j++)
    {
        neighList.insert(x_ins[j], radius_ins[j]);
    }

    return nspheres;
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

    vectorScalarMult3D(x_bound_ins,r_scale);
}
