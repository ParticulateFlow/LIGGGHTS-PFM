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

#ifndef LMP_MULTISPHERE_I_H
#define LMP_MULTISPHERE_I_H

/* ---------------------------------------------------------------------- */

inline double Multisphere::max_r_bound()
{
    double max_r_bound = 0.;

    for(int i = 0; i < nbody_; i++)
        max_r_bound = MathExtraLiggghts::max(max_r_bound,r_bound_(i));

    MPI_Max_Scalar(max_r_bound,world);

    return max_r_bound;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::copy_body(int from_local, int to_local)
{
    int tag_from =  id_(from_local);

    customValues_.copyElement(from_local, to_local);

    mapArray_[tag_from] = to_local;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::remove_body(int ilocal)
{
    /*NL*/// fprintf(screen,"deleting body tag %d ilocal %d\n",tag,ilocal);
    mapArray_[id_(ilocal)] = -1;

    //NP do not copy the body if the last body in the array is deleted
    /*if(ilocal < nbody_-1)
        copy_body(nbody_-1,ilocal);*/
    customValues_.deleteElement(ilocal);

    nbody_--;
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::calc_nbody_all()
{
   MPI_Sum_Scalar(nbody_,nbody_all_,world);
}

/* ---------------------------------------------------------------------- */

inline void Multisphere::reset_forces(bool extflag)
{
    /*for(int i = 0; i < nbody_; i++)
    {
        vectorZeroize3D(fcm_(i));
        vectorZeroize3D(torquecm_(i));
        if(extflag) vectorZeroize3D(dragforce_cm_(i));
    }*/

    fcm_.setAll(nbody_,0.);
    torquecm_.setAll(nbody_,0.);
    if(extflag) dragforce_cm_.setAll(nbody_,0.);
}

/* ---------------------------------------------------------------------- */

inline int Multisphere::calc_n_steps(int iatom,int body,double *p_ref,double *normalvec,double *v_normal)
{
    double pos_rel[3],dt,dist_normal;
    int ibody,timestep,n_steps;

    //NP skip if atom not in rigid body
    if(body < 0)
        return -1;

    //NP body ID stored in atom is global
    //NP need to know where stored in my data
    ibody = map(body);

    //NP get data
    timestep = update->ntimestep;
    dt = update->dt;

    //NP error if body not owned by this proc
    if(ibody < 0)
        error->one(FLERR,"Illegal situation in FixMultisphere::calc_n_steps");

    //NP return if step already set
    if(start_step_(ibody) >= 0)
        return (start_step_(ibody) - timestep);

    //NP calculate number of steps
    vectorSubtract3D(p_ref,xcm_(ibody),pos_rel);
    dist_normal = vectorDot3D(pos_rel,normalvec);
    n_steps = static_cast<int>(dist_normal/(vectorMag3D(v_normal)*dt));
    start_step_(ibody) = n_steps + timestep;
    v_integrate_.set(ibody,v_normal);

    /*NL*///fprintf(screen,"set step for body %d to %d\n",ibody,start_step_(ibody));
    /*NL*///error->all(FLERR,"set step");

    return n_steps;
}

#endif
