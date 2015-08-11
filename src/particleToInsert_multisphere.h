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

#ifndef LMP_PARTICLE_TO_INSERT_MULTISPHERE_H
#define LMP_PARTICLE_TO_INSERT_MULTISPHERE_H

#include "particleToInsert.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsertMultisphere : public ParticleToInsert
    {
        public:

           ParticleToInsertMultisphere(LAMMPS* lmp,int ns);
           virtual ~ParticleToInsertMultisphere();

           // per-particle displace in body coordinates
           double **displace;

           // vector to center of bounding sphere in body coos
           double xcm_to_xbound[3];

           // center of mass, should be 0/0/0
           double xcm_ins[3];

           double quat_ins[4];

           // body principal axes in space coords
           // = coordinate axis for body co sys
           double ex_space[3],ey_space[3],ez_space[3];

           double inertia[3];

           int type_ms;

           int insert();
           int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
           int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList);
           int set_x_v_omega(double *,double *,double *, double *);

           void random_rotate(double,double,double);

           virtual void scale_pti(double r_scale);
    };
}

#endif
