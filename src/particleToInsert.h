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

#ifndef LMP_PARTICLE_TO_INSERT_H
#define LMP_PARTICLE_TO_INSERT_H

#include "memory.h"
#include "pointers.h"
#include "region_neighbor_list.h"

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsert : protected Pointers
    {
     public:

        ParticleToInsert(LAMMPS* lmp,int ns = 1);

        //NP need to make this virtual, so delete invoked by FixTemplateSphere
        //NP destroys the right type(ParticleToInsert or ParticleToInsertMultisphere)
        virtual ~ParticleToInsert();

        // insertion properties
        int nspheres;
        int groupbit;
        int atom_type;
        double density_ins;
        double volume_ins;
        double mass_ins;
        double r_bound_ins;

        // per-sphere radius, position
        // if atom_type_vector exists, each sphere has different type
        double *radius_ins;
        double **x_ins;
        bool atom_type_vector_flag;
        int *atom_type_vector;

        // center of bounding sphere

        double x_bound_ins[3];

        // velocity and omega at insertion
        //NP is per-body
        //NP per-sphere vel and omega calculated from
        //NP   rigid body constraint afterwards
        double v_ins[3];
        double omega_ins[3];

        // value of a fix property/atoms at insertion
        class FixPropertyAtom **fix_property;
        int n_fix_property;
        int *fix_property_nentry;
        double **fix_property_value;

        virtual int insert();
        virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        virtual int check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList);
        virtual int check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear);
        virtual int check_near_set_x_v_omega_ms(double *x,double *v, double *omega, double *quat, LIGGGHTS::RegionNeighborList & neighList);

        virtual int set_x_v_omega(double *,double *,double *,double *);

        virtual void scale_pti(double r_scale);
    };

}

#endif
