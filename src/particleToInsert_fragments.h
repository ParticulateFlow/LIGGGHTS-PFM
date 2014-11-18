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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_PARTICLE_TO_INSERT_FRAGMENTS_H
#define LMP_PARTICLE_TO_INSERT_FRAGMENTS_H

#include "particleToInsert.h"
#include <string>

using namespace LAMMPS_NS;

namespace LAMMPS_NS {
    class ParticleToInsertFragments : public ParticleToInsert
    {
     public:

        ParticleToInsertFragments(LAMMPS* lmp,int ns = 1);
        virtual ~ParticleToInsertFragments();

        virtual int insert();
        double collision_factor;
        std::string fix_property_atom_id;
    };
}

#endif
