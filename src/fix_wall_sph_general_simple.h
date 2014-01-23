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

#ifdef FIX_CLASS

FixStyle(wall/sph/general/simple,FixWallSphGeneralSimple)

#else

#ifndef LMP_FIX_WALL_SPH_GENERAL_SIMPLE_H
#define LMP_FIX_WALL_SPH_GENERAL_SIMPLE_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallSphGeneralSimple : public FixWallGran {

   public:
      FixWallSphGeneralSimple(class LAMMPS *, int, char **);
      ~FixWallSphGeneralSimple();

      void compute_force(CollisionData & cdata, double *vwall);

    protected:

      // SPH parameters
      double r0,D;
};

}

#endif
#endif
