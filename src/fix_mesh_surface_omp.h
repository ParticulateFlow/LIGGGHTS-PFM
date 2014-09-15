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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(mesh/surface/omp,FixMeshSurfaceOMP)
FixStyle(mesh/surface/planar/omp,FixMeshSurfaceOMP)

#else

#ifndef LMP_FIX_SURFACE_MESH_OMP_H
#define LMP_FIX_SURFACE_MESH_OMP_H

#include "fix_mesh.h"
#include "tri_mesh.h"
#include "fix_contact_history.h"
#include "fix_neighlist_mesh_omp.h"
#include "custom_value_tracker.h"
#include "fix_mesh_surface.h"

namespace LAMMPS_NS
{
  class FixMeshSurfaceOMP : public FixMeshSurface
  {
      public:

        FixMeshSurfaceOMP(LAMMPS *lmp, int narg, char **arg);
        virtual ~FixMeshSurfaceOMP();

        virtual void createWallNeighList(int igrp);
        virtual class FixNeighlistMesh* createOtherNeighList(int igrp,const char *nId);
  };

} /* namespace LAMMPS_NS */

#endif /* LMP_FIX_MESH_SURFACE_OMPH */
#endif /* FIX_CLASS */
