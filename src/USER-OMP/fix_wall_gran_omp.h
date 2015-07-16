/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-2014 DCS Computing GmbH, Linz
                     Copyright 2013-     JKU Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
FixStyle(wall/gran/omp,FixWallGranOMP)
#else

#ifndef LMP_FIX_WALL_GRAN_OMP_H
#define LMP_FIX_WALL_GRAN_OMP_H

#include "fix_wall_gran.h"
#include <map>

namespace LCM = LIGGGHTS::ContactModels;

namespace LAMMPS_NS {

class FixWallGranOMP : public FixWallGran {
 public:
  FixWallGranOMP(class LAMMPS *, int, char **);
  ~FixWallGranOMP();

  /* INHERITED FROM Fix */

  virtual void setup(int vflag);
  virtual void post_force(int vflag);
  virtual void pre_neighbor();
  virtual void pre_force(int vflag);

  /* INHERITED FROM FixWallGran */

  virtual void post_force_wall(int vflag);
  virtual void add_contactforce_wall(int ip, const LCM::ForceData & i_forces);
  virtual void cwl_add_wall_2(const LCM::CollisionData & cdata, const LCM::ForceData & i_forces);

  // mesh and primitive force implementations
  virtual void post_force_mesh(int);
  virtual void post_force_primitive(int);

 private:
  std::map<int, size_t> num_triangles_processed;
  std::map<int, size_t> num_neighbors_processed;
  std::map<int, size_t> num_contacts_computed;
  std::map<int, double> total_compute_time;
  std::map<int, double> contact_compute_time;
  std::map<int, double> idle_time;

  void addHeatFlux(class TriMesh *mesh, int i, int wall_type, double rsq, double area_ratio);

  inline void post_force_eval_contact(LCM::CollisionData & cdata, double * v_wall, int iMesh = -1, FixMeshSurface *fix_mesh = 0, TriMesh *mesh = 0, int iTri = 0);
};

}

#endif
#endif
