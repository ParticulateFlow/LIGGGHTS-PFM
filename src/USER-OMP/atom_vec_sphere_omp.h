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
   Contributing author:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(sphere/omp,AtomVecSphereOMP)

AtomStyle(granular/omp,AtomVecSphereOMP)

#else

#ifndef LMP_ATOM_VEC_SPHERE_OMP_H
#define LMP_ATOM_VEC_SPHERE_OMP_H

#include "atom_vec_sphere.h"
#include "fix_omp.h"

namespace LAMMPS_NS {

class AtomVecSphereOMP : public AtomVecSphere {
public:
  AtomVecSphereOMP(class LAMMPS *);
  ~AtomVecSphereOMP() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  bigint memory_usage();
protected:
  int * thread;
  FixOMP * fixOMP;
};

}

#endif
#endif
