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

#include "atom_vec_sphere_omp.h"
#include "atom.h"
#include "comm.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecSphereOMP::AtomVecSphereOMP(LAMMPS *lmp) : AtomVecSphere(lmp)
{
  // register fix omp with this class
  int ifix = lmp->modify->find_fix("package_omp");
  if (ifix < 0)
    lmp->error->all(FLERR,"The 'package omp' command is required for /omp styles");
  fixOMP = static_cast<FixOMP *>(lmp->modify->fix[ifix]);
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSphereOMP::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");

  if(fixOMP->use_reduction()) {
    f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");
    torque = memory->grow(atom->torque,nmax*comm->nthreads,3,"atom:torque");
  } else {
    f = memory->grow(atom->f,nmax,3,"atom:f");
    torque = memory->grow(atom->torque,nmax,3,"atom:torque");
  }

  thread = memory->grow(atom->thread,nmax,"atom:thread");

  radius = memory->grow(atom->radius,nmax,"atom:radius");
  density = memory->grow(atom->density,nmax,"atom:density");
  rmass = memory->grow(atom->rmass,nmax,"atom:rmass");
  omega = memory->grow(atom->omega,nmax,3,"atom:omega");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSphereOMP::grow_reset()
{
  AtomVecSphere::grow_reset();
  thread = atom->thread;
}

/* ----------------------------------------------------------------------
   copy atom I info to atom J
------------------------------------------------------------------------- */

void AtomVecSphereOMP::copy(int i, int j, int delflag)
{
  AtomVecSphere::copy(i, j, delflag);
  thread[j] = thread[i];
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSphereOMP::pack_exchange(int i, double *buf)
{
  int m = AtomVecSphere::pack_exchange(i, buf);
  // atom removed from partition, other from tail will be taking its place
  // enforce resorting of data
  atom->dirty = true;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSphereOMP::unpack_exchange(double *buf)
{
  int m = AtomVecSphere::unpack_exchange(buf);
  int nlocal = atom->nlocal - 1; // nlocal was increased in base class
  thread[nlocal] = -1; // don't know where this should go yet
  atom->dirty = true;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSphereOMP::unpack_restart(double *buf)
{
  int m = AtomVecSphere::unpack_restart(buf);
  int nlocal = atom->nlocal - 1; // nlocal was increased in base class
  thread[nlocal] = -1; // don't know where this will be placed yet
  atom->dirty = true;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecSphereOMP::create_atom(int itype, double *coord)
{
  AtomVecSphere::create_atom(itype, coord);
  int nlocal = atom->nlocal - 1; // nlocal was increased in base class
  thread[nlocal] = -1; // don't know where this will be placed yet
  atom->dirty = true;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSphereOMP::data_atom(double *coord, tagint imagetmp, char **values)
{
  AtomVecSphere::data_atom(coord, imagetmp, values);
  int nlocal = atom->nlocal - 1; // nlocal was increased in base class
  thread[nlocal] = -1; // don't know where this will be placed yet
  atom->dirty = true;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecSphereOMP::memory_usage()
{
  bigint bytes = AtomVecSphere::memory_usage();

  if (atom->memcheck("thread")) bytes += memory->usage(thread,nmax);

  return bytes;
}
