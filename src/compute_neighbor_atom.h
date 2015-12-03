/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-2014 DCS Computing GmbH, Linz
                     Copyright 2015-     JKU Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(neighbor/atom,ComputeNeighborAtom)

#else

#ifndef LMP_COMPUTE_NEIGHBOR_ATOM_H
#define LMP_COMPUTE_NEIGHBOR_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeNeighborAtom : public Compute {
 public:
  ComputeNeighborAtom(class LAMMPS *, int, char **);
  ~ComputeNeighborAtom();
  virtual void init();
  void init_list(int, class NeighList *);
  virtual void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 protected:
  int nmax;
  class NeighList *list;
  double * neighs;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute neighbor/atom requires atom style sphere

Self-explanatory.

E: Compute neighbor/atom requires a pair style be defined

Self-explantory.

W: More than one compute neighbor/atom

It is not efficient to use compute neighbor/atom more than once.

*/
