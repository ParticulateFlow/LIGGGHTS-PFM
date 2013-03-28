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

#ifdef COMPUTE_CLASS

ComputeStyle(rigid,ComputeRigid)

#else

#ifndef LMP_COMPUTE_RIGID_H
#define LMP_COMPUTE_RIGID_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRigid : public Compute {
 public:
  ComputeRigid(class LAMMPS *, int, char **);
  ~ComputeRigid();
  void init();
  void compute_local();
  double memory_usage();

 private:
  void update_pointers();

  class MultisphereParallel* multisphere_;
  class ContainerBase *property_;

};

}

#endif
#endif
