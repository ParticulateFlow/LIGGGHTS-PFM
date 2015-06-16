/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2009-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef PARTITIONER_H
#define PARTITIONER_H

#include "pointers.h"
#include <vector>

namespace LAMMPS_NS {

class Partitioner : protected Pointers {
public:
  enum Result {
    NO_CHANGE,
    NEW_PARTITIONS,
    FAILED
  };

  Partitioner(class LAMMPS * lmp) : Pointers(lmp) {}
  virtual ~Partitioner(){}
  virtual bool is_cost_effective() const = 0;
  virtual Result generate_partitions(int * permute, std::vector<int> & thread_offsets) = 0;
};

}

#endif
