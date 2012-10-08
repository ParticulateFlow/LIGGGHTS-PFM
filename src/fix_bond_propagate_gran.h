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

FixStyle(bond/propagate/gran,FixBondPropagateGran)

#else

#ifndef LMP_FIX_BOND_PROPAGATE_GRAN_H
#define LMP_FIX_BOND_PROPAGATE_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondPropagateGran : public Fix {
 public:
  FixBondPropagateGran(class LAMMPS *, int, char **);
  ~FixBondPropagateGran();
  int setmask();
  void pre_exchange();
  void write_restart(FILE *);
  void restart(char *);

 private:
  void remove_bond(int ilocal,int ibond);
};

}

#endif
#endif
