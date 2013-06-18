
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

FixStyle(diam/max,FixDiamMax)

#else

#ifndef LMP_FIX_DIAM_MAX_H
#define LMP_FIX_DIAM_MAX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDiamMax : public Fix {
 public:

  FixDiamMax(class LAMMPS *, int, char **);
  ~FixDiamMax();

  int setmask();
  void init();
  double compute_scalar();

 private:

  // maximum bounding sphere radius
  double maxrbound_;
};

}

#endif
#endif
