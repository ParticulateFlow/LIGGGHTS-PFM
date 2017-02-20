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
   contributing authors
   Philippe Seil (JKU)
   ---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/pack/dense,FixInsertPackDense)

#else

#ifndef LMP_FIX_INSERT_PACK_DENSE_H
#define LMP_FIX_INSERT_PACK_DENSE_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertPackDense : public Fix {
 public:

  FixInsertPackDense(class LAMMPS *, int, char **);
  ~FixInsertPackDense();

  // virtual void restart(char *);

  virtual int setmask() { return 0; }
 protected:

  // region to be used for insertion
  class Region *ins_region;
  double region_volume,region_volume_local;

};

}

#endif /* FIX_INSERT_PACK_DENSE */
#endif /* FIX_CLASS */
