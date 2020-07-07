/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations
   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com
   Copyright 2015-     JKU Linz
   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   This software is distributed under the GNU General Public License.
   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(meanfreetime,FixMeanFreeTime)

#else

#ifndef LMP_FIX_MEAN_FREE_TIME_H
#define LMP_FIX_MEAN_FREE_TIME_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixMeanFreeTime : public Fix {

 public:

  FixMeanFreeTime(class LAMMPS *, int, char **);
  ~FixMeanFreeTime();

  void init();
  void init_list(int, class NeighList *);
  int setmask();
  void end_of_step();

  virtual double compute_scalar();

 protected:

  class FixPropertyAtom *fix_meanfreetime_;
  class NeighList *list;

  double check_every_;
  double t_start_;
};

}
#endif
#endif

