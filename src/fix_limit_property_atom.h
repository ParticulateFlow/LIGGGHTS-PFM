/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

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

FixStyle(limit/property/atom,FixLimitPropertyAtom)

#else

#ifndef LMP_FIX_LIMIT_PROPERTY_ATOM_H
#define LMP_FIX_LIMIT_PROPERTY_ATOM_H

#include "fix.h"
#include "input.h"
#include "fix_property_atom.h"

namespace LAMMPS_NS {

class FixLimitPropertyAtom : public Fix {
 friend class Set;
 public:
  FixLimitPropertyAtom(class LAMMPS *, int, char **,bool parse = true);
  ~FixLimitPropertyAtom();
  void init();
  virtual int setmask();
  virtual void end_of_step();

 protected:
  virtual void parse_args(int narg, char **arg);

  char *targetfixname;   // name of the fix property atom to be limited
  class FixPropertyAtom* fix_target_;

  int nvalues;
  double *maxvalues;
  double *minvalues;
}; //end class

}
#endif
#endif
