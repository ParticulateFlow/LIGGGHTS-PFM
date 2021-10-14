/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2021- Eindhoven University of Technology

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Tim Nijssen (TU/e)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/dissolve,FixCfdCouplingDissolve)

#else

#ifndef LMP_FIX_CFD_COUPLING_DISSOLVE_H
#define LMP_FIX_CFD_COUPLING_DISSOLVE_H

#include "fix_cfd_coupling.h"
#include <vector>

namespace LAMMPS_NS {

class FixCfdCouplingDissolve : public Fix {

 public:
  FixCfdCouplingDissolve(class LAMMPS *, int, char **);
  ~FixCfdCouplingDissolve();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void post_force(int);
  void pre_exchange();

 protected:
  class FixCfdCoupling* fix_coupling_;
  class FixPropertyAtom* fix_convectiveFlux_;
  double rmin_;

  std::vector<int> atom_tags_delete_;

  void delete_atoms();
};

}

#endif
#endif
