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

FixStyle(couple/cfd/fluidproperties,FixCfdCouplingFluidproperties)

#else

#ifndef LMP_FIX_CFD_COUPLING_FLUIDPROPERTIES_H
#define LMP_FIX_CFD_COUPLING_FLUIDPROPERTIES_H

#include "fix_cfd_coupling.h"
#include <vector>

namespace LAMMPS_NS {

class FixCfdCouplingFluidproperties : public Fix {

 public:
  FixCfdCouplingFluidproperties(class LAMMPS *, int, char **);
  ~FixCfdCouplingFluidproperties();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();

 protected:
  class FixCfdCoupling* fix_coupling_;
  class FixPropertyAtom* fix_fluid_density_;
  class FixPropertyAtom* fix_fluid_viscosity_;

};

}

#endif
#endif
