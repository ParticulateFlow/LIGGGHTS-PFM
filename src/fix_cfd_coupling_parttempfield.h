/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particulate Flow Modelling
   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/parttempfield,FixCfdCouplingPartTempField)

#else

#ifndef LMP_FIX_CFD_COUPLING_PARTTEMPFIELD_H
#define LMP_FIX_CFD_COUPLING_PARTTEMPFIELD_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingPartTempField : public Fix {

 public:
  FixCfdCouplingPartTempField(class LAMMPS *, int, char **);
  ~FixCfdCouplingPartTempField();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  class FixCfdCoupling* fix_coupling;
  class FixPropertyAtom* fix_temp;
  bool created_fix_temp;
  double T0;
};

}

#endif
#endif
