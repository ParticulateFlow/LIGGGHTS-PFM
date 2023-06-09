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

FixStyle(couple/cfd/convection,FixCfdCouplingConvection)

#else

#ifndef LMP_FIX_CFD_COUPLING_CONVECTION_H
#define LMP_FIX_CFD_COUPLING_CONVECTION_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingConvection : public Fix {

 public:
  FixCfdCouplingConvection(class LAMMPS *, int, char **);
  ~FixCfdCouplingConvection();
  void post_create();
  void pre_delete(bool unfixflag);

  virtual int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  class FixCfdCoupling* fix_coupling;
  class FixPropertyAtom* fix_conductiveFlux;
  class FixPropertyAtom* fix_convectiveFlux;
  class FixPropertyAtom* fix_heatFlux;
  double T0;
  bool gran_field_conduction;
  
  bool limit_change;
  double max_change;
};

}

#endif
#endif
