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
  void delete_atoms();

 protected:
  class FixCfdCoupling* fix_coupling;
  //class FixPropertyAtom* fix_conductiveFlux;
  class FixPropertyAtom* fix_convectiveFlux;
  //class FixPropertyAtom* fix_heatFlux;
  double rmin;
  //bool gran_field_conduction;

  std::vector<int> atom_tags_delete_;

};

}

#endif
#endif
