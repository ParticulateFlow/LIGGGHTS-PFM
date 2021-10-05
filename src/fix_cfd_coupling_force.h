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

FixStyle(couple/cfd/force,FixCfdCouplingForce)
FixStyle(couple/cfd/dragforce,FixCfdCouplingForce)

#else

#ifndef LMP_FIX_CFD_COUPLING_FORCE_H
#define LMP_FIX_CFD_COUPLING_FORCE_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingForce : public Fix  {
 public:
  FixCfdCouplingForce(class LAMMPS *, int, char **);
  ~FixCfdCouplingForce();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void post_force(int);
  double compute_vector(int n);

 protected:

  void dont_use_force()
  { use_force_ = false; }

  double dragforce_total[3];
  double hdtorque_total[3]; // superquadric
  class FixCfdCoupling* fix_coupling_;
  class FixPropertyAtom* fix_dragforce_;
  class FixPropertyAtom* fix_hdtorque_; // hdtorque = hydrodynamic torque
  class FixPropertyAtom* fix_volumeweight_;

 private:
  bool use_force_, use_torque_, use_dens_, use_type_, use_property_, use_molecule_;
  bool apply_force_, apply_torque_;
// superquadric start
  bool use_superquadric_;
// superquadric end

  char property_name[200];
  char property_type[200];
};

}

#endif
#endif
