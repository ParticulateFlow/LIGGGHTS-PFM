/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nve/sphere,FixNVESphere)

#else

#ifndef LMP_FIX_NVE_SPHERE_H
#define LMP_FIX_NVE_SPHERE_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixNVESphere : public FixNVE {
 public:
  FixNVESphere(class LAMMPS *, int, char **);
  virtual ~FixNVESphere() {}
  void init();
  void post_create();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  bool implicitIntegration() {return implicitIntegration_;}

 protected:
  int extra;

  bool implicitIntegration_;
  class FixPropertyAtom* fix_Ksl_;
  class FixCfdCouplingForceImplicit* fix_cfd_coupling_force_implicit_;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix nve/sphere requires atom style sphere

Self-explanatory.

E: Fix nve/sphere requires atom attribute mu

An atom style with this attribute is needed.

E: Fix nve/sphere requires extended particles

This fix can only be used for particles of a finite size.

*/
