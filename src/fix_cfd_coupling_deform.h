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

FixStyle(couple/cfd/deform,FixCfdCouplingDeform)

#else

#ifndef LMP_FIX_CFD_COUPLING_DEFORM_H
#define LMP_FIX_CFD_COUPLING_DEFORM_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingDeform : public Fix  {
 public:
  FixCfdCouplingDeform(class LAMMPS *, int, char **);
  ~FixCfdCouplingDeform();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  void init();
  void initial_integrate(int);
  void post_force(int);

  double compute_vector(int n);

 protected:

  class FixCfdCoupling* fix_coupling_;
  class FixPropertyAtom* fix_partdeformations_;
  class FixPropertyAtomPolydispParcel* fix_effvolfactors_;

 private:
  bool verbose_;

  int compress_flag_;

  int igroup_fully_deformed_;

  int igroup_fully_deformed_bit_;

  int particles_removed_;

  double mass_removed_;

  double rmin_;

  const double fmax_;

  double* default_effvolfactors_;

  bool monitor_heat_;

  double heat_removed_;

  class FixPropertyAtom *fix_temp_;

  class FixPropertyGlobal* fix_capacity_;

  bool use_latent_heat_;

  double latent_heat_per_mass_;

  double latent_heat_transferred_;

  class FixPropertyAtom *fix_latent_heat_;

  double *latent_heat_;

  void delete_particle(int);
};

}

#endif
#endif
