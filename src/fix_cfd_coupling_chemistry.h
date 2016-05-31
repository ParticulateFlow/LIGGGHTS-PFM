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

FixStyle(couple/cfd/chemistry,FixCfdCouplingChemistry)

#else

#ifndef LMP_FIX_CFD_COUPLING_CHEMISTRY_H
#define LMP_FIX_CFD_COUPLING_CHEMISTRY_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingChemistry : public Fix  {
 public:
  FixCfdCouplingChemistry(class LAMMPS *, int, char **);
  ~FixCfdCouplingChemistry();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void post_force(int);

 protected:
  int num_species_;                          // # of species
  char **species_names_;                      // list of species names
  int iarg_, narg_;

  class FixCfdCoupling* fix_coupling_;
<<<<<<< HEAD
  class FixPropertyAtom* fix_tgas_;         // data pulled from cfdemCoupling - partTemp_
  class FixPropertyAtom* fix_rhogas_;       // data pulled from cfdemCoupling - partRho_
  class FixPropertyAtom** fix_massfrac_;    // data pulled from cfdemCoupling - concentrations_
  class FixPropertyAtom** fix_masschange_;  // data pushed to cfdemCoupling - changeOfSpeciesMass_
  class FixPropertyAtom* fix_reactionheat_; // data pushed to cfdemCoupling - reactionHeat_
=======
//  class FixPropertyAtom* fix_tgas_;
//  class FixPropertyAtom* fix_rhogas_;
//  class FixPropertyAtom** fix_massfrac_;
//  class FixPropertyAtom** fix_masschange_;
//  class FixPropertyAtom* fix_reactionheat_;
>>>>>>> bbd2d2ca3900e8e2a0b6d14d8e32ffe63aaa1475

 private:
  bool use_Re_;

};

}

#endif
#endif
