﻿/* ----------------------------------------------------------------------
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
   M.Efe Kinaci (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/cfd/chemistry,FixCfdCouplingChemistry)

#else

#ifndef LMP_FIX_CFD_COUPLING_CHEMISTRY_H
#define LMP_FIX_CFD_COUPLING_CHEMISTRY_H

#include "fix.h"
#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingChemistry : public Fix  {

 public:
  FixCfdCouplingChemistry(class LAMMPS *, int, char **);
  void post_create();
  ~FixCfdCouplingChemistry();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void initial_integrate(int);

 protected:
  int num_species;                              // # of species
  char **species_names_;                        // list of species names
  int iarg_;                                    // int narg_
  char **mod_spec_names_;
  char **diffusant_names_, **molarfraction_names;
  int num_diffusant;

  class FixCfdCoupling*     fix_coupling_;
  class FixPropertyAtom*    fix_tgas_;                 // data pulled from cfdemCoupling - partTemp_
  class FixPropertyAtom*    fix_rhogas_;               // data pulled from cfdemCoupling - partRho_
  class FixPropertyAtom**   fix_masschange_;           // data pushed to cfdemCoupling - changeOfSpeciesMass_
  class FixPropertyAtom*    fix_reactionheat_;         // data pushed to cfdemCoupling - reactionHeat_
  class FixPropertyAtom*    fix_totalmole_;            // data pulled from cfdemCoupling - partN_;
  class FixPropertyAtom**   fix_diffusionCoeff_;
  class FixPropertyAtom*    fix_nufField_;
  class FixPropertyAtom*    fix_partReynolds_;
  class FixPropertyAtom**   fix_molarfraction_;
  class FixPropertyAtom*    fix_partPressure_;
  class FixPropertyAtom*    fix_reactantPerParticle_;
  class FixPropertyAtom*    fix_relLayerRadii_;

 private:
  bool init_layer_radii_,use_Re_, use_reactant_;
};

}

#endif
#endif
