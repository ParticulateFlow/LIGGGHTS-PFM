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

/* ----------------------------------------------------------------------
   Contributing authors:
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */


#ifdef FIX_CLASS

FixStyle(couple/cfd/recurrence,FixCfdCouplingRecurrence)

#else

#ifndef LMP_FIX_CFD_COUPLING_RECURRENCE_H
#define LMP_FIX_CFD_COUPLING_RECURRENCE_H

#include "fix_cfd_coupling.h"

namespace LAMMPS_NS {

class FixCfdCouplingRecurrence : public Fix  {
 public:
  FixCfdCouplingRecurrence(class LAMMPS *, int, char **);
  ~FixCfdCouplingRecurrence();
  void post_create();
  void pre_delete(bool unfixflag);

  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void post_force(int);

 protected:

  void dont_use_force()
  { use_force_ = false; }

  double dragforce_total[3];
  class FixCfdCoupling* fix_coupling_;
  class FixPropertyAtom* fix_vrec_;
  class FixPropertyAtom* fix_vfluc_;    // fluctuations
  class FixPropertyAtom* fix_dragforce_;  // drag due to recurrence fields
  class FixPropertyAtom* fix_volumeweight_;
  class FixPropertyAtom* fix_tracerconcentration_;

 private:
  bool use_force_, use_fluc_, use_dens_, use_type_, use_tracer_;

  bool use_property_;
  char property_name[200];
  char property_type[200];

  bool limit_fluc;
  bool relative_limit;
  double maxvfluc;

  int iregion;
  char *idregion;

  void limit_vfluc(double*, double*, double);
};

}

#endif
#endif
