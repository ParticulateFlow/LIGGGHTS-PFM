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

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(couple/lb/onetoone,FixLbCouplingOnetoone)

#else

#ifndef FIX_LB_COUPLING_ONETOONE
#define FIX_LB_COUPLING_ONETOONE

#include "fix.h"

namespace LAMMPS_NS {
  class FixLbCouplingOnetoone : public Fix {
  public:
    FixLbCouplingOnetoone(class LAMMPS * lmp, int narg, char ** arg);
    ~FixLbCouplingOnetoone();
    
    virtual int setmask();

    virtual void post_create();
    virtual void init();

    virtual void post_force(int);

    double **get_force_ptr();
    double **get_torque_ptr();
    void comm_force_torque();
  private:
    class FixPropertyAtom* fix_dragforce_;
    class FixPropertyAtom* fix_hdtorque_; // hdtorque = hydrodynamic torque
  }; /* class FixLbCouplingOnetoone */

}; /* LAMMPS_NS */

#endif /* FIX_LB_COUPLING_ONETOONE */
#endif /* FIX_CLASS */




