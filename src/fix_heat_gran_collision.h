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
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heat/gran/collision,FixHeatGranColl)

#else

#ifndef LMP_FIX_HEATGRAN_COLLISION_H
#define LMP_FIX_HEATGRAN_COLLISION_H

#include "fix_heat_gran.h"

namespace LAMMPS_NS {

  class FixHeatGranColl : public FixHeatGran {
  public:
    FixHeatGranColl(class LAMMPS *, int, char **);
    ~FixHeatGranColl();
    virtual void post_create();
    void pre_delete(bool);

    int setmask();
    void init();
    virtual void post_force(int);

    void cpl_evaluate(class ComputePairGranLocal *);
    void register_compute_pair_local(ComputePairGranLocal *);
    void unregister_compute_pair_local(ComputePairGranLocal *);

  protected:
    int iarg_;

  private:
    template <int,int> void post_force_eval(int,int);

    class FixPropertyGlobal* fix_conductivity_;
    double *conductivity_;
    
    class FixPropertyGlobal* fix_capacity_;
    double *capacity_;

  };

}

#endif
#endif

