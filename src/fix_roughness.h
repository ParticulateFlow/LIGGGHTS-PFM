/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Stefan Radl, radl@tugraz.at
   Copyright 2013-  Graz University of Technology (TU Graz)

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.

   Description
   Basic fix to account for particle roughness in
   particle-particle, as well as particle-wall collisions.
------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(roughness,FixRoughness)

#else

#ifndef LMP_FIX_ROUGHNESS_H
#define LMP_FIX_ROUGHNESS_H

#include "fix.h"

#define SMALL 1e-8

namespace LAMMPS_NS {

  class FixRoughness : public Fix {
  public:
    FixRoughness(class LAMMPS *, int, char **);
    ~FixRoughness();
    void post_create();
    void pre_delete(bool unfixflag);

    double compute_scalar();
    int setmask();
    void init();
    void initial_integrate(int);
    void post_force(int);

    void register_compute_pair_local(class ComputePairGranLocal *);
    void unregister_compute_pair_local(class ComputePairGranLocal *);
    void updatePtrs();

    double generateDeltaGamma(double gamma);
    double generatePsi();

    inline bool    haveNonZeroDeltaGamma(){return haveNonZeroDeltaGamma_;}

    int     n_history_extra() const;
    bool   history_args(char** args) const;

  protected:
    template <int> void post_force_eval(int,int);
    class ComputePairGranLocal *cpl;

    //Arrays Holding particle information
    double *liqFlux;
    double *liqSource;
    double *liqOnParticle;
    double **condLiqFlux;

    double transfer_ratio;

    //Mean particle roughness and deltaGammaMean
    bool    haveNonZeroDeltaGamma_;
    double roughMean;
    double deltaGammaMean;

    bool FHG_init_flag;

    class PairGran *pair_gran;
    int dnum,dnum_mine;
    int history_flag;

    int seed;
    class RanPark *random_equal;

  };

}

#endif
#endif
