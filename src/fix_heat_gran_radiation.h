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

FixStyle(heat/gran/radiation,FixHeatGranRad)

#else

#ifndef LMP_FIX_HEATGRAN_RADIATION_H
#define LMP_FIX_HEATGRAN_RADIATION_H

#include "fix_heat_gran.h"
#include "random_mars.h"

namespace LAMMPS_NS {

  class FixHeatGranRad : public FixHeatGran {

  public:
    FixHeatGranRad(class LAMMPS *, int, char **);
    ~FixHeatGranRad();
    void pre_delete(bool){}

    double extend_cut_ghost() const;
    int setmask();
    void init();
    void post_force(int);
    void setup(int);

  private:
    bool intersectRaySphere(const double *, const double *, const double *, double, double &, double *);
    int nextBin(int, const double *, const double *, double *, int &, int &, int &);
    int trace(int, int, const double *, const double *, double *, double *);
    void randDir(const double *, double *);
    void randOnSphere(const double *, double, double *, double *);
    void reflect(int, int, int, const double *, const double *, double, double, int, double *);
    void updateQr();
    void createStencils();

    // model parameters
    double Qr;      // energy of one ray
    int avgNRays;   // average number of rays per particle per timestep
    int maxBounces; // maximum number of bounces
    int nTimesteps, updateCounter;
    double cutGhost, cutGhostsq;
    int seed;

    // physical parameters
    double TB; // background temperature

    class FixPropertyGlobal* fix_emissivity;
    double *emissivity;

    double Sigma;        // stefan bolzmann constant
    double Qtot;         // total radiative energy in the system
    class RanMars *RGen; // random number generator

    // pre-allocate these for optimization
    double *raypoint;

    int *stencilLength;
    int *binStencildx;
    int *binStencilmdx;
    int *binStencildy;
    int *binStencilmdy;
    int *binStencildz;
    int *binStencilmdz;
  };
}

#endif
#endif
