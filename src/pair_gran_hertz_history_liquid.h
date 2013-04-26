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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gran/hertz/history/liquid,PairGranHertzHistoryLiquid)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_LIQUID_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_LIQUID_H

#include "pair_gran_hertz_history.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryLiquid : public PairGranHertzHistory {

 friend class FixWallGranHertzHistory;
 friend class FixCheckTimestepGran;

 public:

  PairGranHertzHistoryLiquid(class LAMMPS *);
  //~PairGranHertzHistoryLiquid();

  virtual void settings(int, char **);
  virtual void init_granular();

  virtual void compute_force(int, int,int);


 protected:

  virtual void history_args(char**);

  class FixPropertyGlobal* liquidVolume1; //NP Modified D.N.
  class FixPropertyGlobal* surfaceTension1; //NP Modified D.N.

  double liquidVolume,surfaceTension;

  virtual void addCapillaryForce(double,double,double &,double &, double, double, double);      //NP modified D.N.

  int capillaryflag;

};

}

#endif
#endif
