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

PairStyle(gran/hooke/history/cohesion,PairGranHookeHistoryCohesion)

#else

#ifndef LMP_PAIR_GRAN_HOOKE_HISTORY_COHESION_H
#define LMP_PAIR_GRAN_HOOKE_HISTORY_COHESION_H

#include "pair_gran_hooke_history.h"

namespace LAMMPS_NS {

class PairGranHookeHistoryCohesion : public PairGranHookeHistory {

 friend class FixWallGranHookeHistory;
 friend class FixCheckTimestepGran;

 public:

  PairGranHookeHistoryCohesion(class LAMMPS *);
  ~PairGranHookeHistoryCohesion();

  virtual void settings(int, char **);
  virtual void init_granular(); //NP makes inits specific to certain gran style     //NP modified C.K.

  virtual void compute(int, int,int);

 protected:

  virtual void history_args(char**);
  void allocate_properties(int);

  double kn2k2Max_, kn2kc_, phiF_;
};

}

#endif
#endif
