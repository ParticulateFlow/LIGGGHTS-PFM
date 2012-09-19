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

PairStyle(gran/hertz/history/hamaker,PairGranHertzHistoryHamaker)

#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_HAMAKER_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_HAMAKER_H

#include "pair_gran_hertz_history.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryHamaker : public PairGranHertzHistory {

 friend class FixWallGranHertzHistory;

 public:
  PairGranHertzHistoryHamaker(class LAMMPS *);
  ~PairGranHertzHistoryHamaker();

  virtual void settings(int, char **);
  virtual void init_granular();

  virtual void compute(int, int,int);

 protected:
  void allocate_properties(int);

  class FixPropertyGlobal* aHamaker_; // Hamaker Constant
  class FixPropertyGlobal* hCut_; // min. distance for cohesion force

  double **aHamakerEff, **hCutEff, **hMaxEff;

  double addCohesionForce(int,int,double);

};

}

#endif
#endif
