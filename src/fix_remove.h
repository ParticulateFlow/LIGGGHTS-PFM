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

FixStyle(remove,FixRemove)

#else

#ifndef LMP_FIX_REMOVE_H
#define LMP_FIX_REMOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRemove : public Fix {
 public:

  FixRemove(class LAMMPS *, int, char **);
  ~FixRemove();

  int setmask();
  void init();
  void pre_exchange();
  void write_restart(FILE *fp);
  void restart(char *buf);

 private:

  void delete_particle(int);

  class RanPark *random;

  int type_remove; //NP atom type of particles to remove
  int iregion, *inflag;
  int style;
  double delete_below;
  double rate_remove; //NP in kg per sec
  int nremoved_this_me;
  int seed;

  double mass_removed,mass_to_remove;
  int time_origin;
  double dt;

};

}

#endif
#endif
