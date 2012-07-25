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
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(setSph,SetSph)

#else

#ifndef LMP_SET_SPH_H
#define LMP_SET_SPH_H

#include "pointers.h"

namespace LAMMPS_NS {

class SetSph : protected Pointers {
 public:
  SetSph(class LAMMPS *);
  void command(int, char **);

 private:
  char *id;
  int *select;
  int style,ivalue,count;
  double dvalue;

  int kernel_id;
  char *kernel_style;

  double PI;

  void selection(int);
  void set(int);
};

}

#endif
#endif
