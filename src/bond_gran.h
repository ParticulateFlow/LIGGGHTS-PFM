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

#ifdef BOND_CLASS

BondStyle(gran,BondGran)

#else

#ifndef LMP_BOND_GRAN_H
#define LMP_BOND_GRAN_H

#include <stdio.h>
#include "bond.h"

namespace LAMMPS_NS {

class BondGran : public Bond {
 public:
  BondGran(class LAMMPS *);
  ~BondGran();
  void init_style();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *fp);
  double single(int, double, int, int,double&);

 protected:
  int breakmode;
  double *rb,*Sn,*St;
  double *r_break,*sigman_break,*tau_break,*T_break;
  void allocate();

  class FixPropertyAtom *fix_Temp;
  double *Temp;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
