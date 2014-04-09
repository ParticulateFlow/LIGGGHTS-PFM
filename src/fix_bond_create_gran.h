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

FixStyle(bond/create/gran,FixBondCreateGran)

#else

#ifndef LMP_FIX_BOND_CREATE_GRAN_H
#define LMP_FIX_BOND_CREATE_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondCreateGran : public Fix {
 public:
  FixBondCreateGran(class LAMMPS *, int, char **);
  ~FixBondCreateGran();
  void post_create();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();
  void post_integrate_respa(int, int);
  int modify_param(int,char**);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  bool already_bonded(int,int);

  int me;
  int iatomtype,jatomtype;
  int btype,seed;
  int imaxbond,jmaxbond;
  int inewtype,jnewtype;
  double cutsq,fraction;

  int createcount,createcounttotal;   // bond formation stats

  int nmax;
  int newperts;
  int *bondcount;        // count of created bonds this atom is part of
  int *npartner;           //# of preferred atoms for this atom to bond to //NP modified C.K.
  int **partner;          // IDs of preferred atoms for this atom to bond to //NP modified C.K.
  double *probability;   // random # to use in decision to form bond

  class RanMars *random;
  class NeighList *list;
  int countflag,commflag;
  int nlevels_respa;
};

}

#endif
#endif
