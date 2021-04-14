/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com


LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov


Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.


See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(bond/gran,AtomVecBondGran)

#else

#ifndef LMP_ATOM_VEC_BOND_GRAN_H
#define LMP_ATOM_VEC_BOND_GRAN_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecBondGran : public AtomVec {
 public:
  AtomVecBondGran(class LAMMPS *);
  void settings(int narg, char **arg);
  void init();
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  int pack_border_one(int, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int unpack_border_one(int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  int data_atom_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  bigint memory_usage();

 private:
  int *tag,*type,*mask,*image;
  double **x,**v,**f;
  double **x_mol,**v_mol;   // modified A.N.
  int *molecule;
  int **nspecial,**special;
  int *num_bond;
  int **bond_type,**bond_atom;
  int num_bondhist;
  double ***bond_hist;

  class FixBondPropagateGran *fbpg;
};

}

#endif
#endif
