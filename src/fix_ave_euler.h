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
   Contributing author for triclinic: Andreas Aigner (JKU)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/euler,FixAveEuler)
FixStyle(ave/euler/stress,FixAveEuler)

#else

#ifndef LMP_FIX_AVE_EULER_H
#define LMP_FIX_AVE_EULER_H

#include <stdio.h>
#include "fix.h"
#include "vector_liggghts.h"

namespace LAMMPS_NS {

class FixAveEuler : public Fix {

 public:

  FixAveEuler(class LAMMPS *, int, char **);
  virtual ~FixAveEuler();

  virtual void post_create();
  virtual int setmask();
  virtual void init();
  void setup(int vflag);
  void post_integrate();
  void pre_force(int);

  void end_of_step();

  virtual double compute_array(int i, int j);

  int ncells_pack() const;

  inline int ncells() const
  { return ncells_; }

  virtual double cell_volume(int) const
  { return cell_volume_; }

  // access functions for cell based values

  virtual double cell_center(int i, int j) const
  { return center_[i][j]; }

  inline double cell_v_av(int i, int j) const
  { return v_av_[i][j]; }

  inline double cell_vol_fr(int i) const
  { return vol_fr_[i]; }

  inline double cell_mass(int i) const
  { return mass_[i]; }

  inline double cell_radius(int i) const
  { return radius_[i]; }

  inline double cell_pressure(int i) const
  { return stress_[i][0]; }

  inline double cell_stress(int i,int j) const
  { return stress_[i][j+1]; }

  inline int cell_count(int i) const
  { return ncount_[i]; }

  inline int cell_head(int i)
  { if (dirty_) lazy_bin_atoms(i); return cellhead_[i]; }

  inline int cell_ptr(int i) const
  { return cellptr_[i]; }

  inline double cell_weight(int i) const
  { return weight_[i]; }

  inline void set_cell_weight(int i, double w)
  { weight_[i] = w; }

 protected:
  inline int ntry_per_cell() const
  { return 50; }

 private:
  virtual void setup_bins();
  virtual void bin_atoms();
  virtual void lazy_bin_atoms(int i) { bin_atoms(); }
  virtual void calculate_eu();
  void allreduce();
  inline int coord2bin(double *x); //NP modified A.A.

 protected:
  bool parallel_;
  bool dirty_;

  int exec_every_;
  bool box_change_size_, box_change_domain_;
  int triclinic_; //NP modified A.A.

  // desired cell size over max particle diameter
  double cell_size_ideal_rel_[3];

  // desired cell size
  double cell_size_ideal_[3];
  double cell_size_ideal_lamda_[3];

  // number of cells, either globally or locally on each proc
  int ncells_, ncells_dim_[3];

  // extent of grid in xyz, either globally or locally on each proc
  double lo_[3],hi_[3];
  double lo_lamda_[3],hi_lamda_[3]; //NP modified A.A.

  // cell size and inverse size in xyz, cell and volume
  double cell_size_[3];
  double cell_size_inv_[3];
  double cell_volume_;
  double cell_size_lamda_[3]; //NP modified A.A.
  double cell_size_lamda_inv_[3]; //NP modified A.A.

  // length of cellhead_, center_, v_av_, vol_fr_ arrays
  int ncells_max_;

  // length of cellptr_ array
  int ncellptr_max_;

  // atom - cell mapping
  int *cellhead_;    // ptr to 1st atom in each cell
  int *cellptr_;       // ptr to next atom in each bin

  // region
  char *idregion_;
  class Region *region_;

  /* ---------  DATA  --------- */

  // cell center
  double **center_;

  // cell-based averaged velocity
  double **v_av_;

  // cell-based volume fraction
  double *vol_fr_;

  // cell-based weight for each cell

  double *weight_;

  // cell-based average radius
  double *radius_;

  // cell-based number of particles
  int *ncount_;

  // cell-based mass
  double *mass_;

  // cell-based stress
  // [0]: pressure
  // [1-3]: 00-11-22
  // [4-6]: 01-02-12
  double **stress_;

  // stress computation
  class ComputeStressAtom *compute_stress_;

  class RanPark *random_;
};

}

#endif
#endif
