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

FixStyle(particletemplate/multisphere,FixTemplateMultisphere)

#else

#ifndef LMP_FIX_TEMPLATE_MULTISPHERE_H
#define LMP_FIX_TEMPLATE_MULTISPHERE_H

#include "fix_template_multiplespheres.h"

namespace LAMMPS_NS {

class FixTemplateMultisphere : public FixTemplateMultiplespheres {
 public:
  FixTemplateMultisphere(class LAMMPS *, int, char **);
  ~FixTemplateMultisphere();
  virtual void post_create();

  void init();

  // called at single insertion
  virtual void randomize_single();

  // called at multi insertion
  void init_ptilist(int);
  void delete_ptilist();
  void randomize_ptilist(int ,int );

  void finalize_insertion();

  int type() {return type_;}

  //NP seed, random generator: inherited

 protected:

  void calc_volumeweight();
  void calc_inertia();
  void calc_eigensystem();
  void calc_displace_xcm_x_body();
  void print_info();

  // type of clump
  int type_;

  // flags
  bool mass_set_, moi_set_;
  int use_density_;

  // inertia of clump
  double moi_[3][3];      // 3x3 inertia tensor
  double inertia_[3]; // 3 principal components of inertia
  double ex_space_[3],ey_space_[3],ez_space_[3]; // principal axes in global coords

  // sphere coordinates in ex,ey,ez frame
  double **displace_;

  // vector from center of mass (which is 0 0 0) to x_bound in body coordinates
  double xcm_to_xb_body_[3];

  // volume weight of each sphere
  // used for volume fraction calculation
  // 1 for spherical or non-overlapping multisphere
  // < 1 for overlapping multisphere
  double *volumeweight_;
};

}

#endif
#endif
