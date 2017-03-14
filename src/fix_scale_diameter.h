/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(scale/diameter,FixScaleDiameter)

#else

#ifndef LMP_FIX_SCALE_DIAMETER_H
#define LMP_FIX_SCALE_DIAMETER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixScaleDiameter : public Fix {
 public:

  FixScaleDiameter(class LAMMPS *, int, char **);
  ~FixScaleDiameter();
  int setmask();
  void post_create();
  void pre_delete(bool unfixflag);
  void init();
  void setup_pre_force(int);
  void pre_force(int);

  void set_scale_mass(bool scale) { scale_mass_ = scale; }
  void set_scale(double scale);
  double get_scale() const { return scale_to_; }

 private:

  void change_settings();

  class FixPropertyAtom *fix_property_;
  class Region *scale_region;
  char *idregion;
  int region_style;
  double center_[3];
  double radius_;
  double radius_inner_;
  double scale_to_;
  int scale_to_var_, scale_to_style_;
  char *scale_to_str_;
  double scale_range_;
  double scale_width_;
  bool scale_mass_;
};

}

#endif
#endif
