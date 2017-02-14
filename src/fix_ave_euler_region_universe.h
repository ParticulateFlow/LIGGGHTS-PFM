/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2017- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(ave/euler/region/universe,FixAveEulerRegionUniverse)

#else

#ifndef LMP_FIX_AVE_EULER_REGION_UNIVERSE_H
#define LMP_FIX_AVE_EULER_REGION_UNIVERSE_H

#include "fix_ave_euler_region.h"

namespace LAMMPS_NS {

class FixAveEulerRegionUniverse : public FixAveEulerRegion {

 public:

  FixAveEulerRegionUniverse(class LAMMPS *, int, char **);
  ~FixAveEulerRegionUniverse();

 private:
  unsigned int id_hash_;
  int send_to_world_;

  virtual void send_post_create_data();
  virtual void send_coupling_data();
};

}

#endif
#endif
