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

#ifdef COMPUTE_CLASS

ComputeStyle(crosssection,ComputeCrosssection)

#else

#ifndef LMP_COMPUTE_CROSSSECTION_H
#define LMP_COMPUTE_CROSSSECTION_H

#include "compute_contact_atom.h"

namespace LAMMPS_NS {

class ComputeCrosssection : public ComputeContactAtom
{
 public:

  ComputeCrosssection(class LAMMPS *, int, char **);
  ~ComputeCrosssection();

  void compute_vector();
  void compute_peratom();

  void write();

 private:

  void setup_cuts();
  void compute_convex_hull();
  inline int mod(double coo);

  class ModifiedAndrew **mod_andrew_;
  int dim_, n_cuts_;
  double min_,max_;
  double cut_thickness_half_;
  double cut_dist_;

  FILE *file_;
};

}

#endif
#endif
