/* ----------------------------------------------------------------------
   LIGGGHTS academic

   Copyright 2005- JKU Linz

   LIGGGHTS academic is based on LIGGGHTS
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations
   Christoph Kloss, christoph.kloss@cfdem.com

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Thomas Lichtenegger (JKU Linz)
   thomas.lichtenegger@jku.at
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(property/atom/cumulativetracer,FixPropertyAtomCumulativeTracer)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_CUMULATIVETRACER_H
#define LMP_FIX_PROPERTY_ATOM_CUMULATIVETRACER_H

#include "fix_property_atom.h"

namespace LAMMPS_NS {


class FixPropertyAtomCumulativeTracer : public FixPropertyAtom {

 public:

  FixPropertyAtomCumulativeTracer(class LAMMPS *, int, char **, bool parse = true);
  ~FixPropertyAtomCumulativeTracer();

  virtual void init();
  virtual int setmask();
  void end_of_step();
  double compute_scalar();

 protected:

  int iarg_;

  char *tracer_name_;
  
  double source_strength_;
  double accumulated_source_strength_;
  double absorbed_strength_;
  double begin_time_;
  double end_time_;
  
  int tot_n_marked_;

  int check_every_;

  // params for region marker
  int iregion_;
  char *idregion_;


}; //end class

}
#endif
#endif
