/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/pack,FixInsertPack)

#else

#ifndef LMP_FIX_INSERT_PACK_H
#define LMP_FIX_INSERT_PACK_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertPack : public FixInsert {
 public:

  FixInsertPack(class LAMMPS *, int, char **);
  ~FixInsertPack();

  void init();
  virtual void restart(char *);

 protected:

  virtual void calc_insertion_properties();
  void init_defaults();

  void calc_region_volume_local();

  virtual int calc_ninsert_this();
  virtual int calc_maxtry(int);
  void x_v_omega(int,int&,int&,double&);
  double insertion_fraction();

  int is_nearby(int);
  int is_nearby_body(int);
  virtual BoundingBox getBoundingBox() const;

  // region to be used for insertion
  class Region *ins_region;
  char *idregion;
  double region_volume,region_volume_local;
  int ntry_mc;
  double px_;
  double py_;
  double pz_;
  bool insert_at;

  // target that region should fulfil after each insertion
  double volumefraction_region;
  int ntotal_region;
  double masstotal_region;

  // ratio how many particles have been inserted
  double insertion_ratio;

  // warn if region extends outside box
  bool warn_region;

};

}

#endif
#endif
