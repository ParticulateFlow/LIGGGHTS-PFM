/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/pack/face,FixInsertPackFace)

#else

#ifndef LMP_FIX_INSERT_PACK_FACE_H
#define LMP_FIX_INSERT_PACK_FACE_H

#include "fix_insert.h"
#include <map>
#include <deque>

namespace LAMMPS_NS {

class FixInsertPackFace : public FixInsert {
 public:

  FixInsertPackFace(class LAMMPS *, int, char **);
  ~FixInsertPackFace();

  virtual void post_create();
  virtual void init();
  virtual void restart(char *);

 protected:

  virtual void calc_insertion_properties();
  virtual void init_defaults();

  virtual void calc_region_volume_local();

  virtual int calc_ninsert_this();
  virtual int calc_maxtry(int);
  virtual void x_v_omega(int,int&,int&,double&);
  virtual double insertion_fraction();

  virtual double insertion_fraction_face(int face_id);

  virtual int is_nearby(int);
  virtual BoundingBox getBoundingBox() const;

  // massflow to be used
  class FixMassflowMeshFace *massflowface;
  char *idmassflowface;
  // region to be used for insertion
  class Region *ins_region;
  class RegHexMesh *ins_region_mesh_hex;
  char *idregion;
  double region_volume,region_volume_local;
  int ntry_mc;

  // ratio how many particles have been inserted
  double insertion_ratio;

  // warn if region extends outside box
  bool warn_region;
  double cg_, cg3_;
  int type_offset;
  std::map<int, double> remaining_ptis;
  std::vector<double> fraction_face_local;
  std::vector<double> min_face_extent_local;
  std::vector<double> volume_face_absolut;
  std::deque<double> rounding_all;

 private:

  virtual int distribute_ninsert_this(int);
  virtual void post_create_per_face_data();
  virtual void receive_ninsert_this(int&);
  virtual int get_face_index_check(int);
};

}

#endif
#endif
