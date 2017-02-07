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

FixStyle(insert/pack/face/universe,FixInsertPackFaceUniverse)

#else

#ifndef LMP_FIX_INSERT_PACK_FACE_UNIVERSE_H
#define LMP_FIX_INSERT_PACK_FACE_UNIVERSE_H

#include "fix_insert_pack_face.h"
#include <vector>

namespace LAMMPS_NS {

class FixInsertPackFaceUniverse : public FixInsertPackFace {
 public:

  FixInsertPackFaceUniverse(class LAMMPS *, int, char **);
  ~FixInsertPackFaceUniverse();

 protected:

  void x_v_omega(int,int&,int&,double&);
  double insertion_fraction();

  unsigned int idmassflowface_hash;
  std::map<int, int> faceid2index_;

  int receive_from_world_;
  std::vector<double> nparticles_face_;
  std::vector<double> average_vx_face_out_;
  std::vector<double> average_vy_face_out_;
  std::vector<double> average_vz_face_out_;

private:

  virtual int distribute_ninsert_this(int);
  virtual void post_create_per_face_data();
  virtual void receive_ninsert_this(int&);
  virtual int get_face_index_check(int);
};

}

#endif
#endif
