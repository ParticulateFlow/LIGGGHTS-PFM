/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2017-     JKU Linz

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

FixStyle(massflow/mesh/face/universe,FixMassflowMeshFaceUniverse)

#else

#ifndef LMP_FIX_MASSFLOW_MESH_FACE_UNIVERSE_H
#define LMP_FIX_MASSFLOW_MESH_FACE_UNIVERSE_H

#include "fix_massflow_mesh_face.h"

namespace LAMMPS_NS {

class FixMassflowMeshFaceUniverse : public FixMassflowMeshFace {

 public:

  FixMassflowMeshFaceUniverse(class LAMMPS *lmp, int narg, char ** arg);
  ~FixMassflowMeshFaceUniverse();

 protected:
  int couple_every_;
  unsigned int id_hash_;
  int send_to_world_;
  bigint next_couple_;

 private:
  std::map<ConstantParticleTemplateSphere, std::vector<double> > distributions_faces_;
  virtual void reset_distributions(int);
  virtual void increment_distribution(const ConstantParticleTemplateSphere&, int);
  virtual void send_post_create_data();
  virtual void send_coupling_data();
}; //end class

}
#endif
#endif
