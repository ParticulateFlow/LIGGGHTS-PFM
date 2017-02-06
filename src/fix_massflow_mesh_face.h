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

FixStyle(massflow/mesh/face,FixMassflowMeshFace)

#else

#ifndef LMP_FIX_MASSFLOW_MESH_FACE_H
#define LMP_FIX_MASSFLOW_MESH_FACE_H

#include <vector>
#include <map>
#include <set>
#include "fix.h"
#include "constParticleTemplateSphere.h"

namespace LAMMPS_NS {

class FixMassflowMeshFace : public Fix {

  friend class FixParticledistributionDiscreteFace;
  friend class FixInsertPackFace;
 public:

  FixMassflowMeshFace(class LAMMPS *lmp, int narg, char ** arg);
  ~FixMassflowMeshFace();

  void post_create();
  void pre_delete(bool unfixflag);

  void init();
  void setup(int vflag);
  int setmask();

  void post_integrate();
  void pre_exchange();

  void write_restart(FILE *fp);
  void restart(char *buf);

  double compute_scalar();
  double compute_vector(int index);
  double compute_array(int i, int j);

  double compute_array_by_id(int face_id, int j);

  const std::map<int, int>& get_face_ids() const
  {
    return faceid2index_;
  }

  int get_face_ids_size() const
  {
    return faceid2index_.size();
  }

 protected:

  // true if any given particle is
  // counted only once
  bool once_;

  // in case particles counted should be deleted or transferred
  bool delete_atoms_;
  std::vector<int> atom_tags_delete_;
  double mass_deleted_;
  double nparticles_deleted_;

  class FixPropertyAtom* fix_orientation_;

  const DiscreteParticleDistribution& get_distribution_by_id(int face_id);
  const std::vector<DiscreteParticleDistribution>& get_distributions()
  { return distributions_face_; }

 private:

  bool confirm_classification(int ibody, int iPart, int tri, int current_side);

  class FixMeshSurface *fix_mesh_;
  class FixPropertyAtom *fix_counter_;
  class FixPropertyAtom *fix_prevx_;
  char fixid_[200];
  class FixNeighlistMesh *fix_neighlist_;
  double nvec_[3];
  double pref_[3];

  bool   havePointAtOutlet_;
  bool   insideOut_;
  double pointAtOutlet_[3];

  // mass and particles which was counted
  double mass_;
  int nparticles_;
  double average_vx_out_;
  double average_vy_out_;
  double average_vz_out_;
  std::vector<double> mass_face_;
  std::vector<int> nparticles_face_;
  std::vector<double> average_vx_face_out_;
  std::vector<double> average_vy_face_out_;
  std::vector<double> average_vz_face_out_;

  // additional property to check
  char *property_check_name_;
  class FixPropertyAtom *fix_property_check_;
  double d_property_check_;
  int i_property_check_;
  int property_check_index_;
  bool property_check_int_;
  // additional property to sum
  class FixPropertyAtom *fix_property_;
  double property_sum_;

  // data write
  bool screenflag_;
  FILE *fp_;

  // data for particle and mass flow calculation
  double mass_last_;
  int nparticles_last_;
  std::vector<double> mass_face_last_;
  std::vector<int> nparticles_face_last_;
  double t_count_, delta_t_;
  bool reset_t_count_;

  class FixMultisphere* fix_ms_;
  class MultisphereParallel *ms_;
  class ScalarContainer<int> *ms_counter_;

  std::map<int, int> faceid2index_;
  double cg_, cg3_;

  std::vector<DiscreteParticleDistribution> distributions_face_;

  // containers holding per timestep values
  // using member variables to avoid extensive memory reallocation
  std::vector<std::vector<double> > radius_dist_face_local_this;
  std::vector<std::vector<double> > mass_dist_face_local_this;
  std::vector<std::vector<int> > atomtype_dist_face_local_this;
  std::vector<double> mass_face_this;
  std::vector<int> nparticles_face_this;
  std::vector<double> average_vx_face_out_this;
  std::vector<double> average_vy_face_out_this;
  std::vector<double> average_vz_face_out_this;

  std::map<int, int> classified_particles_this;
  std::map<int, int> crossing_particles_this;
  std::set<int> ignore_this;
  std::set<int> once_this;
  std::vector<bool> defined_this; // are in neighlist of a triangle

  // containers to gather data from all procs
  std::vector<double> radius_dist_this_recv;
  std::vector<double> mass_dist_this_recv;
  std::vector<int> atomtype_dist_this_recv;
  std::vector<int> nparticles_face_this_all;

}; //end class

}
#endif
#endif
