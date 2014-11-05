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

FixStyle(mesh/surface/stress/6dof,FixMeshSurfaceStress6DOF)

#else

#ifndef LMP_FIX_MESH_SURFACE_STRESS_6DOF_H
#define LMP_FIX_MESH_SURFACE_STRESS_6DOF_H

#include "fix.h"
#include "input.h"
#include "math.h"
#include "fix_mesh_surface_stress.h"

namespace LAMMPS_NS {

class FixMeshSurfaceStress6DOF : public FixMeshSurfaceStress {

    public:

      FixMeshSurfaceStress6DOF(class LAMMPS *, int, char **);
      virtual ~FixMeshSurfaceStress6DOF();

      virtual void post_create_pre_restart();
      virtual void post_create();

      void init();
      virtual void setup(int vflag);
      virtual void setup_pre_force(int vflag);

      int setmask();
      void initial_integrate(int vflag);
      void final_integrate();

      double compute_vector(int n);

    private:

      void init_defaults();
      void error_checks();

      void init_rotation_props();
      void calc_displace();

      void set_vel();
      void add_gravity();
      void add_suspension_force();
      void rot_flip(); //NP hidden feature

      // properties of 6 dof body

      VectorContainer<double,3> &xcm_;
      VectorContainer<double,4> &quat_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &omega_;
      VectorContainer<double,3> &angmom_;
      ScalarContainer<double> &mass_;
      MultiVectorContainer<double,3,3> &moi_;
      VectorContainer<double,3> &inertia_;
      VectorContainer<double,3> &ex_space_;
      VectorContainer<double,3> &ey_space_;
      VectorContainer<double,3> &ez_space_;

      // original position and orientation

      VectorContainer<double,3> &xcm_orig_;
      VectorContainer<double,4> &quat_orig_;

      // timesteps and flags for integration

      double dtf_,dtv_,dtfm_,dtq_;
      VectorContainer<bool,3> &fflag_;
      VectorContainer<bool,3> &tflag_;

      // displace for each mesh node

      MultiVectorContainer<double,3,3> &displace_;

      // velocity for each node

      MultiVectorContainer<double,3,3> &v_;

      // additional constraint
      bool suspension_flag_;
      double k_t_, c_t_, k_r_, c_r_;

      //NP hidden feaure
      bool rot_flip_flag_;
      double rot_flip_angle_;

}; //end class

}

#endif
#endif
