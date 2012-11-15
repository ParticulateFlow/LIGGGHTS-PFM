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

FixStyle(mesh/surface/stress/servo,FixMeshSurfaceStressServo)

#else

#ifndef LMP_FIX_MESH_SURFACE_STRESS_SERVO_H
#define LMP_FIX_MESH_SURFACE_STRESS_SERVO_H

#include "fix.h"
#include "input.h"
#include "math.h"
#include "fix_mesh_surface_stress.h"

#include "modified_andrew.h"
#include <vector>

using MODIFIED_ANDREW_AUX::Circle;
using MODIFIED_ANDREW_AUX::Line;
using MODIFIED_ANDREW_AUX::Point;

namespace LAMMPS_NS {

class FixMeshSurfaceStressServo : public FixMeshSurfaceStress {

    public:

      FixMeshSurfaceStressServo(class LAMMPS *, int, char **);
      virtual ~FixMeshSurfaceStressServo();

      virtual void post_create();

      void init();
      int setmask();

      virtual void setup_pre_force(int vflag);
      void initial_integrate(int vflag);
      void add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall);
      void final_integrate();

      void reset_dt();
      double compute_vector(int n);

    private:

      void init_defaults();
      void error_checks();

      void limit_vel();
      void update_mass();
      void set_v_node();

      int dim_;

      // properties of mesh

      VectorContainer<double,3> &xcm_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &xcm_orig_;

      // servo settings

      double vel_max_, f_servo_;
      double f_servo_vec_[3];
      char *fstr_;
      int fvar_, fstyle_;

      // timesteps and flags for integration

      double dtf_,dtv_;
      VectorContainer<bool,3> &fflag_;

      bool int_flag_;
      int modify_param(int, char **);

      // controller
      double kp_,ki_,kd_;
      double sign_servo_vec_[3],error_vec_[3],sum_error_vec_[3],old_f_total_[3];

      // velocity for each node

      MultiVectorContainer<double,3,3> &v_;

      //NP TODO: Is somewhere a better place for this function?!
      // signum function
      template <typename T> int sgn(T val) {
          return (T(0) < val) - (val < T(0));
      }

      // for area calculation
      // vector of touching particles
      vector<Circle> contacts_;
      ModifiedAndrew *mod_andrew_;

}; //end class

}

#endif
#endif
