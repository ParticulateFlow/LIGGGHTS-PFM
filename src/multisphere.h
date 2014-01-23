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

#ifndef LMP_MULTISPHERE_H
#define LMP_MULTISPHERE_H

#include "pointers.h"
#include "custom_value_tracker.h"
#include "mpi_liggghts.h"
#include "update.h"
#include <vector>

namespace LAMMPS_NS {

  class Multisphere : protected Pointers {

    friend class FixMultisphere;

    public:

      void add_body(int nspheres, double *xcm_ins, double *xcm_to_xbound_ins,
                    double r_bound_ins, double *v_ins, double *omega_ins,
                    double mass_ins, double dens_ins, int atomtype_ins, int type_ins,
                    double *inertia_ins, double *ex_space_ins, double *ey_space_ins, double *ez_space_ins,
                    double **displace_ins, int start_step_ins = -1, double *v_integrate_ins = NULL);

      void grow_arrays_per_body_local(int);
      void grow_arrays_per_body_global(int);

      void remove_body(int ilocal);
      void copy_body(int from_local, int to_local);

      void remap_bodies(int *body);

      void clear_map();
      void generate_map();
      void id_extend_body_extend(int *body);

      void calc_nbody_all();
      bool check_lost_atoms(int *body, double *atom_delflag,double *body_existflag);
      int calc_n_steps(int iatom,int body,double *p_ref,double *normalvec,double *v_normal);

      double max_r_bound();

      void reset_forces(bool extflag);

      void* extract(const char *name, int &, int &);

      double *extract_double_scalar(const char *name);
      double **extract_double_vector(const char *name);

      double extract_ke();
      double extract_rke();

      // inline access functions

      inline int n_body()
      { return nbody_; }

      inline int n_body_all()
      { return nbody_all_; }

      inline int tag_max_body()
      { return mapTagMax_; }

      inline int map(int i)
      { return mapArray_?mapArray_[i]:-1; }

      //NP IDs start at 1

      inline int tag(int i)
      { return id_(i); }

      inline bool has_tag(int _tag)
      { return mapArray_[_tag] == -1 ? false : true;}

      inline int atomtype(int i)
      { return atomtype_(i); }

      inline void xcm(double *x_cm,int i)
      { vectorCopy3D(xcm_(i),x_cm); }

      inline void vcm(double *v_cm,int i)
      { vectorCopy3D(vcm_(i),v_cm); }

      inline void x_bound(double *x_bnd,int i)
      {
        vectorZeroize3D(x_bnd);
        MathExtraLiggghts::local_coosys_to_cartesian(x_bnd,xcm_to_xbound_(i),
                            ex_space_(i),ey_space_(i),ez_space_(i));
        vectorAdd3D(xcm_(i),x_bnd,x_bnd);
      }

      inline double r_bound(int i)
      { return r_bound_(i); }

      inline double mass(int i)
      { return masstotal_(i); }

      inline double density(int i)
      { return density_(i); }

      inline void set_v_body(int ibody,double *vel)
      { vcm_.set(ibody,vel); }

      inline class CustomValueTracker& prop()
      { return customValues_; }

    protected:

      Multisphere(LAMMPS *lmp);
      virtual ~Multisphere();

      // class holding fields
      CustomValueTracker &customValues_;

      // # of local rigid bodies, and global # bodies on all procs
      int nbody_, nbody_all_;

      // global-local lookup
      int mapTagMax_;
      int *mapArray_;

      // ID of rigid body
      //NP IDs start with 1
      ScalarContainer<int> &id_;

      // basic body properties

      // center-of-mass coords, vels, forces, torques of each rigid body
      // extra (external) force on center-of-mass of each
      VectorContainer<double,3> &xcm_;
      VectorContainer<double,3> &vcm_;
      VectorContainer<double,3> &fcm_;
      VectorContainer<double,3> &torquecm_;
      VectorContainer<double,3> &dragforce_cm_;

      // angular momentum of each in space coords
      // angular velocity of each in space coords
      // quaternion of each rigid body
      VectorContainer<double,3> &angmom_;
      VectorContainer<double,3> &omega_;
      VectorContainer<double,4> &quat_;

      // density and total mass of each rigid body
      // 3 principal components of inertia of each
      // principal axes of each in space coords
      ScalarContainer<int> &atomtype_;
      ScalarContainer<int> &type_;
      ScalarContainer<double> &density_;
      ScalarContainer<double> &masstotal_;
      VectorContainer<double,3> &inertia_;
      VectorContainer<double,3> &ex_space_;
      VectorContainer<double,3> &ey_space_;
      VectorContainer<double,3> &ez_space_;

      // # of atoms in each rigid body
      ScalarContainer<int> &nrigid_;

      // image flags of xcm of each rigid body
      ScalarContainer<int> &imagebody_;
      VectorContainer<int,4> &remapflag_;

      // flag for on/off of center-of-mass force, torque
      VectorContainer<bool,3> &fflag_;
      VectorContainer<bool,3> &tflag_;

      //NP step to start from for integration
      ScalarContainer<int> &start_step_;
      VectorContainer<double,3> &v_integrate_;

      // bounding radius for each body
      // vector from xcm to center of bound sphere
      ScalarContainer<double> &r_bound_;
      VectorContainer<double,3> &xcm_to_xbound_;
  };

  // *************************************
  #include "multisphere_I.h"
  // *************************************

} //Namespace

#endif
