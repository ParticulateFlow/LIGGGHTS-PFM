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

FixStyle(multisphere,FixMultisphere)

#else

#ifndef LMP_FIX_MULTISPHERE_H
#define LMP_FIX_MULTISPHERE_H

#include "fix.h"
#include "math_extra_liggghts.h"
#include "multisphere_parallel.h"
#include "fix_property_atom.h"
#include "atom.h"
#include "comm.h"

namespace LAMMPS_NS {

enum
{
    MS_COMM_UNDEFINED,
    MS_COMM_FW_BODY,
    MS_COMM_FW_IMAGE_DISPLACE,
    MS_COMM_FW_V_OMEGA,
    MS_COMM_FW_F_TORQUE,
    MS_COMM_REV_X_V_OMEGA,
    MS_COMM_REV_V_OMEGA
};

class FixMultisphere : public Fix
{
     public:

      FixMultisphere(class LAMMPS *, int, char **);
      virtual ~FixMultisphere();

      void post_create();
      int setmask();
      virtual void init();

      virtual void setup(int);
      virtual void setup_pre_neighbor();

      double extend_cut_ghost();

      void initial_integrate(int);
      void final_integrate();
      void calc_force();

      void add_body_finalize();

      double memory_usage();
      void grow_arrays(int);
      void copy_arrays(int i, int j);

      void pre_exchange();
      void pre_neighbor();
      void set_arrays(int);

      // restart
      int pack_restart(int i, double *buf);
      void unpack_restart(int nlocal, int nth);
      int size_restart(int nlocal);
      int maxsize_restart();
      void write_restart(FILE *fp);
      void restart(char *buf);

      // communication
      int pack_exchange(int i, double *buf);
      int unpack_exchange(int nlocal, double *buf);

      void forward_comm();
      void reverse_comm();

      int pack_comm(int, int*, double*, int, int*);
      int pack_comm_body(int, int*, double*, int, int*);
      int pack_comm_image_displace(int, int*, double*, int, int*);
      int pack_comm_v_omega(int, int*, double*, int, int*);
      int pack_comm_f_torque(int, int*, double*, int, int*);
      void unpack_comm(int, int, double*);
      void unpack_comm_body(int, int, double*);
      void unpack_comm_image_displace(int, int, double*);
      void unpack_comm_v_omega(int, int, double*);
      void unpack_comm_f_torque(int, int, double*);

      int pack_reverse_comm(int, int, double*);
      int pack_reverse_comm_x_v_omega(int, int, double*);
      int pack_reverse_comm_v_omega(int, int, double*);
      void unpack_reverse_comm(int, int*, double*);
      void unpack_reverse_comm_x_v_omega(int, int*, double*);
      void unpack_reverse_comm_v_omega(int, int*, double*);

      int dof(int);
      double ** get_dump_ref(int &nb, int &nprop, char* prop);
      double max_r_bound();

      // public inline access

      void *extract(char *name, int &len1, int &len2)
      {
          if (strcmp(name,"body") == 0)
          {
              len1 = atom->tag_max();
              len2 = 1;
              return (void *) body_;
          }
          return multisphere_.extract(name,len1,len2);
      }

      inline class MultisphereParallel& data()
      { return multisphere_;}

      inline int belongs_to(int i)
      { return body_[i]; }

      inline int n_body()
      { return data().n_body(); }

      inline int n_body_all()
      { return data().n_body_all(); }

      inline int tag_max_body()
      { return data().tag_max_body(); }

      inline void set_v_body_from_atom_index(int iatom,double *vel)
      { multisphere_.set_v_body(body_[iatom],vel); }

      inline void set_body_displace(int i,double *_displace,int body_id)
      { body_[i] = body_id; vectorCopy3D(_displace,displace_[i]); }

      int calc_n_steps(int iatom,double *p_ref,double *normalvec,double *v_normal)
      { return multisphere_.calc_n_steps(iatom,body_[iatom],p_ref,normalvec,v_normal); }

     private:

      inline int map(int i)
      { return data().map(i); }

      inline int tag(int i)
      { return data().tag(i); }

      void set_xv();
      void set_xv(int);
      void set_v();
      void set_v(int);

      class MultisphereParallel &multisphere_;
      class FixPropertyAtom *fix_corner_ghost_;
      class FixPropertyAtom *fix_delflag_;
      class FixPropertyAtom *fix_existflag_;
      class FixGravity *fix_gravity_;

      //NP flag stating that image and displace must be communicated to ghosts
      //int comm_di_;
      int fw_comm_flag_;
      int rev_comm_flag_;

      // per-atom properties handled by this fix
      int *body_;                // which body each atom is part of (-1 if none)
      double **displace_;        // displacement of each atom in body coords

      double dtv,dtf,dtq;
};

}

#endif
#endif
