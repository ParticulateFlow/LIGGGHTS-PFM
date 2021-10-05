/* -*- c++ -*- ----------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_H
#define LMP_FIX_H

#include "pointers.h"

namespace LAMMPS_NS {

class Fix : protected Pointers {
 public:
  char *id,*style;
  int igroup,groupbit;

  int restart_global;            // 1 if Fix saves global state, 0 if not
  int restart_peratom;           // 1 if Fix saves peratom state, 0 if not
  int restart_file;              // 1 if Fix writes own restart file, 0 if not
  int force_reneighbor;          // 1 if Fix forces reneighboring, 0 if not
  int box_change_size;           // 1 if Fix changes box size, 0 if not
  int box_change_shape;          // 1 if Fix changes box shape, 0 if not
  int box_change_domain;         // 1 if Fix changes proc sub-domains, 0 if not

  bigint next_reneighbor;        // next timestep to force a reneighboring
  int thermo_energy;             // 1 if fix_modify enabled ThEng, 0 if not
  int nevery;                    // how often to call an end_of_step fix
  int rigid_flag;                // 1 if Fix integrates rigid bodies, 0 if not
  int virial_flag;               // 1 if Fix contributes to virial, 0 if not
  int no_change_box;             // 1 if cannot swap ortho <-> triclinic
  int time_integrate;            // 1 if fix performs time integration, 0 if no
  int time_depend;               // 1 if requires continuous timestepping
  int create_attribute;          // 1 if fix stores attributes that need
                                 //      setting when a new atom is created
  int restart_pbc;               // 1 if fix moves atoms (except integrate)
                                 //      so write_restart must remap to PBC
  int wd_header;                 // # of header values fix writes to data file
  int wd_section;                // # of sections fix writes to data file
  int cudable_comm;              // 1 if fix has CUDA-enabled communication

  int rad_mass_vary_flag;        // 1 if particle radius or mass varied by fix //NP modified C.K.
  int just_created;              // 1 if fix was just created
                                 //   ie no run yet since creation //NP modified C.K.

  int scalar_flag;               // 0/1 if compute_scalar() function exists
  int vector_flag;               // 0/1 if compute_vector() function exists
  int array_flag;                // 0/1 if compute_array() function exists
  int size_vector;               // length of global vector
  int size_array_rows;           // rows in global array
  int size_array_cols;           // columns in global array
  int global_freq;               // frequency s/v data is available at

  int peratom_flag;              // 0/1 if per-atom data is stored
  int size_peratom_cols;         // 0 = vector, N = columns in peratom array
  int peratom_freq;              // frequency per-atom data is available at

  int local_flag;                // 0/1 if local data is stored
  int size_local_rows;           // rows in local vector or array
  int size_local_cols;           // 0 = vector, N = columns in local array
  int local_freq;                // frequency local data is available at

  int extscalar;            // 0/1 if global scalar is intensive/extensive
  int extvector;            // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;             // list of 0/1 int/ext for each vec component
  int extarray;             // 0/1 if global array is intensive/extensive

  double *vector_atom;           // computed per-atom vector
  double **array_atom;           // computed per-atom array
  double *vector_local;          // computed local vector
  double **array_local;          // computed local array

  int comm_forward;              // size of forward communication (0 if none)
  int comm_reverse;              // size of reverse communication (0 if none)
  int comm_border;               // size of border communication (0 if none)

  double virial[6];              // accumlated virial
  double **vatom;                // accumulated per-atom virial

  int recent_restart;            //NP modified C.K. similar than below, but includes both per-atom and global restarts

  int restart_reset;             // 1 if restart just re-initialized fix
  unsigned int datamask;
  unsigned int datamask_ext;

  Fix(class LAMMPS *, int, char **);
  virtual ~Fix();
  void modify_params(int, char **);
  virtual void update_fix(int, char **) {}

  virtual int setmask() = 0;

  virtual void post_create_pre_restart() {} //NP modified C.K.
  virtual void post_create() {} //NP modified C.K.
  virtual void pre_delete(bool) {} //NP modified C.K.
  virtual void box_extent(double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi) {
    UNUSED(xlo); UNUSED(xhi); UNUSED(ylo); UNUSED(yhi); UNUSED(zlo); UNUSED(zhi); //NP modified R.B.
  } //NP modified C.K.
  virtual void init() {}
  virtual void init_list(int, class NeighList *) {}
  virtual void setup(int) {}
  virtual void setup_pre_exchange() {}
  virtual void setup_pre_neighbor() {}
  virtual void setup_pre_force(int) {}
  virtual void min_setup(int) {}
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void pre_force(int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  virtual bool iterate_implicitly() {return false;} //NP modified C.K.
  virtual void end_of_step() {}
  virtual void post_run() {}
  virtual void write_restart(FILE *) {}
  virtual void write_restart_file(char *) {}
  virtual void restart(char *) {}

  virtual void grow_arrays(int) {}
  virtual void copy_arrays(int, int, int) {}
  virtual void pre_set_arrays() {}
  virtual void set_arrays(int) {}
  virtual void update_arrays(int, int) {}
  virtual int pack_border(int, int *, double *) {return 0;}
  virtual int unpack_border(int, int, double *) {return 0;}
  virtual int pack_exchange(int, double *) {return 0;}
  virtual int unpack_exchange(int, double *) {return 0;}
  virtual int pack_restart(int, double *) {return 0;}
  virtual void unpack_restart(int, int) {}
  virtual int size_restart(int) {return 0;}
  virtual int maxsize_restart() {return 0;}

  virtual void setup_pre_force_respa(int, int) {}
  virtual void initial_integrate_respa(int, int, int) {}
  virtual void post_integrate_respa(int, int) {}
  virtual void pre_force_respa(int, int, int) {}
  virtual void post_force_respa(int, int, int) {}
  virtual void final_integrate_respa(int, int) {}

  virtual void min_setup_pre_exchange() {}
  virtual void min_setup_pre_neighbor() {}
  virtual void min_setup_pre_force(int) {}
  virtual void min_pre_exchange() {}
  virtual void min_pre_neighbor() {}
  virtual void min_pre_force(int) {}
  virtual void min_post_force(int) {}

  virtual double min_energy(double *) {return 0.0;}
  virtual void min_store() {}
  virtual void min_clearstore() {}
  virtual void min_pushstore() {}
  virtual void min_popstore() {}
  virtual int min_reset_ref() {return 0;}
  virtual void min_step(double, double *) {}
  virtual double max_alpha(double *) {return 0.0;}
  virtual int min_dof() {return 0;}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual double compute_scalar() {return 0.0;}
  virtual double compute_vector(int) {return 0.0;}
  virtual double compute_array(int,int) {return 0.0;}

  virtual int dof(int) {return 0;}
  virtual void deform(int) {}
  virtual void reset_target(double) {}
  virtual void reset_dt() {}
  virtual void reset_timestep(bigint) {}
  virtual void reset_timestep(bigint newstep,bigint oldstep)
  { reset_timestep(newstep); }

  virtual void read_data_header(char *) {}
  virtual void read_data_section(char *, int, char *) {}
  virtual bigint read_data_skip_lines(char *) {return 0;}

  virtual void write_data_header(FILE *, int) {}
  virtual void write_data_section_size(int, int &, int &) {}
  virtual void write_data_section_pack(int, double **) {}
  virtual void write_data_section_keyword(int, FILE *) {}
  virtual void write_data_section(int, FILE *, int, double **, int) {}

  virtual void zero_momentum() {}
  virtual void zero_rotation() {}

  virtual int modify_param(int, char **) {return 0;}
  virtual void *extract(const char *, int &) {return NULL;}

  virtual double memory_usage() {return 0.0;}

  virtual unsigned int data_mask() {return datamask;}
  virtual unsigned int data_mask_ext() {return datamask_ext;}

  virtual double min_rad(int) const {return 0.0;} //NP modified C.K.
  virtual double max_rad(int) const {return 0.0;} //NP modified C.K.
  virtual int min_type() const {return 0;} //NP modified C.K.
  virtual int max_type() const {return 0;} //NP modified C.K.
  virtual double extend_cut_ghost() const {return 0.0;} //NP modified C.K.
  virtual int n_history_extra() const {return 0;} //NP modified C.K.
  virtual bool history_args(char** args) const { UNUSED(args); return false; } //NP modified C.K.

 protected:
  int evflag;
  int vflag_global,vflag_atom;
  int maxvatom;

  void v_setup(int);
  void v_tally(int, int *, double, double *);

 private:
  // add timing functionality to all fixes
  double recorded_time;
  double previous_time;

 public:
  inline void reset_time_recording() {
    recorded_time = 0.0;
  }

  inline double get_recorded_time() const {
    return recorded_time;
  }

  inline void begin_time_recording() {
    previous_time = MPI_Wtime();
  }

  inline void end_time_recording() {
    double delta_time = MPI_Wtime() - previous_time;
    recorded_time += delta_time;
  }

  union ubuf {  //NP modified R.B.
    double   d;
    int64_t  i;
    ubuf(double arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

  int get_vector_atom_int(int index) {  //NP modified R.B.
    return ubuf(vector_atom[index]).i;
  }

  void set_vector_atom_int(int index, int value) {  //NP modified R.B.
    vector_atom[index] = ubuf(value).d;
  }
};

typedef void (Fix::*FixMethod)();
typedef void (Fix::*FixMethodWithVFlag)(int);
typedef void (Fix::*FixMethodRESPA2)(int,int);
typedef void (Fix::*FixMethodRESPA3)(int,int,int);

namespace FixConst {
  static const int INITIAL_INTEGRATE =       1<<0;
  static const int POST_INTEGRATE =          1<<1;
  static const int PRE_EXCHANGE =            1<<2;
  static const int PRE_NEIGHBOR =            1<<3;
  static const int PRE_FORCE =               1<<4;
  static const int POST_FORCE =              1<<5;
  static const int FINAL_INTEGRATE =         1<<6;
  static const int END_OF_STEP =             1<<7;
  static const int THERMO_ENERGY =           1<<8;
  static const int INITIAL_INTEGRATE_RESPA = 1<<9;
  static const int POST_INTEGRATE_RESPA =    1<<10;
  static const int PRE_FORCE_RESPA =         1<<11;
  static const int POST_FORCE_RESPA =        1<<12;
  static const int FINAL_INTEGRATE_RESPA =   1<<13;
  static const int MIN_PRE_EXCHANGE =        1<<14;
  static const int MIN_PRE_NEIGHBOR =        1<<15;
  static const int MIN_PRE_FORCE =           1<<16;
  static const int MIN_POST_FORCE =          1<<17;
  static const int MIN_ENERGY =              1<<18;
  static const int POST_RUN =                1<<19;
  static const int ITERATE_IMPLICITLY =      1<<20; //NP modified C.K.
  static const int FIX_CONST_LAST =          1<<21; //NP modified C.K.
  static const int PARALLEL_OPENMP =         1<<31; //NP modified R.B.
}

}

#endif

/* ERROR/WARNING messages:

E: Fix ID must be alphanumeric or underscore characters

Self-explanatory.

E: Could not find fix group ID

A group ID used in the fix command does not exist.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
