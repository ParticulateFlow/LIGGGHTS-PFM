/* ----------------------------------------------------------------------
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

#ifndef LMP_NEIGHBOR_H
#define LMP_NEIGHBOR_H

#include "pointers.h"
#include <algorithm> //NP modified C.K.

namespace LIGGGHTS { class RegionNeighborList; }

namespace LAMMPS_NS {

class Neighbor : protected Pointers {
  friend class Cuda;
  //NP modified C.K.
  friend class FixNeighlistMesh;
  friend class FixNeighlistMeshOMP;
  friend class OneLevelGrid;
  /*NL*/ friend class Lbalance;
  //NP modified St.A.
  friend class FixHeatGranRad;
  friend class FixLiquidTrackingInstant;
  friend class LIGGGHTS::RegionNeighborList;

 public:
  int style;                       // 0,1,2 = nsq, bin, multi
  int every;                       // build every this many steps
  int delay;                       // delay build for this many steps
  double contactDistanceFactor;    // contact distance factor to be used for non-contact touch (no radius overlap)
  int dist_check;                  // 0 = always build, 1 = only if 1/2 dist
  int ago;                         // how many steps ago neighboring occurred
  int pgsize;                      // size of neighbor page
  int oneatom;                     // max # of neighbors for one atom
  int includegroup;                // only build pairwise lists for this group
  int build_once;                  // 1 if only build lists once per run
  int cudable;                     // GPU <-> CPU communication flag for CUDA

  double skin;                     // skin distance
  double cutneighmin;              // min neighbor cutoff for all type pairs
  double cutneighmax;              // max neighbor cutoff for all type pairs
                                   //NP this is radius*contactDistanceFactor + skin
                                   //NP contactDistanceFactor is added via PairGran::init_one()
  double *cuttype;                 // for each type, max neigh cut w/ others

  bigint ncalls;                   // # of times build has been called
  bigint ndanger;                  // # of dangerous builds
  bigint lastcall;                 // timestep of last neighbor::build() call

  bigint last_setup_bins_timestep;

  int nrequest;                    // requests for pairwise neighbor lists
  class NeighRequest **requests;   // from Pair, Fix, Compute, Command classes
  int maxrequest;

  int old_style;                   // previous run info to avoid
  int old_nrequest;                // re-creation of pairwise neighbor lists
  int old_triclinic;
  int old_pgsize;
  int old_oneatom;
  class NeighRequest **old_requests;

  int nlist;                       // pairwise neighbor lists
  class NeighList **lists;

  int nbondlist;                   // list of bonds to compute
  int **bondlist;                  //NP
  double **bondhistlist;           //NP modified C.K.
  int nanglelist;                  // list of angles to compute
  int **anglelist;
  int ndihedrallist;               // list of dihedrals to compute
  int **dihedrallist;
  int nimproperlist;               // list of impropers to compute
  int **improperlist;

  Neighbor(class LAMMPS *);
  virtual ~Neighbor();
  virtual void init();
  int request(void *);              // another class requests a neighbor list
  void print_lists_of_lists();      // debug print out
  int decide();                     // decide whether to build or not
  virtual int check_distance();     // check max distance moved since last build
  void setup_bins();                // setup bins based on box and cutoff
  virtual void build(int topoflag=1);  // create all neighbor lists (pair,bond)
  virtual void build_topology();    // create all topology neighbor lists
  void build_one(int);              // create a single neighbor list
  void set(int, char **);           // set neighbor style and skin distance
  void modify_params(int, char**);  // modify parameters that control builds
  bigint memory_usage();
  int exclude_setting();
  int neigh_once(){return build_once;} //NP modified C.K.
  int n_neighs(); //NP modified C.K.
  int n_blist() {return nblist;}

  //NP modified C.K.
  void multi_levels(double &, double &, int &);
  int multi_levels();

 protected:
  int me,nprocs;

  int maxatom;                     // size of atom-based NeighList arrays
  int maxbond,maxangle,maxdihedral,maximproper;   // size of bond lists
  int maxwt;                       // max weighting factor applied + 1

  int must_check;                  // 1 if must check other classes to reneigh
  int restart_check;               // 1 if restart enabled, 0 if no
  int fix_check;                   // # of fixes that induce reneigh
  int *fixchecklist;               // which fixes to check

  double **cutneighsq;             // neighbor cutneigh sq for each type pair
  double **cutneighghostsq;        // neighbor cutnsq for each ghost type pair
  double cutneighmaxsq;            // cutneighmax squared
  double *cuttypesq;               // cuttype squared

  double triggersq;                // trigger = build when atom moves this dist
  int cluster_check;               // 1 if check bond/angle/etc satisfies minimg

  double **xhold;                      // atom coords at last neighbor build
  //NP modified C.K.
  double *rhold;                       // atom radii at last neighbor build
  int maxhold;                         // size of xhold array
  int boxcheck;                        // 1 if need to store box size
  double boxlo_hold[3],boxhi_hold[3];  // box size at last neighbor build
  double corners_hold[8][3];           // box corners at last neighbor build

  int nbinx,nbiny,nbinz;           // # of global bins
  int *bins;                       // ptr to next atom in each bin
  int maxbin;                      // size of bins array

  int *binhead;                    // ptr to 1st atom in each bin
  int maxhead;                     // size of binhead array
  class MultiLevelGrid* mlg;       //NP modified C.K.

  int mbins;                       // # of local bins and offset
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  int binsizeflag;                 // user-chosen bin size
  double binsize_user;

  double binsizex,binsizey,binsizez;  // actual bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  int sx,sy,sz,smax;               // bin stencil extents

  int dimension;                   // 2/3 for 2d/3d
  int triclinic;                   // 0 if domain is orthog, 1 if triclinic
  int newton_pair;                 // 0 if newton off, 1 if on for pairwise

  double *bboxlo,*bboxhi;          // ptrs to full domain bounding box
  double (*corners)[3];            // ptr to 8 corners of triclinic box

  double inner[2],middle[2];       // rRESPA cutoffs for extra lists
  double cut_inner_sq;                   // outer cutoff for inner neighbor list
  double cut_middle_sq;            // outer cutoff for middle neighbor list
  double cut_middle_inside_sq;     // inner cutoff for middle neighbor list

  int special_flag[4];             // flags for 1-2, 1-3, 1-4 neighbors

  int anyghostlist;                // 1 if any non-occasional list
                                   // stores neighbors of ghosts

  int exclude;                     // 0 if no type/group exclusions, 1 if yes

  int nex_type;                    // # of entries in type exclusion list
  int maxex_type;                  // max # in type list
  int *ex1_type,*ex2_type;         // pairs of types to exclude
  int **ex_type;                   // 2d array of excluded type pairs

  int nex_group;                   // # of entries in group exclusion list
  int maxex_group;                 // max # in group list
  int *ex1_group,*ex2_group;       // pairs of group #'s to exclude
  int *ex1_bit,*ex2_bit;           // pairs of group bits to exclude

  int nex_mol;                     // # of entries in molecule exclusion list
  int maxex_mol;                   // max # in molecule list
  int *ex_mol_group;               // molecule group #'s to exclude
  int *ex_mol_bit;                 // molecule group bits to exclude

  int no_build;                    // no neigh lists are built, but exchange
                                   // of particles takes place

  int nblist,nglist,nslist;    // # of pairwise neigh lists of various kinds
  int *blist;                  // lists to build every reneighboring
  int *glist;                  // lists to grow atom arrays every reneigh
  int *slist;                  // lists to grow stencil arrays every reneigh

  void bin_atoms();                     // bin all atoms
  double bin_distance(int, int, int);   // distance between binx
  double bin_largest_distance(int, int, int); //NP modified C.K.
  int coord2bin(const double *);              // mapping atom coord to a bin
  int coord2bin(const double *, int &, int &, int&); // ditto
  void bin_center(int ix, int iy, int iz, double * center);

  void binBorders(int, double &, double &, double &, double &, double &, double &); //NP modified St.A.
  void bin2XYZ(int, int &, int &, int &); //NP modified St.A.
  int XYZ2bin(int, int, int); //NP modified St.A.
  int binHop(int, int, int, int); //NP modified St.A.

  int exclusion(int, int, int,
                int, int *, int *) const;  // test for pair exclusion

  virtual void choose_build(int, class NeighRequest *);
  void choose_stencil(int, class NeighRequest *);

  // pairwise build functions

  typedef void (Neighbor::*PairPtr)(class NeighList *);
  PairPtr *pair_build;

  void half_nsq_no_newton(class NeighList *);
  void half_nsq_no_newton_ghost(class NeighList *);
  void half_nsq_newton(class NeighList *);

  void half_bin_no_newton(class NeighList *);
  void half_bin_no_newton_ghost(class NeighList *);
  void half_bin_newton(class NeighList *);
  void half_bin_newton_tri(class NeighList *);

  void half_multi_no_newton(class NeighList *);
  void half_multi_newton(class NeighList *);
  void half_multi_newton_tri(class NeighList *);

  void full_nsq(class NeighList *);
  void full_nsq_ghost(class NeighList *);
  void full_bin(class NeighList *);
  void full_bin_ghost(class NeighList *);
  void full_multi(class NeighList *);

  void half_from_full_no_newton(class NeighList *);
  void half_from_full_newton(class NeighList *);
  void skip_from(class NeighList *);
  void skip_from_granular(class NeighList *);
  void skip_from_respa(class NeighList *);
  void copy_from(class NeighList *);

  void granular_nsq_no_newton(class NeighList *);
  void granular_nsq_newton(class NeighList *);
  void granular_bin_no_newton(class NeighList *);
  void granular_bin_no_newton_ghost(NeighList *list); //NP modified C.K.
  void granular_bin_newton(class NeighList *);
  void granular_bin_newton_tri(class NeighList *);

  void granular_multi_no_newton(class NeighList *); //NP modified C.K.

  void respa_nsq_no_newton(class NeighList *);
  void respa_nsq_newton(class NeighList *);
  void respa_bin_no_newton(class NeighList *);
  void respa_bin_newton(class NeighList *);
  void respa_bin_newton_tri(class NeighList *);

  // include prototypes for multi-threaded neighbor lists
  // builds or their corresponding dummy versions

#define LMP_INSIDE_NEIGHBOR_H
#include "accelerator_omp.h"
#undef LMP_INSIDE_NEIGHBOR_H

  // pairwise stencil creation functions

  typedef void (Neighbor::*StencilPtr)(class NeighList *, int, int, int);
  StencilPtr *stencil_create;

  void stencil_half_bin_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_ghost_bin_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_bin_3d_no_newton(class NeighList *, int, int, int);
  void stencil_half_ghost_bin_3d_no_newton(class NeighList *, int, int, int);
  void stencil_half_bin_2d_newton(class NeighList *, int, int, int);
  void stencil_half_bin_3d_newton(class NeighList *, int, int, int);
  void stencil_half_bin_2d_newton_tri(class NeighList *, int, int, int);
  void stencil_half_bin_3d_newton_tri(class NeighList *, int, int, int);

  void stencil_half_multi_2d_no_newton(class NeighList *, int, int, int);
  void stencil_half_multi_3d_no_newton(class NeighList *, int, int, int);
  void stencil_half_multi_2d_newton(class NeighList *, int, int, int);
  void stencil_half_multi_3d_newton(class NeighList *, int, int, int);
  void stencil_half_multi_2d_newton_tri(class NeighList *, int, int, int);
  void stencil_half_multi_3d_newton_tri(class NeighList *, int, int, int);

  void stencil_full_bin_2d(class NeighList *, int, int, int);
  void stencil_full_ghost_bin_2d(class NeighList *, int, int, int);
  void stencil_full_bin_3d(class NeighList *, int, int, int);
  void stencil_full_ghost_bin_3d(class NeighList *, int, int, int);
  void stencil_full_multi_2d(class NeighList *, int, int, int);
  void stencil_full_multi_3d(class NeighList *, int, int, int);

  void stencil_gran_multi_3d_no_newton(NeighList *,int, int, int); //NP modified C.K.

  // topology build functions

  typedef void (Neighbor::*BondPtr)();   // ptrs to topology build functions

  BondPtr bond_build;                 // ptr to bond list functions
  void bond_all();                    // bond list with all bonds
  void bond_partial();                // exclude certain bonds
  void bond_check();

  BondPtr angle_build;                // ptr to angle list functions
  void angle_all();                   // angle list with all angles
  void angle_partial();               // exclude certain angles
  void angle_check();

  BondPtr dihedral_build;             // ptr to dihedral list functions
  void dihedral_all();                // dihedral list with all dihedrals
  void dihedral_partial();            // exclude certain dihedrals
  void dihedral_check(int, int **);

  BondPtr improper_build;             // ptr to improper list functions
  void improper_all();                // improper list with all impropers
  void improper_partial();            // exclude certain impropers

  // find_special: determine if atom j is in special list of atom i
  // if it is not, return 0
  // if it is and special flag is 0 (both coeffs are 0.0), return -1
  // if it is and special flag is 1 (both coeffs are 1.0), return 0
  // if it is and special flag is 2 (otherwise), return 1,2,3
  //   for which level of neighbor it is (and which coeff it maps to)

  inline int find_special(const int *list, const int *nspecial,
                          const int tag) const {
    const int n1 = nspecial[0];
    const int n2 = nspecial[1];
    const int n3 = nspecial[2];

    for (int i = 0; i < n3; i++) {
      if (list[i] == tag) {
        if (i < n1) {
          if (special_flag[1] == 0) return -1;
          else if (special_flag[1] == 1) return 0;
          else return 1;
        } else if (i < n2) {
          if (special_flag[2] == 0) return -1;
          else if (special_flag[2] == 1) return 0;
          else return 2;
        } else {
          if (special_flag[3] == 0) return -1;
          else if (special_flag[3] == 1) return 0;
          else return 3;
        }
      }
    }
    return 0;
  };

  //NP modified C.K.
  void register_contact_dist_factor(double cdf)
  { contactDistanceFactor = std::max(contactDistanceFactor,cdf); }
};

}

#endif

/* ERROR/WARNING messages:

E: Neighbor delay must be 0 or multiple of every setting

The delay and every parameters set via the neigh_modify command are
inconsistent.  If the delay setting is non-zero, then it must be a
multiple of the every setting.

E: Neighbor page size must be >= 10x the one atom setting

This is required to prevent wasting too much memory.

E: Invalid atom type in neighbor exclusion list

Atom types must range from 1 to Ntypes inclusive.

W: Neighbor exclusions used with KSpace solver may give inconsistent Coulombic energies

This is because excluding specific pair interactions also excludes
them from long-range interactions which may not be the desired effect.
The special_bonds command handles this consistently by insuring
excluded (or weighted) 1-2, 1-3, 1-4 interactions are treated
consistently by both the short-range pair style and the long-range
solver.  This is not done for exclusions of charged atom pairs via the
neigh_modify exclude command.

E: Neighbor include group not allowed with ghost neighbors

This is a current restriction within LAMMPS.

E: Neighbor multi not yet enabled for ghost neighbors

This is a current restriction within LAMMPS.

E: Neighbor multi not yet enabled for granular

Self-explanatory.

E: Neighbor multi not yet enabled for rRESPA

Self-explanatory.

E: Too many local+ghost atoms for neighbor list

The number of nlocal + nghost atoms on a processor
is limited by the size of a 32-bit integer with 2 bits
removed for masking 1-2, 1-3, 1-4 neighbors.

W: Building an occasional neighobr list when atoms may have moved too far

This can cause LAMMPS to crash when the neighbor list is built.
The solution is to check for building the regular neighbor lists
more frequently.

E: Domain too large for neighbor bins

The domain has become extremely large so that neighbor bins cannot be
used.  Most likely, one or more atoms have been blown out of the
simulation box to a great distance.

E: Cannot use neighbor bins - box size << cutoff

Too many neighbor bins will be created.  This typically happens when
the simulation box is very small in some dimension, compared to the
neighbor cutoff.  Use the "nsq" style instead of "bin" style.

E: Too many neighbor bins

This is likely due to an immense simulation box that has blown up
to a large size.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid group ID in neigh_modify command

A group ID used in the neigh_modify command does not exist.

E: Neigh_modify include group != atom_modify first group

Self-explanatory.

E: Neigh_modify exclude molecule requires atom attribute molecule

Self-explanatory.

*/
