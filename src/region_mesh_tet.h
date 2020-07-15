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

#ifdef REGION_CLASS

RegionStyle(mesh/tet,RegTetMesh)

#else

#ifndef LMP_REGION_TET_MESH_H
#define LMP_REGION_TET_MESH_H

#include "random_park.h"
#include "region.h"
#include "bounding_box.h"
#include <vector>
#include <set>

namespace LAMMPS_NS {

class RegTetMesh : public Region {

 friend class InputMeshTet;

 public:

  RegTetMesh(class LAMMPS *, int, char **);
  ~RegTetMesh();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);
  void rebuild();

  void add_tet(double **n);
  int n_tet();
  double total_vol();
  double tet_vol(int i);
  double tet_acc_vol(int i);

 protected:

   int is_inside_tet(int iTet,double *pos);

   // functions are actually not called at the moment
   virtual void generate_random(double *);
   virtual void generate_random_cut(double *,double);

   void grow_arrays();
   void set_extent();
   double volume_of_tet(double* v0, double* v1, double* v2, double* v3);
   double volume_of_tet(int iTet);

   void mesh_randpos(double *pos);
   int  tet_rand_tri();

   char *filename;
   double scale_fact;
   double off_fact[3], rot_angle[3];

   int nTet,nTetMax;
   double ***node;
   double **center;
   double total_volume;
   double *volume;
   double *acc_volume;
   std::vector<BoundingBox> tet_bbox;

  // icosaedron coordinates for surface_interior and surface_exterior
  double **ico_points;
  static double const phi;
  static int const n_ico_point;
  void precalc_ico_points();

  // search tree stuff

  // this is empirical... turned out that 50 elements per node
  // are a good size because the number of false negatives
  // (tested positions that lie in a bounding box, but not in a
  // tet inside this bounding box) is still rather small.
#define TREE_MIN_ELEMENTS_PER_NODE 50

  typedef std::set<int> TreeBin;

  std::vector<TreeBin> tree_data;
  std::vector<BoundingBox> tree_key;

  int tree_max_depth;

  void build_tree();
  void tree_populate_node(int iTreeNode);
  void extend_bb(BoundingBox &box, TreeBin const &data);

  BoundingBox split_bbox_largest_extent(BoundingBox &orig,bool lower);
  void tree_create_children(int current);

  bool tree_is_inside_bin(double *x, TreeBin const &data);

  bool tree_is_inside(double *x);
  int tree_is_inside(double *x, double r);

  int tree_left(int const i) {return 2*i+1;}
  int tree_right(int const i) {return 2*i+2;}
  int tree_parent(int const i) {return (i-1)/2;}

  bool tree_is_leaf(int const i)
  {
    return tree_left(i) >= tree_size()
      || (tree_data[tree_left(i)].empty() && tree_data[tree_right(i)].empty());
  }

  int tree_size() { return tree_key.size(); }
  int tree_level(int i)
  { int level = 0; while(i>0){ i = i >> 1; level++; } return level; }


#include "region_mesh_tet_I.h"
 private:
  double domain_sublo[3];
  double domain_subhi[3];
};

}

#endif
#endif
