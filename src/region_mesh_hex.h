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

#if defined(LAMMPS_VTK) //NP do not use #ifdef here (VS C++ bug)
#ifdef REGION_CLASS

RegionStyle(mesh/hex,RegHexMesh)

#else

#ifndef LMP_REGION_HEX_MESH_H
#define LMP_REGION_HEX_MESH_H

#include "random_park.h"
#include "region.h"
#include "custom_value_tracker.h"
#include <vector>
#include <vtkHexahedron.h>
#include <vtkPoints.h>
#include <vtkNew.h>

namespace LAMMPS_NS {

class RegHexMesh : public Region {

  friend class InputMeshHex;

 public:

  RegHexMesh(class LAMMPS *, int, char **);
  ~RegHexMesh();
  int inside(double, double, double);
  int surface_interior(double *, double);
  int surface_exterior(double *, double);

  void add_hex(double **n);
  int n_hex() const;
  double total_vol() const;
  double hex_vol(int i) const;
  double hex_acc_vol(int i) const;
  void hex_bounds(int iHex, double bounds[6]) const;
  void hex_points(int iHex, double points[24]) const;
  const double* hex_center(int iHex) const
  { return center[iHex]; }

  inline CustomValueTracker& prop()
  { return customValues_; }

  int get_hex(double *pos);
  int get_hex(const char* property, int value);
  int match_hex_cut(int iHex,double *pos,double cut);
  int is_inside_hex(int iHex,double *pos);

 protected:


   // functions are actually not called at the moment
   virtual void generate_random(double *);
   virtual void generate_random_cut(double *,double);

   void grow_arrays();
   void set_extent();
   double volume_of_hex(double** v);
   double volume_of_hex(int iHex);

   void mesh_randpos(double *pos);
   int  mesh_randcell();
   bool is_cell_orthogonal(int iHex) const;

   char *filename;
   double scale_fact;
   double off_fact[3], rot_angle[3];
   bool read_cell_data_;

   int nHex,nHexMax;
   double ***node;
   double **center;
   double total_volume;
   double *volume;
   double *acc_volume;
   bool *orthogonal;

   // class holding fields
   CustomValueTracker &customValues_;

   class AABBTree *tree_;
   std::vector<int> potential_cells;
   vtkNew<vtkHexahedron> hexahedron;
   vtkNew<vtkTetra> tetra;
   vtkNew<vtkIdList> ptIds;
   vtkNew<vtkPoints> pts;
   double pCoords[3];
   double weights[8];

 public:
   #include "region_mesh_hex_I.h"
};

}

#endif
#endif
#endif
