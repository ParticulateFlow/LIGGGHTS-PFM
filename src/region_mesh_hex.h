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
#include <vtkHexahedron.h>
#include <vtkPoints.h>

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
  int n_hex();
  double total_vol();
  double hex_vol(int i);
  double hex_acc_vol(int i);

 protected:

   int is_inside_hex(int iHex,double *pos);

   // functions are actually not called at the moment
   virtual void generate_random(double *);
   virtual void generate_random_cut(double *,double);

   void grow_arrays();
   void set_extent();
   double volume_of_hex(double* v0, double* v1, double* v2, double* v3,
                        double* v4, double* v5, double* v6, double* v7);
   double volume_of_hex(int iHex);

   void mesh_randpos(double *pos);
   int  mesh_randcell();

   char *filename;
   double scale_fact;
   double off_fact[3], rot_angle[3];

   int nHex,nHexMax;
   double ***node;
   double **center;
   double total_volume;
   double *volume;
   double *acc_volume;

   #include "region_mesh_hex_I.h"
};

}

#endif
#endif
#endif
