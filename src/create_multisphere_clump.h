/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2019- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) // do not use #ifdef here (VS C++ bug)
#ifdef COMMAND_CLASS

CommandStyle(create_multisphere_clump,CreateMultisphereClump)

#else

#ifndef LMP_CREATE_MULTISPHERE_CLUMP_H
#define LMP_CREATE_MULTISPHERE_CLUMP_H

#include "pointers.h"
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vector>

namespace LAMMPS_NS {

  /**
   * @brief CreateMultisphereClump class
   *        create a multi-sphere clump from a vtk mesh file and write the clump data to a file.
   *
   * Uses the vtk library to read an unstructured grid or polydata dataset from vtk simple legacy
   * or xml format and triangulate the polygons. Optionally, subdivide the surface.
   * Finally, use the overlapping multi-sphere clump method (OMCM) to generate a clump configuration
   * and write it to a file.
   * Cf. J.F. Ferellec and G.R. McDowell, Granular Matter (2010) 12:459–467
   * usage:
   * create_multisphere_clump dmin (absolute)|(radius_ratio) 0.0001 rmin 0.002 pmax 1.0 seed 1321
   *                          surfacefile infile.vtk [invert_normals no]
   *                          [subdivide (linear)|(loop)|(butterfly) 3 [subdivisionfile subdiv.vtk]]
   *                          (clumpfile clump.dat [clump.vtk])|(datafile data.clump [atom_type density])
   */
class CreateMultisphereClump : protected Pointers {
 public:
  CreateMultisphereClump(class LAMMPS*);
  ~CreateMultisphereClump();
  void command(int, char**);

 private:
  /**
   * @brief read input mesh
   */
  template<class TReader> vtkDataSet* read_file(const char* filename);
  /**
   * @brief triangulate input mesh
   */
  vtkSmartPointer<vtkPolyData> triangulate(vtkPolyData* pd);
  /**
   * @brief subdivide input mesh
   */
  vtkSmartPointer<vtkPolyData> subdivide(int narg, char** arg, vtkPolyData* pd);
  /**
   * @brief output subdivided mesh
   */
  void write_subdiv_file(int /*narg*/, char** arg, vtkPolyData* pd);
  /**
   * @brief generate the spheres that make up the multisphere clump
   */
  void generate_spheres(vtkPolyData* dset);
  /**
   * @brief output the sphere configuration to a LAMMPS data file
   */
  void write_data_file(const char* filename);
  /**
   * @brief output the sphere configuration to a clump file
   */
  void write_clump_file(const char* filename);
  /**
   * @brief output the sphere configuration to a vtk file
   */
  void write_clump_file_debug(const char*filename);
  /**
   * @brief check the distance between generated spheres and the given point
   */
  bool check_sphere_distance(const double* x);
  /**
   * @brief check if the given sphere contains the given point
   */
  bool is_point_in_sphere(const double* center, double radius, const double* x);
  /**
   * @brief check if the given box contains the given sphere
   */
  bool is_sphere_in_bounds(const double* bounds, const double *center, double radius);

  int me;
  int iarg;
  int atom_type;
  static int tag;
  bool binary;
  bool absolute_dmin;
  double density;
  double dmin;
  double rmin;
  double pmax;
  double sign_normals;
  int seed;
  class RanPark *random;
  std::vector<double> radii;
  std::vector<double> sx;
  std::vector<double> sy;
  std::vector<double> sz;
};

}

#endif
#endif
#endif
