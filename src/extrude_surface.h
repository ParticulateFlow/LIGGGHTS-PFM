/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2016- JKU Linz

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

#if defined(LAMMPS_VTK) //NP do not use #ifdef here (VS C++ bug)
#ifdef COMMAND_CLASS

CommandStyle(extrude_surface,ExtrudeSurface)

#else

#ifndef LMP_EXTRUDE_SURFACE_H
#define LMP_EXTRUDE_SURFACE_H

#include "pointers.h"
#include <vtkDataSet.h>
#include <vtkDataArray.h>

namespace LAMMPS_NS {

  /**
   * @brief ExtrudeSurface class
   *        extrude surface from unstructured grid or polydata vtk file and write data to vtk file.
   *
   * Uses the vtk library to read an unstructured grid or polydata dataset from vtk simple legacy
   * or xml format, triangulate the polygons and write the data to a vtk file of type polydata.
   * IDs for all cell faces are added if not present so all triangles forming a face will have
   * the same face ID.
   */
class ExtrudeSurface : protected Pointers {
 public:
  ExtrudeSurface(class LAMMPS *);
  ~ExtrudeSurface();
  void command(int, char **);

 private:
  int me;
  bool binary;
  template<class TReader> vtkDataSet *read_file(const char*filename);
  void triangulate(int narg, char **arg, vtkDataSet* dset);
  void extrude(int narg, char **arg, vtkDataSet* dset);
  bool collinear(double *a, double *b, double *c);
  void extrude_point_via_normal(double x[3], vtkIdType id, vtkDataArray *n, double scale);
};

}

#endif
#endif
#endif
