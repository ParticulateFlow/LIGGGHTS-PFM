/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2015- JKU Linz

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

CommandStyle(extract_surface,ExtractSurface)

#else

#ifndef LMP_EXTRACT_SURFACE_H
#define LMP_EXTRACT_SURFACE_H

//#include "stdio.h"
#include "pointers.h"
#include <vtkDataSet.h>
#include <vtkDataArray.h>

namespace LAMMPS_NS {

  /**
   * @brief ExtractSurface class
   *        extract surface from unstructured grid vtk file and write data to vtk file.
   *
   * Uses the vtk library to read an unstructured grid dataset composed of hexahedra from
   * vtk simple legacy or xml format, extract the surface, triangulate the polygons and
   * write the data to a vtk file of type polydata.
   * IDs for all hexahedral cells and all cell faces are added so all triangles forming a
   * face will have the same face ID and all triangles belonging to the same cell will
   * have the same cell ID.
   */
class ExtractSurface : protected Pointers {
 public:
  ExtractSurface(class LAMMPS *);
  ~ExtractSurface();
  void command(int, char **);

 private:
  int me;
  template<class TReader> vtkDataSet *read_file(const char*filename);
  void extrude_point_via_normal(double x[3], vtkIdType id, vtkDataArray *n, double scale);
};

}

#endif
#endif
#endif
