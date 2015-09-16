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
#ifndef LMP_INPUT_MESH_HEX_H
#define LMP_INPUT_MESH_HEX_H

#include "stdio.h"
#include "input.h"

namespace LAMMPS_NS {

class InputMeshHex : protected Input {
 public:

  InputMeshHex(class LAMMPS *, int, char **);
  ~InputMeshHex();

  void meshhexfile(class RegHexMesh *);
  void meshhexfile(const char *, class RegHexMesh *, bool verbose, bool read_cell_data=false); // analogon to file(const char *filename)

 private:
  void meshhexfile_vtk(class RegHexMesh *);
  bool verbose_;
  bool read_cell_data_;
};

}

#endif
#endif
