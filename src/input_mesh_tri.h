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

#ifndef LMP_INPUT_MESH_TRI_H
#define LMP_INPUT_MESH_TRI_H

#include <stdio.h>
#include "input.h"

namespace LAMMPS_NS {

class InputMeshTri : protected Input
{
  public:

    InputMeshTri(class LAMMPS *, int, char **);
    ~InputMeshTri();

    void meshtrifile(const char *,class TriMesh *,bool verbose,
                     const int size_exclusion_list,int *exclusion_list,
                     bool read_cell_data=false, bool restart=false);

  private:

    bool verbose_;
    int i_exclusion_list_;
    int size_exclusion_list_;
    int *exclusion_list_;
    bool read_cell_data_;
    bool restart_;

    void meshtrifile_vtk(class TriMesh *);
    void meshtrifile_stl(class TriMesh *);
    inline void addTriangle(class TriMesh *mesh,
         double *a, double *b, double *c,int lineNumber);

};

}

#endif
