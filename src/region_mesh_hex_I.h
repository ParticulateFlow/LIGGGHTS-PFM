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


#ifndef LMP_REGION_MESH_HEX_I_H
#define LMP_REGION_MESH_HEX_I_H

/* ---------------------------------------------------------------------- */

inline void hex_randpos(int iHex, double *pos)
{
    double pcoords[3];

    pcoords[0] = random->uniform();
    pcoords[1] = random->uniform();
    pcoords[2] = random->uniform();

    pcoords_to_cart(iHex, pcoords, pos);
}

/* ---------------------------------------------------------------------- */

inline void pcoords_to_cart(int iHex, double *pcoords, double *pos)
{
    vtkHexahedron *hexahedron = vtkHexahedron::New();

    hexahedron->GetPointIds()->SetNumberOfIds(8);
    for(int i=0; i<8; ++i) {
        hexahedron->GetPointIds()->SetId(i,i);
        hexahedron->GetPoints()->SetPoint(i, node[iHex][i][0], node[iHex][i][1], node[iHex][i][2]);
    }

    int subId;
    double weights[8]={};
    // parametric coordinates to global coordinates
    /*vtkHexahedron::EvaluateLocation(int& vtkNotUsed(subId), double pcoords[3],
                                       double x[3], double *weights)*/
    hexahedron->EvaluateLocation(subId, pcoords, pos, weights);
    hexahedron->Delete();
}



/* ---------------------------------------------------------------------- */

#endif
