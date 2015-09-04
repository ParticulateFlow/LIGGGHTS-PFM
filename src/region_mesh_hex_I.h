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
    hexahedron->GetPointIds()->SetId(0,0);
    hexahedron->GetPointIds()->SetId(1,1);
    hexahedron->GetPointIds()->SetId(2,2);
    hexahedron->GetPointIds()->SetId(3,3);
    hexahedron->GetPointIds()->SetId(4,4);
    hexahedron->GetPointIds()->SetId(5,5);
    hexahedron->GetPointIds()->SetId(6,6);
    hexahedron->GetPointIds()->SetId(7,7);

    hexahedron->GetPoints()->SetPoint(0, node[iHex][0][0], node[iHex][0][1], node[iHex][0][2]);
    hexahedron->GetPoints()->SetPoint(1, node[iHex][1][0], node[iHex][1][1], node[iHex][1][2]);
    hexahedron->GetPoints()->SetPoint(2, node[iHex][2][0], node[iHex][2][1], node[iHex][2][2]);
    hexahedron->GetPoints()->SetPoint(3, node[iHex][3][0], node[iHex][3][1], node[iHex][3][2]);
    hexahedron->GetPoints()->SetPoint(4, node[iHex][4][0], node[iHex][4][1], node[iHex][4][2]);
    hexahedron->GetPoints()->SetPoint(5, node[iHex][5][0], node[iHex][5][1], node[iHex][5][2]);
    hexahedron->GetPoints()->SetPoint(6, node[iHex][6][0], node[iHex][6][1], node[iHex][6][2]);
    hexahedron->GetPoints()->SetPoint(7, node[iHex][7][0], node[iHex][7][1], node[iHex][7][2]);

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
