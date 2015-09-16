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

#ifdef LAMMPS_VTK
#include "stdlib.h"
#include "string.h"
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkMeshQuality.h>
#include "region_mesh_hex.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "input_mesh_hex.h"
#define DELTA_HEX 1000
#define BIG 1.e20

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegHexMesh::RegHexMesh(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg),
  read_cell_data_(false),
  customValues_(*(new CustomValueTracker(lmp)))
{
  if(narg < 14) error->all(FLERR,"Illegal region mesh/hex command");

  if(strcmp(arg[2],"file"))
    error->all(FLERR,"Illegal region mesh/hex command, expecting keyword 'file'");
  char *filename = arg[3];

  if(strcmp(arg[4],"scale"))
    error->all(FLERR,"Illegal region mesh/hex command, expecting keyword 'scale'");
  scale_fact = atof(arg[5]);
  if(strcmp(arg[6],"move"))
    error->all(FLERR,"Illegal region mesh/hex command, expecting keyword 'move'");
  off_fact[0] = atof(arg[7]);
  off_fact[1] = atof(arg[8]);
  off_fact[2] = atof(arg[9]);
  if(strcmp(arg[10],"rotate"))
    error->all(FLERR,"Illegal region mesh/hex command, expecting keyword 'rotate'");
  rot_angle[0] = atof(arg[11]);
  rot_angle[1] = atof(arg[12]);
  rot_angle[2] = atof(arg[13]);

  int iarg = 14;
  if(narg > 15 && strcmp(arg[14],"cell_data") == 0) {
    if(strcmp(arg[15],"yes") == 0)
      read_cell_data_ = true;
    else if(strcmp(arg[15],"no"))
      error->all(FLERR,"Illegal region mesh/hex command, expecting 'yes' or 'no' for 'cell_data'");
    iarg = 16;
  }

  options(narg-iarg,&arg[iarg]);

  if(scaleflag) error->all(FLERR,"Lattice scaling not implemented for region mesh/hex, please use 'units box'");

  node = NULL;
  center = NULL;
  volume = NULL;
  acc_volume = NULL;
  nHex = 0;
  nHexMax = 0;
  total_volume = 0.;

  // manage input
  InputMeshHex *my_input = new InputMeshHex(lmp, 0, NULL);
  my_input->meshhexfile(filename,this,true,read_cell_data_);
  delete my_input;

  // extent of hexmesh

  if (interior) {
    bboxflag = 1;
    set_extent();
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
}

/* ---------------------------------------------------------------------- */

RegHexMesh::~RegHexMesh()
{
  delete [] contact;
  delete &customValues_;

  memory->destroy(node);
  memory->destroy(center);
  memory->sfree(volume);
  memory->sfree(acc_volume);
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegHexMesh::inside(double x, double y, double z)
{
   double pos[3];
   pos[0] = x; pos[1] = y; pos[2] = z;

   // check subdomain
   if(!domain->is_in_subdomain(pos)) return 0;

   // check bbox, only if exists
   if(bboxflag)
   {
       if(pos[0] < extent_xlo || pos[0] > extent_xhi) return 0;
       if(pos[1] < extent_ylo || pos[1] > extent_yhi) return 0;
       if(pos[2] < extent_zlo || pos[2] > extent_zhi) return 0;
   }

   // brute force naive search
   for(int i = 0; i < nHex; i++)
   {
      if(is_inside_hex(i,pos) > 0) return 1;
   }

   //fprintf(screen,"checking pos %f %f %f, result %d; ntet %d\n",x,y,z,inside_mesh,nTet);

   return 0;
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::surface_interior(double *x, double cutoff)
{
  error->one(FLERR,"This feature is not available for hex mesh regions");
  return 0;
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::surface_exterior(double *x, double cutoff)
{
  error->one(FLERR,"This feature is not available for hex mesh regions");
  return 0;
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::generate_random(double *pos)
{
    if(!interior) error->all(FLERR,"Impossible to generate random points on hex mesh region with side = out");
    mesh_randpos(pos);
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::generate_random_cut(double *pos,double cut)
{
    // function actually not called at the moment
    if(!interior) error->all(FLERR,"Impossible to generate random points on hex mesh region with side = out");
    error->all(FLERR,"This feature is not available for hex mesh regions");
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::add_hex(double **n)
{
    double ctr[3];

    if(nHex == nHexMax) grow_arrays();

    vectorZeroize3D(ctr);
    for(int i=0; i<8; ++i)
    {
        vectorCopy3D(n[i],node[nHex][i]);
        vectorAdd3D(ctr,node[nHex][i],ctr);
    }
    vectorScalarDiv3D(ctr,8.);
    vectorCopy3D(ctr,center[nHex]);

    volume[nHex] = volume_of_hex(nHex);
    total_volume += volume[nHex];
    acc_volume[nHex] = volume[nHex];
    if(nHex > 0) acc_volume[nHex] += acc_volume[nHex-1];
    ++nHex;
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::grow_arrays()
{
    nHexMax += DELTA_HEX;
    node = (double***)(memory->grow(node,nHexMax, 8, 3, "vtk_hex_node"));
    center = (double**)(memory->grow(center,nHexMax, 3, "vtk_hex_center"));
    volume = (double*)(memory->srealloc(volume,nHexMax*sizeof(double),"vtk_hex_volume"));
    acc_volume = (double*)(memory->srealloc(acc_volume,nHexMax*sizeof(double),"vtk_hex_acc_volume"));
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::n_hex()
{
    return nHex;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::total_vol()
{
    return total_volume;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::hex_vol(int i)
{
    return volume[i];
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::hex_acc_vol(int i)
{
    return acc_volume[i];
}

/* ---------------------------------------------------------------------- */

inline double RegHexMesh::volume_of_hex(int iHex)
{
    return volume_of_hex(node[iHex][0],node[iHex][1],node[iHex][2],node[iHex][3],
                         node[iHex][4],node[iHex][5],node[iHex][6],node[iHex][7]);
}

/* ---------------------------------------------------------------------- */

inline int RegHexMesh::is_inside_hex(int iHex,double *pos)
{
    /*double vol1,vol2,vol3,vol4;

    vol1 = volume_of_tet(node[iHex][0], node[iHex][1], node[iHex][2], pos          );
    vol2 = volume_of_tet(node[iHex][0], node[iHex][1], pos,           node[iHex][3]);
    vol3 = volume_of_tet(node[iHex][0], pos,           node[iHex][2], node[iHex][3]);
    vol4 = volume_of_tet(pos          , node[iHex][1], node[iHex][2], node[iHex][3]);

    if(vol1 > 0. && vol2 > 0. && vol3 > 0. && vol4 > 0.) return 1;*/
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

    double hexahedronCoords[3], hexahedronWeights[8];//, hexahedronClosest[3];
    int subId;
    double dist2;

    int result = hexahedron->EvaluatePosition(pos, NULL, subId, hexahedronCoords, dist2, hexahedronWeights);
    hexahedron->Delete();
    if(result > 0) return 1;
    return 0;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::volume_of_hex(double* v0, double* v1, double* v2, double* v3,
                                 double* v4, double* v5, double* v6, double* v7)
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

    hexahedron->GetPoints()->SetPoint(0, v0[0], v0[1], v0[2]);
    hexahedron->GetPoints()->SetPoint(1, v1[0], v1[1], v1[2]);
    hexahedron->GetPoints()->SetPoint(2, v2[0], v2[1], v2[2]);
    hexahedron->GetPoints()->SetPoint(3, v3[0], v3[1], v3[2]);
    hexahedron->GetPoints()->SetPoint(4, v4[0], v4[1], v4[2]);
    hexahedron->GetPoints()->SetPoint(5, v5[0], v5[1], v5[2]);
    hexahedron->GetPoints()->SetPoint(6, v6[0], v6[1], v6[2]);
    hexahedron->GetPoints()->SetPoint(7, v7[0], v7[1], v7[2]);

    double volume = vtkMeshQuality::HexVolume(hexahedron);
    hexahedron->Delete();
    return volume;
}

/* ---------------------------------------------------------------------- */

inline void RegHexMesh::set_extent()
{
    extent_xlo = extent_ylo = extent_zlo =  BIG;
    extent_xhi = extent_yhi = extent_zhi = -BIG;

    for(int i = 0; i < nHex; i++)
        for(int j=0;j<8;j++)
        {
            if(node[i][j][0] < extent_xlo) extent_xlo = node[i][j][0];
            if(node[i][j][1] < extent_ylo) extent_ylo = node[i][j][1];
            if(node[i][j][2] < extent_zlo) extent_zlo = node[i][j][2];

            if(node[i][j][0] > extent_xhi) extent_xhi = node[i][j][0];
            if(node[i][j][1] > extent_yhi) extent_yhi = node[i][j][1];
            if(node[i][j][2] > extent_zhi) extent_zhi = node[i][j][2];
        }
}

/* ---------------------------------------------------------------------- */

inline void RegHexMesh::mesh_randpos(double *pos)
{
    hex_randpos(mesh_randcell(), pos);
    //if(pos[0] == 0. && pos[1] == 0. && pos[2] == 0.)
    //    error->one(FLERR,"illegal RegHexMesh::mesh_randpos");
}

/* ---------------------------------------------------------------------- */

inline int RegHexMesh::mesh_randcell()
{

    double rd = total_volume * random->uniform();
    int chosen = 0;
    while (rd > acc_volume[chosen] && chosen < nHex-1) chosen++;
    return chosen;
}

#endif
