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
#include <vtkNew.h>
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
  // NOTE: no distinction between internal and external faces is made!

  // check subdomain
  if(!domain->is_in_subdomain(x)) return 0;

  // check bbox, only if exists
  if(bboxflag) {
    if(x[0] < extent_xlo || x[0] > extent_xhi) return 0;
    if(x[1] < extent_ylo || x[1] > extent_yhi) return 0;
    if(x[2] < extent_zlo || x[2] > extent_zhi) return 0;
  }

  // brute force naive search
  vtkNew<vtkHexahedron> hexahedron;
  hexahedron->GetPointIds()->SetNumberOfIds(8);

  for(int iHex=0; iHex<nHex; ++iHex) {
    for(int j=0; j<8; ++j) {
      hexahedron->GetPointIds()->SetId(j,j);
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double hexahedronCoords[3], hexahedronWeights[8];
    int subId;
    double dist2 = 0.;

    int result = hexahedron->EvaluatePosition(x, NULL, subId, hexahedronCoords, dist2, hexahedronWeights);
    if(result > 0) {
      // if inside hexahedron, dist2 is set to 0; need to check each face for distance
      double quadCoords[3], quadWeights[4], quadClosest[3];
      const int nFaces = hexahedron->GetNumberOfFaces();
      for(int i=0; i<nFaces; ++i) {
        vtkCell *quad = hexahedron->GetFace(i);
        result = quad->EvaluatePosition(x, quadClosest, subId, quadCoords, dist2, quadWeights);
        if(result > 0) {
          if(dist2 < cutoff*cutoff) {
            add_contact(0,x,quadClosest[0],quadClosest[1],quadClosest[2]);
            return 1;
          }
        }
      }
      return 0;
    }
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::surface_exterior(double *x, double cutoff)
{
  // NOTE: no distinction between internal and external faces is made!

  // check subdomain
  if(!domain->is_in_subdomain(x)) return 0;

  int nearby = 0;

  // brute force naive search
  vtkNew<vtkHexahedron> hexahedron;
  hexahedron->GetPointIds()->SetNumberOfIds(8);

  for(int iHex=0; iHex<nHex; ++iHex) {
    for(int j=0; j<8; ++j) {
      hexahedron->GetPointIds()->SetId(j,j);
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double hexahedronCoords[3], hexahedronWeights[8], hexahedronClosest[3];
    int subId;
    double dist2 = 0.;

    int result = hexahedron->EvaluatePosition(x, hexahedronClosest, subId, hexahedronCoords, dist2, hexahedronWeights);
    if(result == 0) {
      if(dist2 < cutoff*cutoff) {
        add_contact(0,x,hexahedronClosest[0],hexahedronClosest[1],hexahedronClosest[2]);
        nearby = 1;
      }
    } else {
      return 0; // inside
    }
  }

  return nearby;
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
  return volume_of_hex(node[iHex]);
}

/* ---------------------------------------------------------------------- */

inline int RegHexMesh::is_inside_hex(int iHex,double *pos)
{
  vtkNew<vtkHexahedron> hexahedron;
  hexahedron->GetPointIds()->SetNumberOfIds(8);

  for(int i=0; i<8; ++i) {
    hexahedron->GetPointIds()->SetId(i,i);
    hexahedron->GetPoints()->SetPoint(7, node[iHex][7][0], node[iHex][7][1], node[iHex][7][2]);
  }

  double hexahedronCoords[3], hexahedronWeights[8];//, hexahedronClosest[3];
  int subId;
  double dist2 = 0.;

  if(hexahedron->EvaluatePosition(pos, NULL, subId, hexahedronCoords, dist2, hexahedronWeights) > 0)
    return 1;
  return 0;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::volume_of_hex(double** v)
{
  vtkNew<vtkHexahedron> hexahedron;
  hexahedron->GetPointIds()->SetNumberOfIds(8);

  for(int i=0; i<8; ++i) {
    hexahedron->GetPointIds()->SetId(i,i);
    hexahedron->GetPoints()->SetPoint(0, v[i][0], v[i][1], v[i][2]);
  }

  return vtkMeshQuality::HexVolume(hexahedron.GetPointer());
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
