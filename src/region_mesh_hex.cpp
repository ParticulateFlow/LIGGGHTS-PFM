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
  if(bboxflag) {
    if(pos[0] < extent_xlo || pos[0] > extent_xhi) return 0;
    if(pos[1] < extent_ylo || pos[1] > extent_yhi) return 0;
    if(pos[2] < extent_zlo || pos[2] > extent_zhi) return 0;
  }

  // brute force naive search
  for(int i=0; i<nHex; ++i) {
    if(is_inside_hex(i,pos) > 0) return 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

// test if point inside hex AND within a minimum distance from surface
int RegHexMesh::match_hex_cut(int iHex, double *pos, double cutoff)
{
  double x[3];
  vectorCopy3D(pos,x);

  if(interior) {
    vtkNew<vtkHexahedron> hexahedron;
    hexahedron->GetPointIds()->SetNumberOfIds(8);

    for(int j=0; j<8; ++j) {
      hexahedron->GetPointIds()->SetId(j,j);
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double *bounds = hexahedron->GetBounds();
    if(x[0] < bounds[0] || x[0] > bounds[1]) return 0; // outside bb of hex
    if(x[1] < bounds[2] || x[1] > bounds[3]) return 0; // outside bb of hex
    if(x[2] < bounds[4] || x[2] > bounds[5]) return 0; // outside bb of hex


    double hexahedronCoords[3], hexahedronWeights[8];
    int subId;
    double dist2 = 0.;

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, NULL, subId, hexahedronCoords, dist2, hexahedronWeights);

    if(result > 0) {
      // if inside hexahedron, dist2 is set to 0; need to check each face for distance
      double quadCoords[3], quadWeights[4], quadClosest[3];
      bool undefinedQuad = false;
      const int nFaces = hexahedron->GetNumberOfFaces();
      for(int i=0; i<nFaces; ++i) {
        vtkCell *quad = hexahedron->GetFace(i);
        result = quad->EvaluatePosition(x, quadClosest, subId, quadCoords, dist2, quadWeights);
        if(result > 0) {
          if(dist2 < cutoff*cutoff) {
            add_contact(0,x,quadClosest[0],quadClosest[1],quadClosest[2]);
            return 1;
          }
        } else if (result < 0) {
          undefinedQuad = true;
        }
      }
      if(undefinedQuad) {
        add_contact(0,x,x[0],x[1],x[2]);
        return 1;
      }

      return 0;
    } else if (result < 0) {
      add_contact(0,x,x[0],x[1],x[2]);
      return 1;
    }

    return 0;
  }
  else return surface_exterior(x,cutoff);
}

/* ---------------------------------------------------------------------- */
// test if point inside mesh AND within a minimum distance from surface
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

  bool undefinedHex = false;

  for(int iHex=0; iHex<nHex; ++iHex) {
    for(int j=0; j<8; ++j) {
      hexahedron->GetPointIds()->SetId(j,j);
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double hexahedronCoords[3], hexahedronWeights[8];
    int subId;
    double dist2 = 0.;

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, NULL, subId, hexahedronCoords, dist2, hexahedronWeights);
    if(result > 0) {
      // if inside hexahedron, dist2 is set to 0; need to check each face for distance
      double quadCoords[3], quadWeights[4], quadClosest[3];
      bool undefinedQuad = false;
      const int nFaces = hexahedron->GetNumberOfFaces();
      for(int i=0; i<nFaces; ++i) {
        vtkCell *quad = hexahedron->GetFace(i);
        result = quad->EvaluatePosition(x, quadClosest, subId, quadCoords, dist2, quadWeights);
        if(result > 0) {
          if(dist2 < cutoff*cutoff) {
            add_contact(0,x,quadClosest[0],quadClosest[1],quadClosest[2]);
            return 1;
          }
        } else if (result < 0) {
          undefinedQuad = true;
        }
      }
      if(undefinedQuad) {
        add_contact(0,x,x[0],x[1],x[2]);
        return 1;
      }

      return 0;
    } else if (result < 0) {
      undefinedHex = true;
    }
  }

  if(undefinedHex) {
    add_contact(0,x,x[0],x[1],x[2]);
    return 1;
  }

  return 0;
}

/* ---------------------------------------------------------------------- */

// test if point outside mesh AND within a minimum distance from surface
int RegHexMesh::surface_exterior(double *x, double cutoff)
{
  // NOTE: no distinction between internal and external faces is made!

  // check subdomain
  if(!domain->is_in_subdomain(x)) return 0;

  if(bboxflag) {
    if(x[0] < extent_xlo-cutoff || x[0] > extent_xhi+cutoff) return 0;
    if(x[1] < extent_ylo-cutoff || x[1] > extent_yhi+cutoff) return 0;
    if(x[2] < extent_zlo-cutoff || x[2] > extent_zhi+cutoff) return 0;
  }

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

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, hexahedronClosest, subId, hexahedronCoords, dist2, hexahedronWeights);
    if(result == 0) {
      if(dist2 < cutoff*cutoff) {
        add_contact(0,x,hexahedronClosest[0],hexahedronClosest[1],hexahedronClosest[2]);
        nearby = 1;
      }
    } else if(result < 0) {
      double *bounds = hexahedron->GetBounds();
      if(x[0] < bounds[0]-cutoff || x[0] > bounds[1]+cutoff) continue;
      if(x[1] < bounds[2]-cutoff || x[1] > bounds[3]+cutoff) continue;
      if(x[2] < bounds[4]-cutoff || x[2] > bounds[5]+cutoff) continue;

      add_contact(0,x,x[0],x[1],x[2]);
      nearby = 1;
    } else { // result == 1
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
  for(int i=0; i<8; ++i) {
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

void RegHexMesh::hex_bounds(int iHex, double bounds[6])
{
  vtkNew<vtkHexahedron> hexahedron;
  hexahedron->GetPointIds()->SetNumberOfIds(8);

  for(int i=0; i<8; ++i) {
    hexahedron->GetPointIds()->SetId(i,i);
    hexahedron->GetPoints()->SetPoint(i, node[iHex][i][0], node[iHex][i][1], node[iHex][i][2]);
  }

  hexahedron->GetBounds(bounds);
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::get_hex(const char* property, int value)
{
  ScalarContainer<int> *values = customValues_.getElementProperty<ScalarContainer<int> >(property);
  if(values) {
    for(int iHex=0; iHex<nHex; ++iHex) {
      if((*values)(iHex) == value)
        return iHex;
    }
  }
  return -1;
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
    hexahedron->GetPoints()->SetPoint(i, node[iHex][i][0], node[iHex][i][1], node[iHex][i][2]);
  }

  double *bounds = hexahedron->GetBounds();
  if(pos[0] < bounds[0] || pos[0] > bounds[1]) return 0; // outside bb of hex
  if(pos[1] < bounds[2] || pos[1] > bounds[3]) return 0; // outside bb of hex
  if(pos[2] < bounds[4] || pos[2] > bounds[5]) return 0; // outside bb of hex

  double hexahedronCoords[3], hexahedronWeights[8];
  int subId;
  double dist2 = 0.;

  // inside(=1), outside(=0) cell, or (-1) computational problem encountered
  int result = hexahedron->EvaluatePosition(pos, NULL, subId, hexahedronCoords, dist2, hexahedronWeights);
  if(result > 0)
    return 1;
  else if(result < 0) // if in doubt, assume position is inside, cause this function is used to determine
                      // particles that potentially overlap with newly inserted particles in this region
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
    hexahedron->GetPoints()->SetPoint(i, v[i][0], v[i][1], v[i][2]);
  }

  return vtkMeshQuality::HexVolume(hexahedron.GetPointer());
}

/* ---------------------------------------------------------------------- */

inline void RegHexMesh::set_extent()
{
  extent_xlo = extent_ylo = extent_zlo =  BIG;
  extent_xhi = extent_yhi = extent_zhi = -BIG;

  for(int i=0; i<nHex; ++i) {
    for(int j=0; j<8; ++j) {
      if(node[i][j][0] < extent_xlo) extent_xlo = node[i][j][0];
      if(node[i][j][1] < extent_ylo) extent_ylo = node[i][j][1];
      if(node[i][j][2] < extent_zlo) extent_zlo = node[i][j][2];

      if(node[i][j][0] > extent_xhi) extent_xhi = node[i][j][0];
      if(node[i][j][1] > extent_yhi) extent_yhi = node[i][j][1];
      if(node[i][j][2] > extent_zhi) extent_zhi = node[i][j][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

inline void RegHexMesh::mesh_randpos(double *pos)
{
  hex_randpos(mesh_randcell(), pos);
}

/* ---------------------------------------------------------------------- */

inline int RegHexMesh::mesh_randcell()
{
  double rd = total_volume * random->uniform();
  int chosen = 0;
  while (rd > acc_volume[chosen] && chosen < nHex-1) ++chosen;
  return chosen;
}

#endif
