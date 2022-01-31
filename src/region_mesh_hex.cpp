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
#include <stdlib.h>
#include <string.h>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkMeshQuality.h>
#include <vtkTetra.h>
#include <float.h>
#include "region_mesh_hex.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "math.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "input_mesh_hex.h"
#define DELTA_HEX 1000
#define BIG 1.e20

namespace LAMMPS_NS {

class AABBNode {
 public:
  AABBNode() :
    left_(NULL), right_(NULL),
    xLo(0.), xHi(0.),
    yLo(0.), yHi(0.),
    zLo(0.), zHi(0.),
    cell_(-1) {}

  AABBNode(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, int cell=-1) :
    left_(NULL), right_(NULL),
    xLo(xlo), xHi(xhi),
    yLo(ylo), yHi(yhi),
    zLo(zlo), zHi(zhi),
    cell_(cell) {}

  AABBNode(double bounds[6], int cell=-1) :
    left_(NULL), right_(NULL),
    xLo(bounds[0]), xHi(bounds[1]),
    yLo(bounds[2]), yHi(bounds[3]),
    zLo(bounds[4]), zHi(bounds[5]),
    cell_(cell) {}

  ~AABBNode() { delete left_; delete right_; }


  void build(RegHexMesh *mesh, const std::vector<int>& subnodes)
  {
    int splitAxis = 0;
    double splitPlane[3] = { 0.5*(xHi + xLo), 0.5*(yHi + yLo), 0.5*(zHi + zLo) };
    if(xHi - xLo < yHi - yLo) {
      splitAxis = 1;
    }
    if(yHi - yLo < zHi - zLo) {
      splitAxis = 2;
    }

    std::vector<int> leftsubnodes;
    std::vector<int> rightsubnodes;
    double center[3];
    double bounds[6];

    do {
      leftsubnodes.clear();
      rightsubnodes.clear();
      for(unsigned int i=0; i<subnodes.size(); ++i) {
        mesh->hex_bounds(subnodes[i], bounds);
        getCenter(bounds, center);

        if(center[splitAxis] < splitPlane[splitAxis]) {
          leftsubnodes.push_back(subnodes[i]);
        } else {
          rightsubnodes.push_back(subnodes[i]);
        }
      }
      splitAxis = (splitAxis + 1) % 3; // may have to try a different axis
    } while(leftsubnodes.empty() || rightsubnodes.empty());

    if(!leftsubnodes.empty()) {
      mesh->hex_bounds(leftsubnodes[0], bounds);
      left_ = new AABBNode(bounds);
      if(leftsubnodes.size() < 2) {
        left_->cell_ = leftsubnodes[0];
      } else {
        for(unsigned int i=1; i<leftsubnodes.size(); ++i) {
          mesh->hex_bounds(leftsubnodes[i], bounds);
          left_->grow(bounds);
        }
        left_->build(mesh, leftsubnodes);
      }
    }

    if(!rightsubnodes.empty()) {
      mesh->hex_bounds(rightsubnodes[0], bounds);
      right_ = new AABBNode(bounds);
      if(rightsubnodes.size() < 2) {
        right_->cell_ = rightsubnodes[0];
      } else {
        for(unsigned int i=1; i<rightsubnodes.size(); ++i) {
          mesh->hex_bounds(rightsubnodes[i], bounds);
          right_->grow(bounds);
        }
        right_->build(mesh, rightsubnodes);
      }
    }
  }

  void getCenter(double bounds[6], double center[3])
  {
    center[0] = 0.5*(bounds[0] + bounds[1]);
    center[1] = 0.5*(bounds[2] + bounds[3]);
    center[2] = 0.5*(bounds[4] + bounds[5]);
  }

  void findCell(const double *p, std::vector<int>& potential_cells)
  {
    if(isInside(p)) {
      if(isLeaf()) {
        potential_cells.push_back(cell_);
      } else {
        if(left_)  left_->findCell(p, potential_cells);
        if(right_) right_->findCell(p, potential_cells);
      }
    }
  }

  bool isInside(const double *p)
  {
    // test for >= and < as in Domain class
    return (p[0] >= xLo && p[0] < xHi &&
            p[1] >= yLo && p[1] < yHi &&
            p[2] >= zLo && p[2] < zHi);
  }

  bool isLeaf()
  {
    return (cell_ > -1);
  }

  void grow(double bounds[6])
  {
    if(bounds[0] < xLo) xLo = bounds[0];
    if(bounds[1] > xHi) xHi = bounds[1];
    if(bounds[2] < yLo) yLo = bounds[2];
    if(bounds[3] > yHi) yHi = bounds[3];
    if(bounds[4] < zLo) zLo = bounds[4];
    if(bounds[5] > zHi) zHi = bounds[5];
  }

 private:
  AABBNode *left_;
  AABBNode *right_;
  double xLo, xHi, yLo, yHi, zLo, zHi;
  int cell_;
};


class AABBTree {
 public:
  AABBTree() : root_(NULL) {}
  ~AABBTree() { delete root_; }

  void build(RegHexMesh *mesh)
  {
    int nodes = mesh->n_hex();
    root_ = new AABBNode(mesh->extent_xlo, mesh->extent_xhi,
                         mesh->extent_ylo, mesh->extent_yhi,
                         mesh->extent_zlo, mesh->extent_zhi,
                         (nodes > 1) ? -1 : 0);

    if(nodes > 1) {
      std::vector<int> subnodes(nodes, 0);
      for(int i=1; i<nodes; ++i) subnodes[i] = i;
      root_->build(mesh, subnodes);
    }
  }

  void findCell(const double *p, std::vector<int>& potential_cells)
  {
    if(root_)
      root_->findCell(p, potential_cells);
  }

 private:
  AABBNode *root_;
};

}

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegHexMesh::RegHexMesh(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg),
  read_cell_data_(false),
  customValues_(*(new CustomValueTracker(lmp))),
  tree_(NULL)
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
  orthogonal = NULL;
  nHex = 0;
  nHexMax = 0;
  total_volume = 0.;

  // manage input
  InputMeshHex *my_input = new InputMeshHex(lmp, 0, NULL);
  my_input->meshhexfile(filename,this,true,read_cell_data_);
  delete my_input;

  // extent of hexmesh

  set_extent();
  tree_ = new AABBTree();
  tree_->build(this);

  if(interior) {
    bboxflag = 1;
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];

  for(int i=0; i<8; ++i) {
    hexahedron->GetPointIds()->SetId(i,i);
  }
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
  memory->sfree(orthogonal);
  delete tree_;
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

  potential_cells.clear();
  tree_->findCell(pos, potential_cells);
  for(unsigned int i=0; i < potential_cells.size(); ++i) {
    if(is_inside_hex(potential_cells[i],pos) > 0) {
      return 1;
    }
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
    for(int j=0; j<8; ++j) {
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double *bounds = hexahedron->GetBounds();
    if(x[0] < bounds[0] || x[0] > bounds[1]) return 0; // outside bb of hex
    if(x[1] < bounds[2] || x[1] > bounds[3]) return 0; // outside bb of hex
    if(x[2] < bounds[4] || x[2] > bounds[5]) return 0; // outside bb of hex

    int subId;
    double dist2 = 0.;

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, NULL, subId, pCoords, dist2, weights);

    if(result > 0) {
      // if inside hexahedron, dist2 is set to 0; need to check each face for distance
      double quadClosest[3];
      bool undefinedQuad = false;
      const int nFaces = hexahedron->GetNumberOfFaces();
      for(int i=0; i<nFaces; ++i) {
        vtkCell *quad = hexahedron->GetFace(i);
        result = quad->EvaluatePosition(x, quadClosest, subId, pCoords, dist2, weights);
        if(result > 0) {
          if(dist2 < cutoff*cutoff) {
            add_contact(0,x,quadClosest[0],quadClosest[1],quadClosest[2]);
            return 1;
          }
        } else if(result < 0) {
          undefinedQuad = true;
        }
      }
      if(undefinedQuad) {
        add_contact(0,x,x[0],x[1],x[2]);
        return 1;
      }

      return 0;
    } else if(result < 0) {
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

  bool undefinedHex = false;

  potential_cells.clear();
  tree_->findCell(x, potential_cells);
  for(unsigned int i=0; i < potential_cells.size(); ++i) {
    for(int j=0; j<8; ++j) {
      hexahedron->GetPoints()->SetPoint(j, node[potential_cells[i]][j][0], node[potential_cells[i]][j][1], node[potential_cells[i]][j][2]);
    }

    int subId;
    double dist2 = 0.;

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, NULL, subId, pCoords, dist2, weights);
    if(result > 0) {
      // if inside hexahedron, dist2 is set to 0; need to check each face for distance
      double quadClosest[3];
      bool undefinedQuad = false;
      const int nFaces = hexahedron->GetNumberOfFaces();
      for(int i=0; i<nFaces; ++i) {
        vtkCell *quad = hexahedron->GetFace(i);
        result = quad->EvaluatePosition(x, quadClosest, subId, pCoords, dist2, weights);
        if(result > 0) {
          if(dist2 < cutoff*cutoff) {
            add_contact(0,x,quadClosest[0],quadClosest[1],quadClosest[2]);
            return 1;
          }
        } else if(result < 0) {
          undefinedQuad = true;
        }
      }

      if(undefinedQuad) {
        add_contact(0,x,x[0],x[1],x[2]);
        return 1;
      }

      return 0;
    } else if(result < 0) {
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
  for(int iHex=0; iHex<nHex; ++iHex) {
    for(int j=0; j<8; ++j) {
      hexahedron->GetPoints()->SetPoint(j, node[iHex][j][0], node[iHex][j][1], node[iHex][j][2]);
    }

    double hexahedronClosest[3];
    int subId;
    double dist2 = 0.;

    // inside(=1), outside(=0) cell, or (-1) computational problem encountered
    int result = hexahedron->EvaluatePosition(x, hexahedronClosest, subId, pCoords, dist2, weights);
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
  orthogonal[nHex] = is_cell_orthogonal(nHex);
  ++nHex;
}

/* ---------------------------------------------------------------------- */

bool RegHexMesh::is_cell_orthogonal(int iHex) const
{
  /*
    7        6
    *--------*
   /|       /|
 4*-+------*5|
  | |      | |
  |3*------+-*2
  |/       |/
  *--------*
  0        1
   */
  double edge[12][3];
  edge[0][0] = node[iHex][1][0] - node[iHex][0][0];
  edge[0][1] = node[iHex][1][1] - node[iHex][0][1];
  edge[0][2] = node[iHex][1][2] - node[iHex][0][2];

  edge[1][0] = node[iHex][2][0] - node[iHex][1][0];
  edge[1][1] = node[iHex][2][1] - node[iHex][1][1];
  edge[1][2] = node[iHex][2][2] - node[iHex][1][2];

  edge[2][0] = node[iHex][2][0] - node[iHex][3][0];
  edge[2][1] = node[iHex][2][1] - node[iHex][3][1];
  edge[2][2] = node[iHex][2][2] - node[iHex][3][2];

  edge[3][0] = node[iHex][3][0] - node[iHex][0][0];
  edge[3][1] = node[iHex][3][1] - node[iHex][0][1];
  edge[3][2] = node[iHex][3][2] - node[iHex][0][2];


  edge[4][0] = node[iHex][5][0] - node[iHex][4][0];
  edge[4][1] = node[iHex][5][1] - node[iHex][4][1];
  edge[4][2] = node[iHex][5][2] - node[iHex][4][2];

  edge[5][0] = node[iHex][6][0] - node[iHex][5][0];
  edge[5][1] = node[iHex][6][1] - node[iHex][5][1];
  edge[5][2] = node[iHex][6][2] - node[iHex][5][2];

  edge[6][0] = node[iHex][6][0] - node[iHex][7][0];
  edge[6][1] = node[iHex][6][1] - node[iHex][7][1];
  edge[6][2] = node[iHex][6][2] - node[iHex][7][2];

  edge[7][0] = node[iHex][7][0] - node[iHex][4][0];
  edge[7][1] = node[iHex][7][1] - node[iHex][4][1];
  edge[7][2] = node[iHex][7][2] - node[iHex][4][2];


  edge[8][0] = node[iHex][4][0] - node[iHex][0][0];
  edge[8][1] = node[iHex][4][1] - node[iHex][0][1];
  edge[8][2] = node[iHex][4][2] - node[iHex][0][2];

  edge[9][0] = node[iHex][5][0] - node[iHex][1][0];
  edge[9][1] = node[iHex][5][1] - node[iHex][1][1];
  edge[9][2] = node[iHex][5][2] - node[iHex][1][2];

  edge[10][0] = node[iHex][7][0] - node[iHex][3][0];
  edge[10][1] = node[iHex][7][1] - node[iHex][3][1];
  edge[10][2] = node[iHex][7][2] - node[iHex][3][2];

  edge[11][0] = node[iHex][6][0] - node[iHex][2][0];
  edge[11][1] = node[iHex][6][1] - node[iHex][2][1];
  edge[11][2] = node[iHex][6][2] - node[iHex][2][2];

#define SMALL_COS 0.001

  // corner 0
  MathExtra::norm3(edge[0]);
  MathExtra::norm3(edge[3]);
  if (fabs(MathExtra::dot3(edge[0],edge[3])) > SMALL_COS)
    return false;
  MathExtra::norm3(edge[8]);
  if (fabs(MathExtra::dot3(edge[0],edge[8])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[3],edge[8])) > SMALL_COS)
    return false;

  // corner 2
  MathExtra::norm3(edge[1]);
  MathExtra::norm3(edge[2]);
  if (fabs(MathExtra::dot3(edge[1],edge[2] )) > SMALL_COS)
    return false;
  MathExtra::norm3(edge[11]);
  if (fabs(MathExtra::dot3(edge[1],edge[11])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[2],edge[11])) > SMALL_COS)
    return false;

  // corner 3
  if (fabs(MathExtra::dot3(edge[2],edge[3] )) > SMALL_COS)
    return false;
  MathExtra::norm3(edge[10]);
  if (fabs(MathExtra::dot3(edge[2],edge[10])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[3],edge[10])) > SMALL_COS)
    return false;


  // corner 4
  MathExtra::norm3(edge[4]);
  MathExtra::norm3(edge[7]);
  if (fabs(MathExtra::dot3(edge[4],edge[7])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[4],edge[8])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[7],edge[8])) > SMALL_COS)
    return false;

  // corner 5
  MathExtra::norm3(edge[5]);
  if (fabs(MathExtra::dot3(edge[4],edge[5])) > SMALL_COS)
    return false;
  MathExtra::norm3(edge[9]);
  if (fabs(MathExtra::dot3(edge[4],edge[9])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[5],edge[9])) > SMALL_COS)
    return false;

  // corner 6
  MathExtra::norm3(edge[6]);
  if (fabs(MathExtra::dot3(edge[5],edge[6]) ) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[5],edge[11])) > SMALL_COS)
    return false;
  if (fabs(MathExtra::dot3(edge[6],edge[11])) > SMALL_COS)
    return false;

  // corner 1 and 7 must be orthogonal
  return true;
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::grow_arrays()
{
  nHexMax += DELTA_HEX;
  node = (double***)(memory->grow(node,nHexMax, 8, 3, "vtk_hex_node"));
  center = (double**)(memory->grow(center,nHexMax, 3, "vtk_hex_center"));
  volume = (double*)(memory->srealloc(volume,nHexMax*sizeof(double),"vtk_hex_volume"));
  acc_volume = (double*)(memory->srealloc(acc_volume,nHexMax*sizeof(double),"vtk_hex_acc_volume"));
  orthogonal = (bool*)(memory->srealloc(orthogonal,nHexMax*sizeof(bool),"vtk_hex_ortho"));
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::n_hex() const
{
  return nHex;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::total_vol() const
{
  return total_volume;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::hex_vol(int i) const
{
  return volume[i];
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::hex_acc_vol(int i) const
{
  return acc_volume[i];
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::hex_bounds(int iHex, double bounds[6]) const
{
  bounds[0] = node[iHex][0][0];
  bounds[2] = node[iHex][0][1];
  bounds[4] = node[iHex][0][2];

  bounds[1] = node[iHex][0][0];
  bounds[3] = node[iHex][0][1];
  bounds[5] = node[iHex][0][2];

  for(int j=1; j<8; ++j) {
    if(node[iHex][j][0] < bounds[0]) bounds[0] = node[iHex][j][0];
    if(node[iHex][j][1] < bounds[2]) bounds[2] = node[iHex][j][1];
    if(node[iHex][j][2] < bounds[4]) bounds[4] = node[iHex][j][2];

    if(node[iHex][j][0] > bounds[1]) bounds[1] = node[iHex][j][0];
    if(node[iHex][j][1] > bounds[3]) bounds[3] = node[iHex][j][1];
    if(node[iHex][j][2] > bounds[5]) bounds[5] = node[iHex][j][2];
  }
}

/* ---------------------------------------------------------------------- */

void RegHexMesh::hex_points(int iHex, double points[24]) const
{
  for(int i=0; i<8; ++i) {
    for(int j=0; j<3; ++j) {
      points[i*3+j] = node[iHex][i][j];
    }
  }
}

/* ---------------------------------------------------------------------- */

int RegHexMesh::get_hex(double *pos)
{
  potential_cells.clear();
  tree_->findCell(pos, potential_cells);
  for(unsigned int i=0; i < potential_cells.size(); ++i) {
    if(is_inside_hex(potential_cells[i],pos) > 0) {
      return potential_cells[i];
    }
  }

  return -1;
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

int RegHexMesh::is_inside_hex(int iHex,double *pos)
{
  double bounds[6];
  hex_bounds(iHex, bounds);
  if(pos[0] < bounds[0] || pos[0] > bounds[1]) return 0; // outside bb of hex
  if(pos[1] < bounds[2] || pos[1] > bounds[3]) return 0; // outside bb of hex
  if(pos[2] < bounds[4] || pos[2] > bounds[5]) return 0; // outside bb of hex

  if(orthogonal[iHex])
    return 1;

  for(int i=0; i<8; ++i) {
    hexahedron->GetPoints()->SetPoint(i, node[iHex][i][0], node[iHex][i][1], node[iHex][i][2]);
  }

  int subId;
  double dist2 = 0.;

  // inside(=1), outside(=0) cell, or (-1) computational problem encountered
  int result = hexahedron->EvaluatePosition(pos, NULL, subId, pCoords, dist2, weights);
  if(result > 0)
    return 1;
  else if(result < 0) {
    hexahedron->Triangulate(0, ptIds.GetPointer(), pts.GetPointer()); // generates 5 tetrahedra

    for(int t=0; t<5; ++t) {
      for(int i=0; i<4; ++i) {
        tetra->GetPointIds()->SetId(i,ptIds->GetId(4*t+i));
        tetra->GetPoints()->SetPoint(i,pts->GetPoint(4*t+i));
      }
      if(vtkMeshQuality::TetVolume(tetra.GetPointer()) < DBL_MIN)
        continue;
      result = tetra->EvaluatePosition(pos, NULL, subId, pCoords, dist2, weights);
      if(result > 0)
        return 1;
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double RegHexMesh::volume_of_hex(double** v)
{
  for(int i=0; i<8; ++i) {
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
