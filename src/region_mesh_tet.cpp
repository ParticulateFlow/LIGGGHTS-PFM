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

#include <stdlib.h>
#include <string.h>
#include <list>
#include <stack>

#include "region_mesh_tet.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include <math.h>
#include "math_extra_liggghts.h"
#include "input_mesh_tet.h"

#define DELTA_TET 1000
#define BIG 1.e20

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;

/* ---------------------------------------------------------------------- */

double const RegTetMesh::phi = (1.+sqrt(5.))/2.;
int const RegTetMesh::n_ico_point = 12;

RegTetMesh::RegTetMesh(LAMMPS *lmp, int narg, char **arg) :
  Region(lmp, narg, arg)
{
  if(narg < 14) error->all(FLERR,"Illegal region mesh/tet command");
  options(narg-14,&arg[14]);

  if(scaleflag) error->all(FLERR,"Lattice scaling not implemented for region mesh/tet, please use 'units box'");

  if(strcmp(arg[2],"file"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'scale'");
  char *filename = arg[3];

  if(strcmp(arg[4],"scale"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'scale'");
  scale_fact = atof(arg[5]);
  if(strcmp(arg[6],"move"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'move'");
  off_fact[0] = atof(arg[7]);
  off_fact[1] = atof(arg[8]);
  off_fact[2] = atof(arg[9]);
  if(strcmp(arg[10],"rotate"))
    error->all(FLERR,"Illegal region mesh/tet command, expecting keyword 'rotate'");
  rot_angle[0] = atof(arg[11]);
  rot_angle[1] = atof(arg[12]);
  rot_angle[2] = atof(arg[13]);

  node = NULL;
  center = NULL;
  volume = NULL;
  acc_volume = NULL;
  nTet = 0;
  nTetMax = 0;
  total_volume = 0.;

  // manage input
  InputMeshTet *my_input = new InputMeshTet(lmp, 0, NULL);
  my_input->meshtetfile(filename,this,true);
  delete my_input;

  // extent of sphere

  if (interior) {
    bboxflag = 1;
    set_extent();
  } else bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];

  precalc_ico_points();

  build_tree();

}

void RegTetMesh::rebuild()
{
  if (compDouble(domain->sublo[0],domain_sublo[0]) &&
      compDouble(domain->sublo[1],domain_sublo[1]) &&
      compDouble(domain->sublo[2],domain_sublo[2]) &&
      compDouble(domain->subhi[0],domain_subhi[0]) &&
      compDouble(domain->subhi[1],domain_subhi[1]) &&
      compDouble(domain->subhi[2],domain_subhi[2]))
    return;

  build_tree();
}

void RegTetMesh::precalc_ico_points()
{
  // icosaedron point
  ico_points = memory->create<double>(ico_points,n_ico_point,3,"icosaeder points");

  double const coord2 = sqrt(3)/phi;
  double const coord1 = coord2/phi;

  ico_points[0][0] = 0.;
  ico_points[0][1] = coord1;
  ico_points[0][2] = coord2;

  ico_points[1][0] = 0.;
  ico_points[1][1] = -coord1;
  ico_points[1][2] = -coord2;

  ico_points[2][0] = 0.;
  ico_points[2][1] = -coord1;
  ico_points[2][2] = coord2;

  ico_points[3][0] = 0.;
  ico_points[3][1] = coord1;
  ico_points[3][2] = -coord2;

  ico_points[4][0] = coord2;
  ico_points[4][1] = 0.;
  ico_points[4][2] = coord1;

  ico_points[5][0] = -coord2;
  ico_points[5][1] = 0.;
  ico_points[5][2] = -coord1;

  ico_points[6][0] = -coord2;
  ico_points[6][1] = 0.;
  ico_points[6][2] = coord1;

  ico_points[7][0] = coord2;
  ico_points[7][1] = 0.;
  ico_points[7][2] = -coord1;

  ico_points[8][0] = coord1;
  ico_points[8][1] = coord2;
  ico_points[8][2] = 0.;

  ico_points[9][0] = -coord1;
  ico_points[9][1] = -coord2;
  ico_points[9][2] = 0.;

  ico_points[10][0] = -coord1;
  ico_points[10][1] = coord2;
  ico_points[10][2] = 0.;

  ico_points[11][0] = coord1;
  ico_points[11][1] = -coord2;
  ico_points[11][2] = 0.;
}

/* ---------------------------------------------------------------------- */

RegTetMesh::~RegTetMesh()
{
  delete [] contact;

  memory->destroy(node);
  memory->destroy(center);
  memory->sfree(volume);
  memory->sfree(acc_volume);
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegTetMesh::inside(double x, double y, double z)
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

  // return tree_is_inside_bin(pos,tree_data[0]);
  // return tree_is_inside_recursive(pos,0);
  return tree_is_inside(pos);
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::surface_interior(double *x, double cutoff)
{
  // check subdomain
  if(!domain->is_in_subdomain(x)) return 0;

  int n_contact = 0;
  double point[3];

  // instead of solving the full surface/particle problem, check if
  // the 12 points of a surrounding icosaedron are inside the
  // domain.
  for(int i=0;i<n_ico_point;i++){
    vectorAddMultiply3D(x,ico_points[i],cutoff,point);
    if(!inside(point[0],point[1],point[2])) n_contact++;
  }

  return n_contact;
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::surface_exterior(double *x, double cutoff)
{
  // check subdomain
  if(!domain->is_in_subdomain(x)) return 0;

  int n_contact = 0;
  double point[3];

  for(int i=0;i<n_ico_point;i++){
    vectorAddMultiply3D(x,ico_points[i],cutoff,point);
    if(inside(x[0],x[1],x[2])) n_contact++;
  }

  return n_contact;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::generate_random(double *pos)
{
    // function actually not called at the moment
    if(!interior) error->all(FLERR,"Impossible to generate random points on tet mesh region with side = out");
    mesh_randpos(pos);
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::generate_random_cut(double *pos,double cut)
{
    // function actually not called at the moment
    if(!interior) error->all(FLERR,"Impossible to generate random points on tet mesh region with side = out");
    error->all(FLERR,"This feature is not available for tet mesh regions");
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::add_tet(double **n)
{
    double ctr[3];

    if(nTet == nTetMax) grow_arrays();

    vectorZeroize3D(ctr);
    for(int i=0;i<4;i++)
    {
        vectorCopy3D(n[i],node[nTet][i]);
        vectorAdd3D(ctr,node[nTet][i],ctr);
    }
    vectorScalarDiv3D(ctr,4.);
    vectorCopy3D(ctr,center[nTet]);

    double vol = volume_of_tet(nTet);
    if(vol < 0.)
    {
        // flip nodes 0 and 3
        double node0[3];
        vectorCopy3D(node[nTet][0],node0);
        vectorCopy3D(node[nTet][3],node[nTet][0]);
        vectorCopy3D(node0,node[nTet][3]);
    }

    vol = volume_of_tet(nTet);
    if(vol < 0.) error->all(FLERR,"Fatal error: RegTetMesh::add_tet: vol < 0");

    volume[nTet] = vol;
    total_volume += volume[nTet];
    acc_volume[nTet] = volume[nTet];
    if(nTet > 0) acc_volume[nTet] += acc_volume[nTet-1];

    tet_bbox.push_back(BoundingBox());
    for(int i=0;i<4;i++)
      tet_bbox[nTet].extendToContain(n[i]);

    nTet++;
}

/* ---------------------------------------------------------------------- */

void RegTetMesh::grow_arrays()
{
    nTetMax += DELTA_TET;
    node = (double***)(memory->grow(node,nTetMax, 4 , 3, "vtk_tet_node"));
    center = (double**)(memory->grow(center,nTetMax, 3, "vtk_tet_center"));
    volume = (double*)(memory->srealloc(volume,nTetMax*sizeof(double),"vtk_tet_volume"));
    acc_volume = (double*)(memory->srealloc(acc_volume,nTetMax*sizeof(double),"vtk_tet_acc_volume"));
}

/* ---------------------------------------------------------------------- */

int RegTetMesh::n_tet()
{
    return nTet;
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::total_vol()
{
    return total_volume;
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::tet_vol(int i)
{
    return volume[i];
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::tet_acc_vol(int i)
{
    return acc_volume[i];
}

/* ---------------------------------------------------------------------- */

inline double RegTetMesh::volume_of_tet(int iTet)
{
    return volume_of_tet(node[iTet][0],node[iTet][1],node[iTet][2],node[iTet][3]);
}

/* ---------------------------------------------------------------------- */

inline int RegTetMesh::is_inside_tet(int iTet,double *pos)
{

  if(!tet_bbox[iTet].isInside(pos)) return 0;

  double vol1,vol2,vol3,vol4;

  vol1 = volume_of_tet(node[iTet][0], node[iTet][1], node[iTet][2], pos          );
  if(vol1<0) return 0;

  vol2 = volume_of_tet(node[iTet][0], node[iTet][1], pos,           node[iTet][3]);
  if(vol2<0) return 0;

  vol3 = volume_of_tet(node[iTet][0], pos,           node[iTet][2], node[iTet][3]);
  if(vol3<0) return 0;

  vol4 = volume_of_tet(pos          , node[iTet][1], node[iTet][2], node[iTet][3]);
  return static_cast<int>(vol4>0);
}

/* ---------------------------------------------------------------------- */

double RegTetMesh::volume_of_tet(double* v0, double* v1, double* v2, double* v3)
{
   double A[3];
   A[0] = v3[0] - v1[0];
   A[1] = v3[1] - v1[1];
   A[2] = v3[2] - v1[2];

   double B[3];
   B[0] = v2[0] - v1[0];
   B[1] = v2[1] - v1[1];
   B[2] = v2[2] - v1[2];

   double C[3];
   C[0] = v0[0] - v1[0];
   C[1] = v0[1] - v1[1];
   C[2] = v0[2] - v1[2];

   // cross product A x B
   double cp[] = {
       A[1]*B[2] - A[2]*B[1],
      -A[0]*B[2] + A[2]*B[0],
       A[0]*B[1] - A[1]*B[0]
   };

   // dot with C
   double volume = cp[0] * C[0] + cp[1] * C[1] + cp[2] * C[2];
   volume /= 6.;
   return volume;
}

/* ---------------------------------------------------------------------- */

inline void RegTetMesh::set_extent()
{
    extent_xlo = extent_ylo = extent_zlo =  BIG;
    extent_xhi = extent_yhi = extent_zhi = -BIG;

    for(int i = 0; i < nTet; i++)
        for(int j=0;j<4;j++)
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

inline void RegTetMesh::mesh_randpos(double *pos)
{
    tet_randpos(tet_rand_tri(),pos);
    if(pos[0] == 0. && pos[1] == 0. && pos[2] == 0.)
        error->one(FLERR,"illegal RegTetMesh::mesh_randpos");
}

/* ---------------------------------------------------------------------- */

inline int RegTetMesh::tet_rand_tri()
{

    double rd = total_volume * random->uniform();
    int chosen = 0;
    while (rd > acc_volume[chosen] && chosen < nTet-1) chosen++;
    return chosen;
}


void RegTetMesh::build_tree()
{

  tree_key.clear();
  tree_data.clear();

  tree_key.push_back(BoundingBox(extent_xlo,extent_xhi,
                                 extent_ylo,extent_yhi,
                                 extent_zlo,extent_zhi));

  vectorCopy3D(domain->sublo, domain_sublo);
  vectorCopy3D(domain->subhi, domain_subhi);

  tree_key[0].shrinkToSubbox(domain->sublo,domain->subhi);

  tree_data.push_back(TreeBin());

  tree_populate_node(0);

  // this is empirical - on lowest level, assume that
  // half of the nodes hold TREE_MIN_ELEMENTS_PER_NODE
  // and half of the nodes are empty. Could be experimented with
  // a little more if efficiency is an issue
  tree_max_depth = 1;
  while(nTet/pow(2,tree_max_depth) > TREE_MIN_ELEMENTS_PER_NODE/2)
    ++tree_max_depth;

  for(int iNode=0;iNode<pow(2,tree_max_depth)-1;++iNode){
    tree_create_children(iNode);
  }
}

void RegTetMesh::tree_populate_node(int iTreeNode)
{
  BoundingBox &bbox = tree_key[iTreeNode];
  TreeBin &data = tree_data[iTreeNode];

  if(iTreeNode == 0){
    for(int iTet=0;iTet<nTet;iTet++){
      for(int iNode=0;iNode<4;iNode++){
        if(bbox.isInside(node[iTet][iNode])){
          data.insert(iTet);
          break;
        }
      }
    }
    extend_bb(bbox,data);

  } else{
    TreeBin &parent_data = tree_data[tree_parent(iTreeNode)];
    if(parent_data.size() < TREE_MIN_ELEMENTS_PER_NODE)
      return; // nothing to do
    for(TreeBin::iterator it=parent_data.begin();it!=parent_data.end();++it){
      for(int iNode=0;iNode<4;iNode++){
        if(bbox.isInside(node[*it][iNode])){
          data.insert(*it);
          break;
        }
      }

    }
    extend_bb(bbox,data);

  }


  // if data size is equal, then splitting did not result in any gain
  // thus, empty the bin again
}

void RegTetMesh::tree_create_children(int current)
{
  int left_ind = tree_left(current);
  int right_ind = tree_right(current);

  bool const lower = true;
  tree_key.push_back(split_bbox_largest_extent(tree_key[current],lower));
  tree_key.push_back(split_bbox_largest_extent(tree_key[current],!lower));

  tree_data.push_back(TreeBin());
  tree_data.push_back(TreeBin());

  if(!tree_data[left_ind].empty())
    error->one(FLERR,"internal error in region mesh/tet, left child of search tree not empty");

  if(!tree_data[right_ind].empty())
    error->one(FLERR,"internal error in region mesh/tet, right child of search tree not empty");

  tree_populate_node(right_ind);
  tree_populate_node(left_ind);

  // if both new child nodes contain all tets from the parent element,
  // then the parent node should be considered a leaf
  if(tree_data[current].size() == tree_data[left_ind].size() &&
     tree_data[current].size() == tree_data[right_ind].size()){
    tree_data[left_ind].clear();
    tree_data[right_ind].clear();
  }

}

// lower = true --> return lower half
// else return upper half
BoundingBox RegTetMesh::split_bbox_largest_extent(BoundingBox &orig,bool lower)
{
  // build initial split
  double lo[3],hi[3],extent[3];

  orig.getBoxBounds(lo,hi);
  orig.getExtent(extent);

  double max_extent = 0;
  int extent_ind = 0;
  for(int i = 0;i<3;i++){
    if(extent[i] > max_extent){
      max_extent = extent[i];
      extent_ind = i;
    }
  }
  if(lower)
    hi[extent_ind] -= 0.5*max_extent;
  else
    lo[extent_ind] += 0.5*max_extent;

  return BoundingBox(lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
}

void RegTetMesh::extend_bb(BoundingBox &box, TreeBin const &data)
{
  double lo[3],hi[3];
  box.getBoxBounds(lo,hi);

  for(TreeBin::const_iterator it = data.begin(); it != data.end(); ++it){
    box.extendToContain(tet_bbox[*it]);
  }

  box.getBoxBounds(lo,hi);
}

bool RegTetMesh::tree_is_inside(double *x)
{

  std::list<int> nodeList;
  std::stack<int> s;

  s.push(0);
  while(!s.empty()){
    int current = s.top();
    s.pop();
    if(tree_key[current].isInside(x) && !tree_data[current].empty()){
      if(tree_is_leaf(current))
         nodeList.push_back(current);
      else{
        int index_left = tree_left(current);
        int index_right = tree_right(current);
        if(index_left < tree_size()) s.push(index_left);
        if(index_right < tree_size()) s.push(index_right);
      }
    }
  }

  nodeList.sort();

  for(std::list<int>::reverse_iterator it=nodeList.rbegin();it!=nodeList.rend();++it){
    if(tree_is_inside_bin(x,tree_data[*it])){
      return true;
    }
  }
  return false;
}

bool RegTetMesh::tree_is_inside_bin(double *x, TreeBin const &data)
{
  if(data.empty())
    return false;

  for(TreeBin::iterator it=data.begin();it!=data.end();++it)
    if(is_inside_tet(*it,x))
      return true;

  return false;
}

