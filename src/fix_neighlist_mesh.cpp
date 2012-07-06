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

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "atom.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

/*NL*/ #define DEBUGMODE_LMP_FIX_NEIGHLIST_MESH false //(update->ntimestep>15400 && comm->me ==1)
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID 0
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID 60

FixNeighlistMesh::FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp,narg,arg),
  buildNeighList(false),
  movingMesh(false)
{
    if(!modify->find_fix_id(arg[3]) || strcmp(modify->find_fix_id(arg[3])->style,"mesh/surface"))
        error->fix_error(FLERR,this,"illegal caller");

    caller_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[3]));
    mesh_ = caller_->triMesh();
}

/* ---------------------------------------------------------------------- */

FixNeighlistMesh::~FixNeighlistMesh()
{
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_create()
{
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::initializeNeighlist()
{
    // remove old lists, init new ones
    //NP need to initialize values
    //NP loop over owned + ghost tris as ghosts may extend into my subbox
    //NP is called via fix mesh setup()

    for(int i = 0; i < numContacts.size(); i++)
        numContacts.del(i);
    for(int i = 0; i < contactList.size(); i++)
        contactList.del(i);

    int nall = mesh_->sizeLocal()+mesh_->sizeGhost();

    for(int iTri = 0; iTri < nall; iTri++)
        numContacts.add(0);

    contactList.add(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::setup_pre_force(int foo)
{
    //NP initial tri-sphere neighlist build
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_delete(bool unfixflag)
{

}

/* ---------------------------------------------------------------------- */

int FixNeighlistMesh::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_neighbor()
{
    buildNeighList = true;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_force(int vflag)
{
    if(!buildNeighList) return;

    movingMesh = mesh_->isMoving();

    /*NL*/ //fprintf(screen,"building neighbor list at timestep %d\n",update->ntimestep);

    buildNeighList = false;

    contactList.empty();
    numContacts.empty();

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    if(movingMesh)
    {
      skin = neighbor->skin;
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    //NP from neighbor; get the binning
    mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;

    int nall = mesh_->sizeLocal() + mesh_->sizeGhost();


    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID <= atom->get_map_size() && update->ntimestep > 0 &&
    /*NL*/      atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID) >= 0)
    /*NL*/ {
    /*NL*/          if(!r) error->one(FLERR,"debugmode not made for SPH");
    /*NL*/          int iTriDeb = mesh_->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID);
    /*NL*/          int iAtomDeb = atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID);
    /*NL*/          int ixDeb, iyDeb, izDeb;
    /*NL*/          int iBinDeb = neighbor->coord2bin(atom->x[iAtomDeb],ixDeb, iyDeb, izDeb);
    /*NL*/          fprintf(screen, "**step %d, particle id %d at bin %d (indixes %d %d %d) on proc %d, within skin to target tri %s\n",
    /*NL*/                      update->ntimestep,DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID,
    /*NL*/                      iBinDeb,ixDeb, iyDeb, izDeb,comm->me,
    /*NL*/                      mesh_->resolveTriSphereNeighbuild(iTriDeb,atom->radius[iAtomDeb],atom->x[iAtomDeb],skin) ? "true" : "false" );
    /*NL*/ }

    for(int iTri = 0; iTri < nall; iTri++)
      handleTriangle(iTri);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::handleTriangle(int iTri)
{
    // get bounding box of element on this subdomain
    //NP is bounded by sublo, subhi
    BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);

    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    double lo[3],hi[3];

    // extend bbox by cutneighmax and get bin boundaries
    getBinBoundariesFromBoundingBox(b,ixMin,ixMax,iyMin,iyMax,izMin,izMax);

    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
    /*NL*/ {
    /*NL*/          double lo[3],hi[3];
    /*NL*/          b.getBoxBounds(lo,hi);
    /*NL*/          vectorScalarSubtract3D(lo,distmax);
    /*NL*/          vectorScalarAdd3D(hi,distmax);
    /*NL*/          int ixDebLo, iyDebLo, izDebLo, ixDebHi, iyDebHi, izDebHi;
    /*NL*/          neighbor->coord2bin(lo,ixDebLo, iyDebLo, izDebLo);
    /*NL*/          neighbor->coord2bin(hi,ixDebHi, iyDebHi, izDebHi);
    /*NL*/          fprintf(screen, "handleTriangle for tri id %d on proc %d, indices from here  are %d %d //  %d %d //  %d %d , bbox is %f %f //  %f %f //  %f %f  \n",
    /*NL*/                      mesh_->id(iTri),comm->me,ixMin,ixMax,iyMin,iyMax,izMin,izMax,lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
    /*NL*/          fprintf(screen, "handleTriangle for tri id %d on proc %d, indices from neigh are %d %d //  %d %d //  %d %d , bbox is %f %f //  %f %f //  %f %f  \n",
    /*NL*/                      mesh_->id(iTri),comm->me,ixDebLo,ixDebHi,iyDebLo,iyDebHi,izDebLo,izDebHi,lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]);
    /*NL*/ }

    int numContTmp = 0;

    // only do this if I own particles
    if(nlocal)
    {
        for(int ix=ixMin;ix<=ixMax;ix++)
          for(int iy=iyMin;iy<=iyMax;iy++)
            for(int iz=izMin;iz<=izMax;iz++)
            {
              int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;

              /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
              /*NL*/          fprintf(screen, "       handleTriangle tri id %d on proc %d - checking bin %d\n",
              /*NL*/                      mesh_->id(iTri),comm->me, iBin);

              int iAtom = binhead[iBin];
              while(iAtom != -1 && iAtom < nlocal)
              {
                /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri)
                /*NL*/                                     && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID == atom->map(iAtom) )
                /*NL*/          fprintf(screen, "   handleTriangle atom %d for tri id %d on proc %d\n",
                /*NL*/                      atom->map(iAtom),mesh_->id(iTri),comm->me);
                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom] : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  //NP include iAtom in neighbor list
                  numContTmp++;
                  contactList.add(iAtom);
                }
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
              }
            }
    }

    numContacts.add(numContTmp);
    //printf("numTri %d numContacts %d\n",iTri,numContTmp);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesFromBoundingBox(BoundingBox &b,
      int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
    double delta = distmax;
    double tri_xmin[3] = {b.xLo-delta,b.yLo-delta,b.zLo-delta};
    double tri_xmax[3] = {b.xHi+delta,b.yHi+delta,b.zHi+delta};

    int binmin = neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    /*NP OLD
    ixMin = (binmin % (mbiny*mbinx)) % mbinx;
    if(ixMin < 0) ixMin = 0;
    iyMin = static_cast<int>(round(static_cast<double>(((binmin - ixMin) % (mbiny*mbinx)))/static_cast<double>(mbinx)));
    if(iyMin < 0) iyMin = 0;
    izMin = static_cast<int>(round(static_cast<double>((binmin - ixMin - iyMin*mbinx))/static_cast<double>((mbiny*mbinx))));
    if(izMin < 0) izMin = 0;
    */

    int binmax= neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    /*NP OLD
    ixMax = (binmax % (mbiny*mbinx)) % mbinx;
    if(ixMax > mbinx-1) ixMax = mbinx-1;
    iyMax = static_cast<int>(round(static_cast<double>(((binmax - ixMax) % (mbiny*mbinx)))/static_cast<double>(mbinx)));
    if(iyMax > mbiny-1) iyMax = mbiny-1;
    izMax = static_cast<int>(round(static_cast<double>((binmax - ixMax - iyMax*mbinx))/static_cast<double>((mbiny*mbinx))));
    if(izMax > mbinz-1) izMax = mbinz-1;
    */
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getPointers(int *&cList, int *&nContact)
{
    cList = contactList.begin();
    nContact = numContacts.begin();
}
