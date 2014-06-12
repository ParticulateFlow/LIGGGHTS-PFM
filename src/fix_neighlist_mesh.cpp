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
   Richard Berger (JKU Linz)
   Philippe Seil (JKU Linz)
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface.h"
#include "fix_property_atom.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "atom.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>
#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

/*NL*/ #define DEBUGMODE_LMP_FIX_NEIGHLIST_MESH false //(update->ntimestep>15400 && comm->me ==1)
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID 0
/*NL*/ #define DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID 60

FixNeighlistMesh::FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp,narg,arg),
  fix_nneighs_(0),
  fix_nneighs_name_(0),
  buildNeighList(false),
  numAllContacts_(0),
  globalNumAllContacts_(false),
  mbinx(0),
  mbiny(0),
  mbinz(0),
  maxhead(0),
  bins(NULL),
  binhead(NULL),
  skin(0.0),
  distmax(0.0),
  x(NULL),
  r(NULL),
  changingMesh(false),
  changingDomain(false),
  last_bin_update(-1)
{
    if(!modify->find_fix_id(arg[3]) || !dynamic_cast<FixMeshSurface*>(modify->find_fix_id(arg[3])))
        error->fix_error(FLERR,this,"illegal caller");

    caller_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[3]));
    mesh_ = caller_->triMesh();
}

/* ---------------------------------------------------------------------- */

FixNeighlistMesh::~FixNeighlistMesh()
{
    delete [] fix_nneighs_name_;
    last_bin_update = -1;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_create()
{
    // register
    if(!fix_nneighs_)
    {
        const char* fixarg[9];
        delete [] fix_nneighs_name_;
        fix_nneighs_name_ = new char[strlen(mesh_->mesh_id())+1+14];
        sprintf(fix_nneighs_name_,"n_neighs_mesh_%s",mesh_->mesh_id());

        fixarg[0]=fix_nneighs_name_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=fix_nneighs_name_;
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart - REQUIRED!
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_nneighs_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);

        //NP this is to suppress that set_arrays is called in modify->setup
        //NP which would erase data that is written in setup_pre_force() in this class
        fix_nneighs_->just_created = false;
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::initializeNeighlist()
{
    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    // remove old lists, init new ones
    //NP need to initialize values
    //NP loop over owned + ghost tris as ghosts may extend into my subbox
    //NP is called via fix mesh setup()

    const size_t nall = mesh_->sizeLocal()+mesh_->sizeGhost();

    while(triangles.size() > nall) {
        triangles.pop_back();
    }

    while(triangles.size() < nall) {
        triangles.push_back(TriangleNeighlist());
    }

    for(size_t iTri = 0; iTri < nall; iTri++) {
        TriangleNeighlist & triangle = triangles[iTri];
        triangle.contacts.reserve(std::max(triangle.contacts.capacity(), static_cast<size_t>(128)));
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::setup_pre_force(int foo)
{
    //NP initial tri-sphere neighlist build
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_setup_pre_force(int foo)
{
    //NP initial tri-sphere neighlist build
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_delete(bool unfixflag)
{
    if(unfixflag)
    {
        modify->delete_fix(fix_nneighs_->id);
    }
}

/* ---------------------------------------------------------------------- */

int FixNeighlistMesh::setmask()
{
    int mask = 0;
    mask |= MIN_PRE_NEIGHBOR;
    mask |= PRE_NEIGHBOR;
    mask |= MIN_PRE_FORCE;
    mask |= PRE_FORCE;
    mask |= POST_RUN;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_neighbor()
{
    buildNeighList = true;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_pre_force(int vflag)
{
    pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

//NP this is called before FixContactHistoryMesh::pre_force()
//NP because of the order of creation in FixWallGran::post_create()

void FixNeighlistMesh::pre_force(int)
{
    if(!buildNeighList) return;

    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    /*NL*/ //fprintf(screen,"***building neighbor list at timestep "BIGINT_FORMAT"\n",update->ntimestep);

    buildNeighList = false;
    numAllContacts_ = 0;

    // set num_neigh = 0
    memset(fix_nneighs_->vector_atom, 0, atom->nlocal*sizeof(double));

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    //NP cutneighmax includes contactDistanceFactor, thus rmax includes this as well
    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    double prev_skin = skin;
    double prev_distmax = distmax;

    if(changingMesh)
    {
      skin = neighbor->skin;
      //NP cutneighmax includes contactDistanceFactor, thus distmax includes this as well
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      //NP cutneighmax includes contactDistanceFactor, thus distmax includes this as well
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    //NP from neighbor; get the binning
    mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;

    const size_t nall = mesh_->sizeLocal() + mesh_->sizeGhost();

    // update cache if necessary
    if (triangles.size() != nall) {
      initializeNeighlist();
    }

    // update precomputed bins if necessary
    if((skin != prev_skin) || (distmax != prev_distmax) || (neighbor->last_setup_bins_timestep > last_bin_update)) {
      generate_bin_list(nall);
    }

    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID <= atom->get_map_size() && update->ntimestep > 0 &&
    /*NL*/      atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID) >= 0)
    /*NL*/ {
    /*NL*/          if(!r) error->one(FLERR,"debugmode not made for SPH");
    /*NL*/          int iTriDeb = mesh_->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID);
    /*NL*/          int iAtomDeb = atom->map(DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID);
    /*NL*/          int ixDeb, iyDeb, izDeb;
    /*NL*/          int iBinDeb = neighbor->coord2bin(atom->x[iAtomDeb],ixDeb, iyDeb, izDeb);
    /*NL*/          fprintf(screen, "**step "BIGINT_FORMAT", particle id %d at bin %d (indixes %d %d %d) on proc %d, within skin to target tri %s\n",
    /*NL*/                      update->ntimestep,DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID,
    /*NL*/                      iBinDeb,ixDeb, iyDeb, izDeb,comm->me,
    /*NL*/                      mesh_->resolveTriSphereNeighbuild(iTriDeb,atom->radius[iAtomDeb]*neighbor->contactDistanceFactor,atom->x[iAtomDeb],skin) ? "true" : "false" );
    /*NL*/ }

    for(size_t iTri = 0; iTri < nall; iTri++) {
      TriangleNeighlist & triangle = triangles[iTri];
      handleTriangle(iTri);
      numAllContacts_ += triangle.contacts.size();
    }

    /*NL*/ //if(nall > 0) fprintf(screen,"size numContactsSum %d numAllContacts_ %d vs. %d\n",numContactsSum.size(),numAllContacts_,numContacts(nall-1)+numContactsSum(nall-1));

    if(globalNumAllContacts_) {
      MPI_Sum_Scalar(numAllContacts_,world);
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::handleTriangle(int iTri)
{
    TriangleNeighlist & triangle = triangles[iTri];
    std::vector<int> & neighbors = triangle.contacts;
    int & nchecked = triangle.nchecked;
    int *mask = atom->mask;
    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    double contactDistanceFactor = neighbor->contactDistanceFactor;

    neighbors.clear();

    nchecked = 0;

    // only do this if I own particles
    if(nlocal)
    {
      if(changingMesh || changingDomain)
      {
        getBinBoundariesForTriangle(iTri,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
    /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
    /*NL*/ {
    /*NL*/          double lo[3],hi[3];
    /*NL*/          //b.getBoxBounds(lo,hi);
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

        for(int ix=ixMin;ix<=ixMax;ix++) {
          for(int iy=iyMin;iy<=iyMax;iy++) {
            for(int iz=izMin;iz<=izMax;iz++) {
              int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;

              /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
              /*NL*/          fprintf(screen, "       handleTriangle tri id %d on proc %d - checking bin %d\n",
              /*NL*/                      mesh_->id(iTri),comm->me, iBin);

              int iAtom = binhead[iBin];
              //NP only handle local atoms
              while(iAtom != -1 && iAtom < nlocal)
              {
                if(! (mask[iAtom] & groupbit))
                {
                    if(bins) iAtom = bins[iAtom];
                    else iAtom = -1;
                    continue;
                }
                nchecked++;

                /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri)
                /*NL*/                                     && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID == atom->map(iAtom) )
                /*NL*/          fprintf(screen, "   handleTriangle atom %d for tri id %d on proc %d\n",
                /*NL*/                      atom->map(iAtom),mesh_->id(iTri),comm->me);
                if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
                {
                  //NP include iAtom in neighbor list
                  neighbors.push_back(iAtom);
                  fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
                  /*NL*/ //if(377==atom->tag[iAtom]) fprintf(screen,"proc %d, step "BIGINT_FORMAT" adding pair tri tag %d atom ta %d to NEIGHLIST\n",comm->me,update->ntimestep,mesh_->id(iTri),atom->tag[iAtom]);
                }
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
              }
            }
          }
        }
      } else {
        const std::vector<int> & triangleBins = triangle.bins;
        const int bincount = triangleBins.size();
        for(int i = 0; i < bincount; i++) {
          const int iBin = triangleBins[i];
          /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri))
          /*NL*/          fprintf(screen, "       handleTriangle tri id %d on proc %d - checking bin %d\n",
          /*NL*/                      mesh_->id(iTri),comm->me, iBin);

          int iAtom = binhead[iBin];
          while(iAtom != -1 && iAtom < nlocal)
          {
            if(! (mask[iAtom] & groupbit))
            {
                if(bins) iAtom = bins[iAtom];
                else iAtom = -1;
                continue;
            }
            nchecked++;

            /*NL*/ if(DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && DEBUG_LMP_FIX_NEIGHLIST_MESH_M_ID == mesh_->id(iTri)
            /*NL*/                                     && DEBUG_LMP_FIX_NEIGHLIST_MESH_P_ID == atom->map(iAtom) )
            /*NL*/          fprintf(screen, "   handleTriangle atom %d for tri id %d on proc %d\n",
            /*NL*/                      atom->map(iAtom),mesh_->id(iTri),comm->me);
            if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
            {
              //NP include iAtom in neighbor list
              neighbors.push_back(iAtom);
              fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
            }
            if(bins) iAtom = bins[iAtom];
            else iAtom = -1;
          }
        }
      }
    }

    /*NL*/// fprintf(screen,"iTri %d numContacts %d\n",iTri, neighbors.size());
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesFromBoundingBox(BoundingBox &b,
      int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
    double delta = distmax;
    double tri_xmin[3] = {b.xLo-delta,b.yLo-delta,b.zLo-delta};
    double tri_xmax[3] = {b.xHi+delta,b.yHi+delta,b.zHi+delta};

    //int binmin = neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    /*NP OLD
    ixMin = (binmin % (mbiny*mbinx)) % mbinx;
    if(ixMin < 0) ixMin = 0;
    iyMin = static_cast<int>(round(static_cast<double>(((binmin - ixMin) % (mbiny*mbinx)))/static_cast<double>(mbinx)));
    if(iyMin < 0) iyMin = 0;
    izMin = static_cast<int>(round(static_cast<double>((binmin - ixMin - iyMin*mbinx))/static_cast<double>((mbiny*mbinx))));
    if(izMin < 0) izMin = 0;
    */

    //int binmax= neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
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

void FixNeighlistMesh::getBinBoundariesForTriangle(int iTri, int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
  // disable optimization for movingMesh or shrink-wrapped domain
  if(changingMesh || changingDomain) {
    BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);
    // extend bbox by cutneighmax and get bin boundaries
    getBinBoundariesFromBoundingBox(b,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
  } else {
    // use cached boundary information
    const BinBoundary & b = triangles[iTri].boundary;
    ixMin = b.xlo;
    ixMax = b.xhi;
    iyMin = b.ylo;
    iyMax = b.yhi;
    izMin = b.zlo;
    izMax = b.zhi;
  }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_run()
{
  last_bin_update = -1; // reset binning for possible next run
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::generate_bin_list(size_t nall)
{
  // precompute triangle bin boundaries
  // disable optimization for changing mesh or domain
  if (!(changingMesh || changingDomain)) {
    double dx = neighbor->binsizex / 2.0;
    double dy = neighbor->binsizey / 2.0;
    double dz = neighbor->binsizez / 2.0;
    double maxdiag = sqrt(dx * dx + dy * dy + dz * dz);

    for (size_t iTri = 0; iTri < nall; iTri++) {
      TriangleNeighlist & triangle = triangles[iTri];
      std::vector<int> & binlist = triangle.bins;
      binlist.clear();

      BinBoundary& bb = triangle.boundary;
      BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);

      // extend bbox by cutneighmax and get bin boundaries
      getBinBoundariesFromBoundingBox(b, bb.xlo, bb.xhi, bb.ylo, bb.yhi, bb.zlo, bb.zhi);

      // look at bins and exclude unnecessary ones
      double center[3];
      int total = 0;
      for (int ix = bb.xlo; ix <= bb.xhi; ix++) {
        for (int iy = bb.ylo; iy <= bb.yhi; iy++) {
          for (int iz = bb.zlo; iz <= bb.zhi; iz++) {
            int iBin = iz * mbiny * mbinx + iy * mbinx + ix;
            if (iBin < 0 || iBin >= maxhead)
              continue;

            // determine center of bin (ix, iy, iz)
            neighbor->bin_center(ix, iy, iz, center);

            if (mesh_->resolveTriSphereNeighbuild(iTri, maxdiag, center, distmax + skin))
            {
              binlist.push_back(iBin);
            }
            total++;
          }
        }
      }
      /*NL*/ if (DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && comm->me == 0) fprintf(screen, "triangle %lu bins: %lu / %d\n", iTri, binlist.size(), total);
      /*NL*/ if (DEBUGMODE_LMP_FIX_NEIGHLIST_MESH && comm->me == 0) fprintf(logfile, "triangle %lu bins: %lu / %d\n", iTri, binlist.size(), total);
    }
  }

  last_bin_update = update->ntimestep;
}

int FixNeighlistMesh::getSizeNumContacts()
{
  return mesh_->sizeLocal() + mesh_->sizeGhost();
}
