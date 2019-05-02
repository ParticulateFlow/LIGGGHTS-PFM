/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   Copyright 2014-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Philippe Seil (JKU)
   Daniel Queteschiner (JKU)
------------------------------------------------------------------------- */


#include "region_distance_field.h"

#include <math.h>
#include <assert.h>

#include "region.h"
#include "vector_liggghts.h"



namespace LIGGGHTS {

  RegionDistanceField::RegionDistanceField()
    : nx(0),ny(0),nz(0),dx(0.)
  {
    LAMMPS_NS::vectorZeroize3D(xlo);
    LAMMPS_NS::vectorZeroize3D(xhi);
  }

  void RegionDistanceField::build(LAMMPS_NS::Region *region, LAMMPS_NS::BoundingBox &reg_bbox, double const rmax)
  {
    bbox = reg_bbox;
    double bbox_xlo[3],bbox_xhi[3],bbox_extent[3];
    reg_bbox.getBoxBounds(bbox_xlo,bbox_xhi);
    reg_bbox.getExtent(bbox_extent);

    dx = 2.*rmax;
    oneoverdx = 1./dx;
    // Note: test_rad = bounding sphere of bin plus particle radius
    // if the sphere with test_rad radius is completely inside the
    // insertion region, any particle center inside the corresponding
    // bin is also completely inside the insertion region
    test_rad = sqrt(3.0)*rmax+rmax;

    nx = bbox_extent[0]*oneoverdx;
    ny = bbox_extent[1]*oneoverdx;
    nz = bbox_extent[2]*oneoverdx;

    double const x_over = bbox_extent[0] - dx*static_cast<double>(nx);
    x0[0] = bbox_xlo[0] + 0.5*x_over - 0.5*dx;

    double const y_over = bbox_extent[1] - dx*static_cast<double>(ny);
    x0[1] = bbox_xlo[1] + 0.5*y_over - 0.5*dx;

    double const z_over = bbox_extent[2] - dx*static_cast<double>(nz);
    x0[2] = bbox_xlo[2] + 0.5*z_over - 0.5*dx;

    data.insert(data.begin(),nx*ny*nz,OUTSIDE);

    double pos_tmp[3];
    for(int i=0;i<nx;i++){
      for(int j=0;j<ny;j++){
        for(int k=0;k<nz;k++){
          index3ToPos(i,j,k,pos_tmp);

          if(!region->match(pos_tmp[0],pos_tmp[1],pos_tmp[2]))
            continue;

          int const index = index3ToIndex1(i,j,k);
          pos_tmp[0] += 0.5*dx;
          pos_tmp[1] += 0.5*dx;
          pos_tmp[2] += 0.5*dx;

          if(region->match_shrinkby_cut(pos_tmp,test_rad))
            data[index] = INSIDE;
          else
            data[index] = BOUNDARY;
        }
      }
    }
  }

  void RegionDistanceField::reset()
  {
    data.clear();
  }

  bool RegionDistanceField::isInside(double *x)
  {
    int const index = posToIndex(x);

    if(index<0)
      return false;

    return data[index] == INSIDE;
  }

  bool RegionDistanceField::isOutside(double *x)
  {
    int const index = posToIndex(x);

    if(index<0)
      return false;

    return data[index] == OUTSIDE;
  }

  bool RegionDistanceField::isInBoundary(double *x)
  {
    int const index = posToIndex(x);

    if(index<0)
      return false;

    return data[index] == BOUNDARY;
  }


  int RegionDistanceField::index3ToIndex1(int const ix, int const iy, int const iz)
  {
    return (ix + nx*iy + nx*ny*iz);
  }

  int RegionDistanceField::posToIndex(double *x)
  {
    if(!bbox.isInside(x)) return -1;

    int const ix = (x[0]-x0[0])*oneoverdx;
    if(ix < 0 || ix > nx-1) return -1;

    int const iy = (x[1]-x0[1])*oneoverdx;
    if(iy < 0 || iy > ny-1) return -1;

    int const iz = (x[2]-x0[2])*oneoverdx;
    if(iz < 0 || iz > nz-1) return -1;

    return index3ToIndex1(ix,iy,iz);
  }

  void RegionDistanceField::indexToPos(int index, double *x)
  {
    LAMMPS_NS::vectorCopy3D(x0,x);

    int const iz = index / (nx*ny);
    index -= iz*nx*ny;

    int const iy = index/nx;
    int const ix = index % nx;

    x[0] += static_cast<double>(ix)*dx;
    x[1] += static_cast<double>(iy)*dx;
    x[2] += static_cast<double>(iz)*dx;
  }

  void RegionDistanceField::index3ToPos(int ix, int iy, int iz, double *x)
  {
    x[0] = x0[0] + static_cast<double>(ix)*dx;
    x[1] = x0[1] + static_cast<double>(iy)*dx;
    x[2] = x0[2] + static_cast<double>(iz)*dx;
  }


} // namespace LIGGGHTS
