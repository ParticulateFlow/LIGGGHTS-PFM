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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef REGION_NEIGHBOR_LIST_H
#define REGION_NEIGHBOR_LIST_H

#include <vector>
#include "vector_liggghts.h"
#include "bounding_box.h"
#include "pointers.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_superquadric.h"
#endif

namespace LIGGGHTS {

/**
 * @brief A small particle structure
 */
struct Particle {
  double x[3];
  double radius;
  int type;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double shape[3];
  double quaternion[4];
  double blockiness[2];
#endif

  Particle(double * pos, double rad, int type) {
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
    this->type = type;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    quaternion[0] = 1.0;
    quaternion[1] = quaternion[2] = quaternion[3] = 0.0;
    shape[0] = shape[1] = shape[2] = radius;
    blockiness[0] = blockiness[1] = 2.0;
#endif
  }

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  Particle(int _i,double * pos, double rad, int type, double *quaternion_, double *shape_, double *blockiness_, int,int,double,double,double) {
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
    this->type = type;
    LAMMPS_NS::vectorCopy4D(quaternion_, quaternion);
    LAMMPS_NS::vectorCopy3D(shape_, shape);
    LAMMPS_NS::vectorCopy2D(blockiness_, blockiness);
  }
#endif
};


/**
 * @brief A neighbor list of of a certain region
 *
 * This implementation uses the same binning approach used in LAMMPS neighbor lists.
 * Instead of accessing internal data structures directly, manipulations and queries
 * can only occur through the given interface.
 */
class RegionNeighborList : protected LAMMPS_NS::Pointers
{
public:
  typedef std::vector<Particle> ParticleBin;

private:
  std::vector<ParticleBin> bins;  // list of particle bins
  std::vector<int> stencil;       // stencil used to check bins for collisions
  size_t ncount;                  // total number of particles in neighbor list

  double bboxlo[3];               // lowest point of bounding box
  double bboxhi[3];               // highest point of bounding box

  int nbinx,nbiny,nbinz;          // # of global bins
  int mbinx,mbiny,mbinz;          // # of bins in each dimension
  int mbinxlo,mbinylo,mbinzlo;    // offsets of (0,0,0) bin

  double binsizex,binsizey,binsizez;  // bin sizes
  double bininvx,bininvy,bininvz;     // inverse of bin sizes

  double bin_distance(int i, int j, int k);
  int coord2bin(double *x) const;
  bool type_exclusion(int itype, int jtype) const;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  int check_obb_flag;
#endif

public:
    RegionNeighborList(LAMMPS_NS::LAMMPS *lmp);

    bool hasOverlap(Particle &p) const;
    bool hasOverlap(double * x, double radius, int type=-1) const;
    void insert(double * x, double radius, int type);
    size_t count() const;
    void reset();
    bool setBoundingBox(LAMMPS_NS::BoundingBox & bb, double maxrad);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    inline void coord2bin_calc_interpolation_weights(double *x,int ibin,int ix,int iy, int iz,int &quadrant,double &wx,double &wy,double &wz) const;
    int coord2bin(double *x,int &quadrant,double &wx,double &wy,double &wz) const;
    bool hasOverlap_superquadric(double * x, double radius, int type, double *quaternion, double *shape, double *blockiness) const;
    void insert_superquadric(double * x, double radius, int type, double *quaternion, double *shape, double *blockiness, int index = -1);
    void set_obb_flag(int check_obb_flag_) { check_obb_flag = check_obb_flag_; }
    int mbins() const { return mbinx*mbiny*mbinz; }
#endif

    ParticleBin* getParticlesCloseTo(double *x, double cutoff);
};

}

#endif // REGION_NEIGHBOR_LIST_H
