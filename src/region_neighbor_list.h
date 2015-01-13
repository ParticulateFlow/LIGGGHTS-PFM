/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Richard Berger <richard.berger@jku.at> (JKU Linz)
------------------------------------------------------------------------- */
#ifndef REGION_NEIGHBOR_LIST_H
#define REGION_NEIGHBOR_LIST_H

#include <vector>
#include "vector_liggghts.h"
#include "bounding_box.h"

namespace LIGGGHTS {

/**
 * @brief A small particle structure
 */
struct Particle {
  double x[3];
  double radius;

  Particle(double * pos, double rad) {
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
  }
};


/**
 * @brief A neighbor list of of a certain region
 *
 * This implementation uses the same binning approach used in LAMMPS neighbor lists.
 * Instead of accessing internal data structures directly, manipulations and queries
 * can only occur through the given interface.
 */
class RegionNeighborList
{
  typedef std::vector<Particle> ParticleBin;

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

public:
    RegionNeighborList();

    bool hasOverlap(double * x, double radius) const;
    void insert(double * x, double radius);
    size_t count() const;
    void reset();
    bool setBoundingBox(LAMMPS_NS::BoundingBox & bb, double maxrad);
};

}

#endif // REGION_NEIGHBOR_LIST_H
