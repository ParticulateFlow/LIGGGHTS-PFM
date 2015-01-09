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
   Contributing authors:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "particleSpatialDistribution.h"
#include "random_park.h"
#include "atom.h"
#include "primitive_wall.h"
#include "tri_mesh.h"
#include "tri_mesh_contacts.h"

using namespace LAMMPS_NS;

typedef struct
{
  double x;
  double y;
  double z;
} point;

ParticleSpatialDistribution::ParticleSpatialDistribution(RanPark *rp, double overlap, int maxattempt) :
  RNG(rp),
  max_overlap(overlap),
  maxattempt(maxattempt)
{
}


ParticleSpatialDistribution::~ParticleSpatialDistribution()
{
}


bool ParticleSpatialDistribution::isPointInSphere(const std::vector<double> &center, double radius, const std::vector<double> &x, double *dir, double *dist)
{
  // dir points from x to center
  double delta[3];
  delta[0] = center[0] - x[0];
  delta[1] = center[1] - x[1];
  delta[2] = center[2] - x[2];

  double dist_sq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];

  if (dir || dist) {
    if (dist) {
      double _dist = sqrt(dist_sq);
      *dist = _dist;
    }

    if (dir) {
      dir[0] = delta[0];
      dir[1] = delta[1];
      dir[2] = delta[2];
    }
  }

  if (dist_sq > radius*radius)
    return false;

  return true;
}


void ParticleSpatialDistribution::randomPointInSphere(double radius, std::vector<double> &x)
{
  x[0] = RNG->gaussian();
  x[1] = RNG->gaussian();
  x[2] = RNG->gaussian();

  double length_sq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];

  if (length_sq > 0.0) {
    const double inverse_length = 1.0/sqrt(length_sq);
    const double factor =  RNG->uniform() * radius * inverse_length;
    x[0] *= factor;
    x[1] *= factor;
    x[2] *= factor;
  }
}

#define maxGen 5
#define minRad 0.01
#define clans false

bool create(int generation, int clan, double A, double B, double C, double D, double E, std::multimap<double, point>& pos)
{
  // returns false if the radius is less than minRad
  // Convert special coords to inversive coords
  double rt2o2 = sqrt(2.)/2.;
  double a = (   -B + C + D - E) * rt2o2;
  double b = (    B - C + D - E) * rt2o2;
  double c = (    B + C - D - E) * rt2o2;
  double d = (A - B - C - D - E);
  double e = (    B + C + D + E) * sqrt(6.)/2.;

  // Convert inversive coords to Cartesian coords
  // Positive radius required for drawing

  double r = 0.05;
  if (e > d) {
    r = 1./(e - d);
  } else {
    r = 1./(d - e);
  }

  /*if(r < minRad) {
    return false;
  } else*/
  {
    bool inversionCircle = generation < 0;
    bool outerCircle = (generation == 0 && clan == 0);

    if (inversionCircle || outerCircle) {
      return false;
    } else {
      point x = {a * r, b * r, c * r};
      pos.insert(std::pair<double,point>(r,x));
      return true;
    }
  }
}

double generateItem(int item, double itemVal, int i, double A, double B, double C, double D, double E)
{
  switch(i) {
  case 0:
    if(item == i) return -A;
    else return A + itemVal;
  case 1:
    if(item == i) return -B;
    else return B + itemVal;
  case 2:
    if(item == i) return -C;
    else return C + itemVal;
  case 3:
    if(item == i) return -D;
    else return D + itemVal;
  default:
    if(item == i) return -E;
    else return E + itemVal;
  }
}


void generate(int generation, int clan, int i, double A, double B, double C, double D, double E, std::multimap<double, point>& pos)
{
  double A1 = generateItem(0, A, i, A, B, C, D, E);
  double B1 = generateItem(1, B, i, A, B, C, D, E);
  double C1 = generateItem(2, C, i, A, B, C, D, E);
  double D1 = generateItem(3, D, i, A, B, C, D, E);
  double E1 = generateItem(4, E, i, A, B, C, D, E);

  // avoid duplicates
  if ((i == 0 && (B1 < A1 || C1 < A1 || D1 < A1 || E1 < A1)) ||
      (i == 1 && (A1 <=B1 || C1 < B1 || D1 < B1 || E1 < B1)) ||
      (i == 2 && (A1 <=C1 || B1 <=C1 || D1 < C1 || E1 < C1)) ||
      (i == 3 && (A1 <=D1 || B1 <=D1 || C1 <=D1 || E1 < D1)) ||
      (i == 4 && (A1 <=E1 || B1 <=E1 || C1 <=E1 || D1 <=E1))) {
  } else {
    if (create(generation, clan, A1, B1, C1, D1, E1, pos)) {
      if (generation < maxGen) {
        for (int j = 0; j < 5; ++j) {
          if (j != i) {
            generate((generation+1), clan, j, A1, B1, C1, D1, E1, pos);
          }
        }
      }
    }
  }
}


void ParticleSpatialDistribution::apollonianInsertion(
        double radius,
        const std::vector<double>& radii,
        std::vector<double>& x,
        std::vector<double>& y,
        std::vector<double>& z)
{
  const unsigned int nParticles = radii.size();
  if(nParticles == 0) return;

  x.clear();
  y.clear();
  z.clear();

  const double r1 = 1.0;
  std::multimap<double, point> pos;
  // create initial spheres and generate all spheres of their clans
  create(0, 0, r1, 0., 0., 0., 0., pos);
  create(0, 1, 0., r1, 0., 0., 0., pos);
  create(0, 2, 0., 0., r1, 0., 0., pos);
  create(0, 3, 0., 0., 0., r1, 0., pos);
  create(0, 4, 0., 0., 0., 0., r1, pos);
  generate(1, 0, 0, r1, 0., 0., 0., 0., pos);
  generate(1, 1, 1, 0., r1, 0., 0., 0., pos);
  generate(1, 2, 2, 0., 0., r1, 0., 0., pos);
  generate(1, 3, 3, 0., 0., 0., r1, 0., pos);
  generate(1, 4, 4, 0., 0., 0., 0., r1, pos);

  unsigned int idx = 0;
  for (std::multimap<double,point>::reverse_iterator rit = pos.rbegin(); rit != pos.rend(); ++rit, ++idx) {
    x.push_back(rit->second.x*radius);
    y.push_back(rit->second.y*radius);
    z.push_back(rit->second.z*radius);
  }

  while (x.size() < nParticles) {
    x.push_back(0.0);
    y.push_back(0.0);
    z.push_back(0.0);
  }
}


void ParticleSpatialDistribution::randomInsertion(
        double *xp,
        double radius,
        const std::vector<double>& radii,
        std::vector<std::vector<double> >& x,
        const std::vector<std::vector<double> >& ext_atoms,
        const std::vector<PrimitiveWall*>& prim_walls,
        const std::vector<TriMeshContacts*>& meshes)
{
  const unsigned int nParticles = radii.size();
  const unsigned int ext = ext_atoms.size();

  if (nParticles == 0) return;
  x.resize(nParticles);

  for (unsigned int idx = 0; idx < nParticles; ++idx) {
    std::vector<double> x_rand(3, 0.0);
    std::vector<double> most_acceptable_x(4, 0.0);
    bool centerInOccupiedSpace;
    int ntry = 0;

    do {
      randomPointInSphere(radius-radii[idx], x_rand);
      centerInOccupiedSpace = false;
      ++ntry;

      // check if point is inside any external sphere
      for (unsigned int ext_i = 0; ext_i < ext; ++ext_i) {
        std::vector<double> ext_c(3, 0.0);
        ext_c[0] = ext_atoms[ext_i][0] - xp[0];
        ext_c[1] = ext_atoms[ext_i][1] - xp[1];
        ext_c[2] = ext_atoms[ext_i][2] - xp[2];

        centerInOccupiedSpace = isPointInSphere(ext_c, ext_atoms[ext_i][3]+radii[idx], x_rand);
        if (centerInOccupiedSpace)
          break;
      }

      // check if point is inside a primitive wall
      if (!centerInOccupiedSpace) {
        double x_global[3] = {x_rand[0] + xp[0], x_rand[1] + xp[1],x_rand[2] + xp[2]};
        for (std::vector<PrimitiveWall*>::const_iterator it = prim_walls.begin(); it != prim_walls.end(); ++it) {
          double delta[3];
          centerInOccupiedSpace = (!((*it)->resolveSameSide(xp, x_global)) || ((*it)->resolveContact(x_global, radii[idx], delta) < 0.0));
          if (centerInOccupiedSpace)
            break;
        }
      }

      // check if point is inside a mesh wall
      if (!centerInOccupiedSpace) {
        double x_global[3] = {x_rand[0] + xp[0], x_rand[1] + xp[1],x_rand[2] + xp[2]};
        for (std::vector<TriMeshContacts*>::const_iterator it = meshes.begin(); it != meshes.end() && !centerInOccupiedSpace; ++it) {
          bool sameSide = false;
          double delta[3];
          for (std::set<int>::const_iterator it2 = (*it)->contacts.begin(); it2 != (*it)->contacts.end(); ++it2) {
            double deltan = (*it)->mesh_->resolveTriSphereContact(0, *it2, radii[idx], x_global, delta);
            if (deltan < 0.0) {
              centerInOccupiedSpace = true;
              break;
            } else if (!sameSide) {
              // only if particles is smaller than overlap it may be on the wrong side and not collide
              if (2.0*radii[idx] > max_overlap || (*it)->mesh_->resolveSameSide(*it2, x_global, xp)) {
                sameSide = true;
              }
            }
          }
          centerInOccupiedSpace = centerInOccupiedSpace || !sameSide;
        }
      }

      // check if point is inside any sphere already inserted -> max overlap < radius
      if (!centerInOccupiedSpace) {
        bool acceptable = true;
        double dist_min_sq = radii[0];
        for (unsigned int idx2 = 0; idx2 < idx; ++idx2) {
          double dir[3];
          centerInOccupiedSpace = isPointInSphere(x[idx2], radii[idx2], x_rand, dir);

          if (centerInOccupiedSpace) {
            double dist_sq = dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2];
            double dist_acceptable = radii[idx2] - radii[idx];
            acceptable = (dist_sq > dist_acceptable*dist_acceptable);
            if (acceptable) {
              dist_min_sq = std::min(dist_min_sq, dist_sq);
            } else {
              break;
            }
          }
        }
        if (acceptable && dist_min_sq > most_acceptable_x[3]) {
          most_acceptable_x[0] = x_rand[0];
          most_acceptable_x[1] = x_rand[1];
          most_acceptable_x[2] = x_rand[2];
          most_acceptable_x[3] = dist_min_sq;
        }
      }

    } while (centerInOccupiedSpace && ntry < maxattempt);

    if (centerInOccupiedSpace && ntry == maxattempt && most_acceptable_x[3] > 0.0) {
      most_acceptable_x.resize(3);
      x[idx] = most_acceptable_x;
    } else {
      x[idx] = x_rand;
    }
  }
}


void ParticleSpatialDistribution::relax(
        double radius,
        const std::vector<double> &radii,
        std::vector<std::vector<double> > &x,
        const std::vector<double> &ext_radii,
        const std::vector<std::vector<double> > &ext_center)
{
  const unsigned int nParticles = radii.size();

  if (nParticles == 0)
    return;

  std::vector<double> fx(nParticles, 0.0);
  std::vector<double> fy(nParticles, 0.0);
  std::vector<double> fz(nParticles, 0.0);
  for (unsigned int i = 0; i < 500; ++i) {
    calculateForce(radii, x, fx, fy, fz);
    applyForce(radius, radii, x, ext_radii, ext_center, fx, fy, fz);
  }
}


void ParticleSpatialDistribution::calculateForce(
        const std::vector<double> &radii,
        const std::vector<std::vector<double> > &x,
        std::vector<double> &fx,
        std::vector<double> &fy,
        std::vector<double> &fz)
{
  const unsigned int nParticles = radii.size();

  if (nParticles == 0)
    return;

  fx.reserve(nParticles);
  fy.reserve(nParticles);
  fz.reserve(nParticles);
  std::fill(fx.begin(), fx.end(), 0.0);
  std::fill(fy.begin(), fy.end(), 0.0);
  std::fill(fz.begin(), fz.end(), 0.0);

  for(unsigned int i = 0; i < nParticles-1; ++i) {
    for(unsigned int j = i+1; j < nParticles; ++j) {
      double delta[3];
      delta[0] = x[i][0] - x[j][0];
      delta[1] = x[i][1] - x[j][1];
      delta[2] = x[i][2] - x[j][2];

      double dist_sq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];

      double radsum = radii[i] + radii[j];

      if (dist_sq < radsum*radsum) {
        const double dist = sqrt(dist_sq);
        const double inverse_dist = 1.0 / dist;
        const double overlap = radsum - dist;
        const double scale = inverse_dist * overlap;
        delta[0] *= scale;
        delta[1] *= scale;
        delta[2] *= scale;
        fx[i] += delta[0];
        fy[i] += delta[1];
        fz[i] += delta[2];
        fx[j] -= delta[0];
        fy[j] -= delta[1];
        fz[j] -= delta[2];
      }
    }
  }
}


void ParticleSpatialDistribution::applyForce(
        double radius,
        const std::vector<double> &radii,
        std::vector<std::vector<double> > &x,
        const std::vector<double> &ext_radii,
        const std::vector<std::vector<double> > &ext_center,
        const std::vector<double> &fx,
        const std::vector<double> &fy,
        const std::vector<double> &fz)
{
  const unsigned int nParticles = radii.size();
  const unsigned int ext = ext_radii.size();

  if (nParticles == 0)
    return;

  const double alpha = 0.1;
  const std::vector<double> center(3, 0.0);

  for (unsigned int i = 0; i < nParticles; ++i) {
    x[i][0] += fx[i] * alpha;
    x[i][1] += fy[i] * alpha;
    x[i][2] += fz[i] * alpha;

    double dir[3];
    for (unsigned int j = 0; j < ext; ++j) {
      std::vector<double> ext_c(3, 0.0);
      ext_c[0] = ext_center[j][0];
      ext_c[1] = ext_center[j][1];
      ext_c[2] = ext_center[j][2];

      if (isPointInSphere(ext_c, ext_radii[j] + radii[i], x[i], dir)) {
        // moved into contacting sphere
        // dir points from pos to center
        double dist = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        double inverse_dist = 1.0/dist;
        double scale = -inverse_dist * (ext_radii[j] - dist + radii[i]);

        dir[0] *= scale;
        dir[1] *= scale;
        dir[2] *= scale;

        x[i][0] += dir[0];
        x[i][1] += dir[1];
        x[i][2] += dir[2];
      }
    }

    if (!isPointInSphere(center, radius - radii[i], x[i], dir)) {
      // moved outside bounding sphere
      double dist = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
      double inverse_dist = 1.0/dist;
      double scale = inverse_dist * (dist - (radius - radii[i]));

      dir[0] *= scale;
      dir[1] *= scale;
      dir[2] *= scale;

      x[i][0] += dir[0];
      x[i][1] += dir[1];
      x[i][2] += dir[2];
    }
  }
}


double ParticleSpatialDistribution::overlap(
        const std::vector<double>& radii,
        const std::vector<std::vector<double> > &x
        ) const
{
  const unsigned int nParticles = radii.size();

  if (nParticles == 0)
    return 0.0;

  double totalOverlap = 0.0;

  for (unsigned int i = 0; i < nParticles-1; ++i) {
    for (unsigned int j = i+1; j < nParticles; ++j) {
      double delta[3];
      delta[0] = x[i][0] - x[j][0];
      delta[1] = x[i][1] - x[j][1];
      delta[2] = x[i][2] - x[j][2];

      double dist_sq = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];

      double radsum = radii[i] + radii[j];

      if (dist_sq < radsum*radsum)
          totalOverlap += radsum - sqrt(dist_sq);
    }
  }

  return totalOverlap;
}

