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

#ifndef LMP_PARTICLE_SPATIAL_DISTRIBUTION_H
#define LMP_PARTICLE_SPATIAL_DISTRIBUTION_H

namespace LAMMPS_NS {

class RanPark;
class PrimitiveWall;
class TriMeshContacts;

class ParticleSpatialDistribution
{
 public:
  ParticleSpatialDistribution(RanPark *rp, double overlap, int maxattempt);
  ~ParticleSpatialDistribution();

  bool isPointInSphere(const std::vector<double> &center, double radius, const std::vector<double> &x, double *dir=NULL, double *dist=NULL);
  void randomPointInSphere(double radius, std::vector<double> &x);

  void apollonianInsertion(double radius, const std::vector<double> &radii, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z);

  void randomInsertion(double *xp, double radius,
                       const std::vector<double> &radii, std::vector<std::vector<double> > &x,
                       const std::vector<std::vector<double> >& ext_atoms,
                       const std::vector<PrimitiveWall*> &prim_walls,
                       const std::vector<TriMeshContacts*> &meshes);

  void relax(double radius,
             const std::vector<double> &radii, std::vector<std::vector<double> > &x,
             const std::vector<double> &ext_radii, const std::vector<std::vector<double> > &ext_center);

  void calculateForce(const std::vector<double> &radii,
                      const std::vector<std::vector<double> > &x,
                      std::vector<double> &fx, std::vector<double> &fy, std::vector<double> &fz);

  void applyForce(double radius,
                  const std::vector<double> &radii, std::vector<std::vector<double> > &x,
                  const std::vector<double> &ext_radii, const std::vector<std::vector<double> > &ext_center,
                  const std::vector<double> &fx, const std::vector<double> &fy, const std::vector<double> &fz);

  double overlap(const std::vector<double> &radii, const std::vector<std::vector<double> > &x) const;

 private:
  RanPark * RNG;
  double max_overlap;
  int maxattempt;
};

}

#endif // SPATIALPARTICLEDISTRIBUTION_H
