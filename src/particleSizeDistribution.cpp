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

#include <cmath>
#include "math_const.h"
#include "particleSizeDistribution.h"

using namespace LAMMPS_NS;
using namespace MathConst;

ParticleSizeDistribution::ParticleSizeDistribution(double P, double density, double rad_parent, double rad_min, double t10_max) :
  breakage_probability(P),
  density_(density),
  rad_parent_(rad_parent),
  rad_min_(rad_min),
  t10_max_(t10_max),
  t10_(0.0)
{
}


// Note: std::map elements are sorted by their keys.
void ParticleSizeDistribution::mass_fractions(std::map<int, double>& radiiMassFractions)
{
  t10_ = t10();

  std::map<int,double>::iterator it = radiiMassFractions.begin();

  for (; it != radiiMassFractions.end(); ++it) {
    it->second = tn(it->first);
  }
}


void ParticleSizeDistribution::range_mass_fractions(std::map<int, double>& radiiRangeMassFractions)
{
  mass_fractions(radiiRangeMassFractions);

  std::map<int,double>::iterator it1 = radiiRangeMassFractions.begin();
  std::map<int,double>::iterator it2 = radiiRangeMassFractions.begin();

  for (++it2; it2 != radiiRangeMassFractions.end(); ++it1, ++it2) {
    it1->second = it1->second - it2->second;
  }
}


void ParticleSizeDistribution::radii(const std::map<int, double>& radiiRangeMassFractions, std::vector<double>& radii)
{
  const double volume_parent = MY_4PI3 * rad_parent_*rad_parent_*rad_parent_;
  const double mass_parent = volume_parent * density_;
  double totalMassPool = mass_parent;

  std::map<int,double>::const_iterator it1 = radiiRangeMassFractions.begin();
  std::map<int,double>::const_iterator it2 = radiiRangeMassFractions.begin();

  for (++it2; it1 != radiiRangeMassFractions.end(); ++it1, ++it2) {
    double currentParticleRadius = rad_parent_/it1->first;
    double currentMassPool = mass_parent * it1->second;
    double nextParticleRadiusMax = 0.0;
    double nextParticleMassMax = 0.0;

    if (it2 != radiiRangeMassFractions.end()) {
      nextParticleRadiusMax = rad_parent_/it2->first;
      nextParticleMassMax = mass_parent/(it2->first*it2->first*it2->first);
    }

    // as long as we have enough mass left to form a particle with radius larger than the next section ...
    while (currentMassPool       >  nextParticleMassMax   &&
           currentParticleRadius >  nextParticleRadiusMax &&
           currentParticleRadius >= rad_min_) {
      double currentParticleMass = MY_4PI3 * currentParticleRadius*currentParticleRadius*currentParticleRadius * density_;
      while (currentMassPool >= currentParticleMass) {
        radii.push_back(currentParticleRadius);
        currentMassPool -= currentParticleMass;
        totalMassPool -= currentParticleMass;
      }
      currentParticleRadius *= 0.999;
    }
  }

  if (totalMassPool > 0.0) {
    double radius = cbrt(totalMassPool/(MY_4PI3*density_));

    if (radius < rad_min_) {
      radius += radii.back();
      radii.pop_back();
    }

    std::vector<double>::iterator it = radii.begin();
    while (it != radii.end() && *it > radius)
      ++it;
    radii.insert(it, radius); // slow op; use list?
  }
}


double ParticleSizeDistribution::t10()
{
  return t10_max_ * breakage_probability;
}


double ParticleSizeDistribution::tn(int n)
{
  double alpha = 1.0;
  return 1.0 - pow(1.0-t10_, alpha*9.0/(n-1.0));
}
