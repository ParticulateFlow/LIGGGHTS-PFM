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
#include <stdio.h>

using namespace LAMMPS_NS;
using namespace MathConst;

ParticleSizeDistribution::ParticleSizeDistribution(double P, double density, double rad_parent, double rad_min, double rad_max, double t10_max, double rad_omit, bool omit_post) :
  breakage_probability(P),
  density_(density),
  rad_parent_(rad_parent),
  rad_min_(rad_min),
  rad_max_(rad_max),
  t10_max_(t10_max),
  rad_omit_(rad_omit),
  t10_(0.0),
  mass_fraction_omit_(0.0),
  omit_post_(omit_post)
{
  if (rad_omit_ > rad_min_) {
    rad_min_ = rad_omit_;
  }
  if (rad_max_ < rad_min_) {
    fprintf(stdout,"WARNING: max radius < min radius: resetting max radius (%s:%d)\n",__FILE__,__LINE__);
    rad_max_ = rad_min_;
  }
}


// Note: std::map elements are sorted by their keys.
void ParticleSizeDistribution::mass_fractions(std::map<double, double>& radiiMassFractions)
{
  t10_ = t10();

  std::map<double,double>::iterator it = radiiMassFractions.begin();

  for (; it != radiiMassFractions.end(); ++it) {
    if (!omit_post_ && rad_parent_/it->first <= rad_omit_)
      break;
    it->second = tn(it->first);
  }

  radiiMassFractions.erase(it, radiiMassFractions.end());

  if (rad_omit_ >= rad_parent_) {
    mass_fraction_omit_ = 1.0;
  } else if (rad_omit_ > 0.0) {
    if (omit_post_) {
      mass_fraction_omit_ = 0.0;
    } else {
      mass_fraction_omit_ = tn(rad_parent_/rad_omit_);
    }
  }
}


void ParticleSizeDistribution::range_mass_fractions(std::map<double, double>& radiiRangeMassFractions)
{
  mass_fractions(radiiRangeMassFractions);

  if (!radiiRangeMassFractions.empty()) {
    std::map<double,double>::iterator it1 = radiiRangeMassFractions.begin();
    std::map<double,double>::iterator it2 = radiiRangeMassFractions.begin();

    for (++it2; it2 != radiiRangeMassFractions.end(); ++it1, ++it2) {
      it1->second = it1->second - it2->second;
    }

    if (rad_omit_ > 0.0) {
      it1->second = it1->second - mass_fraction_omit_;
    }
  }
}


double ParticleSizeDistribution::radii(const std::map<double, double>& radiiRangeMassFractions, std::vector<double>& radii)
{
  const double volume_parent = MY_4PI3 * rad_parent_*rad_parent_*rad_parent_;
  const double mass_parent = volume_parent * density_;
  double totalMassPool = mass_parent - mass_parent * mass_fraction_omit_;
  double mass_omitted = 0.0;

  std::map<double,double>::const_iterator it1 = radiiRangeMassFractions.begin();
  std::map<double,double>::const_iterator it2 = radiiRangeMassFractions.begin();

  for (++it2; it1 != radiiRangeMassFractions.end(); ++it1, ++it2) {
    double currentParticleRadius = std::min(rad_parent_/it1->first, rad_max_);
    double currentMassPool = mass_parent * it1->second;
    double nextParticleRadiusMax = 0.0;
    double nextParticleMassMax = 0.0;

    if (it2 != radiiRangeMassFractions.end()) {
      nextParticleRadiusMax = rad_parent_/it2->first;
      nextParticleMassMax = mass_parent/(it2->first*it2->first*it2->first);
    }

    if (rad_max_ < nextParticleRadiusMax) {
      fprintf(stdout,"WARNING: skipping radius fraction %f (%s:%d)\n",it1->first,__FILE__,__LINE__);
    }

    // as long as we have enough mass left to form a particle with radius larger than the next section ...
    while (currentMassPool       > nextParticleMassMax   &&
           currentParticleRadius > nextParticleRadiusMax &&
           currentParticleRadius > rad_min_) {
      double currentParticleMass = MY_4PI3 * currentParticleRadius*currentParticleRadius*currentParticleRadius * density_;
      while (currentMassPool >= currentParticleMass) {
        radii.push_back(currentParticleRadius);
        currentMassPool -= currentParticleMass;
        totalMassPool -= currentParticleMass;
      }
      currentParticleRadius *= 0.999;
    }

    if (omit_post_ && currentParticleRadius <= rad_omit_) {
      totalMassPool -= currentMassPool;
      mass_omitted += currentMassPool;
    }
  }

  while (totalMassPool > 0.0) {
    double radius = cbrt(totalMassPool/(MY_4PI3*density_));

    if (radius > rad_max_) {
      radius = rad_max_;
    } else if (radius < rad_min_) {
      radius += radii.back();
      radii.pop_back();
    }

    std::vector<double>::iterator it = radii.begin();
    while (it != radii.end() && *it > radius)
      ++it;
    radii.insert(it, radius); // slow op; use list?
    totalMassPool -= MY_4PI3*density_*radius*radius*radius;
  }

  if (omit_post_) {
    return mass_omitted;
  } else {
    return mass_parent * mass_fraction_omit_;
  }
}


double ParticleSizeDistribution::t10()
{
  return t10_max_ * breakage_probability;
}


double ParticleSizeDistribution::tn(double n)
{
  const double alpha = 1.0;
  return 1.0 - pow(1.0-t10_, alpha*9.0/(n-1.0));
}
