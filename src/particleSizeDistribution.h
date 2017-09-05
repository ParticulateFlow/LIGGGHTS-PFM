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

#ifndef LMP_PARTICLE_SIZE_DISTRIBUTION_H
#define LMP_PARTICLE_SIZE_DISTRIBUTION_H

#include <map>
#include <vector>

namespace LAMMPS_NS {

class ParticleSizeDistribution
{
 public:
  ParticleSizeDistribution(double P, double density, double rad_parent, double rad_min, double rad_max, double t10_max, double rad_omit=0.0, bool omit_post=false);

  void range_mass_fractions(std::map<double, double>& radiiRangeMassFractions);
  double radii(const std::map<double, double>& radiiRangeMassFractions, std::vector<double> &radii);

 private:
  void mass_fractions(std::map<double, double>& radiiMassFractions);
  double t10();
  double tn(double n);

 private:
  double breakage_probability;
  double density_;
  double rad_parent_;
  double rad_min_;
  double rad_max_;
  double t10_max_;
  double rad_omit_;
  double t10_;
  double mass_fraction_omit_;
  bool omit_post_;
};

}

#endif // LMP_PARTICLE_SIZE_DISTRIBUTION_H
