/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   Copyright 2015-     JKU Linz

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

#ifndef LMP_CONST_PARTICLE_TEMPLATE_SPHERE_H
#define LMP_CONST_PARTICLE_TEMPLATE_SPHERE_H

#include <cmath>
#include <map>
#include "math_const.h"

namespace LAMMPS_NS {
  class ConstantParticleTemplateSphere
  {
   public:
    ConstantParticleTemplateSphere() :
      radius_(0.), density_(0.), mass_(0), atomtype_(0) {}
    ConstantParticleTemplateSphere(double radius, double mass, int atomtype) :
      radius_(radius), mass_(mass), atomtype_(atomtype) { density_ = mass/((4./3.)*MathConst::MY_PI*radius*radius*radius); }

    double radius_, density_, mass_;
    int atomtype_;

    bool operator<(const ConstantParticleTemplateSphere &other) const
    {
      //return radius < other.radius; // sort smallest to largest radius
      if(radius_ > other.radius_) return true; // sort largest to smallest radius
      else if(radius_ < other.radius_) return false;

      if(mass_ > other.mass_) return true; // sort largest to smallest mass
      else if(mass_ < other.mass_) return false;

      if(atomtype_ > other.atomtype_) return true;
      /*else if(atomtype_ < other.atomtype_) return false;*/

      return false; // always return false for equivalent templates
    }
  };

  typedef std::map<ConstantParticleTemplateSphere,double> DiscreteParticleDistribution; // template, #particles
}

#endif
