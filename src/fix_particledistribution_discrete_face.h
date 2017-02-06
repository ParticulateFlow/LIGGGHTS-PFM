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


#ifdef FIX_CLASS

FixStyle(particledistribution/discrete/face,FixParticledistributionDiscreteFace)

#else

#ifndef LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_FACE_H
#define LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_FACE_H

#include <vector>
#include "fix_particledistribution.h"
#include "constParticleTemplateSphere.h"


namespace LAMMPS_NS {

class FixParticledistributionDiscreteFace : public FixParticledistribution {
 public:
  friend class FixPourDev;
  FixParticledistributionDiscreteFace(class LAMMPS *, int, char **);
  ~FixParticledistributionDiscreteFace();

  double min_rad(int);
  double max_rad(int);

  double min_rad()
  { return minrad; }
  double max_rad()
  { return maxrad; }
  double max_r_bound()
  { return maxrbound; }

  void random_init_list(int);
  int randomize_list(int,int,int);     // generate a list of random particles

  std::vector<std::vector<ParticleToInsert*> > pti_list_face_local;

  void pre_insert(int n=0, FixPropertyAtom *fp=NULL, double val=0., int idx=-1, int ival=0, int iidx=-1);

  int insert(int n);

  void set_distribution_local(const std::vector<DiscreteParticleDistribution>& distributions, const std::vector<std::vector<int> > & distributions_face_local, double cg, int type_offset);

 protected:

  void delete_pit_list_face_local();
};

}

#endif
#endif
