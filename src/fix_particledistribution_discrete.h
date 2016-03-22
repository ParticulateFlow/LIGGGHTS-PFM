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

FixStyle(particledistribution/discrete,FixParticledistributionDiscrete)

#else

#ifndef LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H
#define LMP_FIX_PARTICLEDISTRIBUTION_DISCRETE_H

#include "fix_particledistribution.h"


namespace LAMMPS_NS {

class FixParticledistributionDiscrete : public FixParticledistribution {
 public:
  friend class FixPourDev;
  FixParticledistributionDiscrete(class LAMMPS *, int, char **);
  ~FixParticledistributionDiscrete();

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

  void pre_insert(int n=0, FixPropertyAtom *fp=NULL, double val=0., int idx=-1, int ival=0);
  int insert(int n);
  void finalize_insertion();

  inline int n_particletemplates()
  { return ntemplates; }
  inline class FixTemplateSphere** particletemplates()
  { return templates; }

protected:

  // particle templates
  int ntemplates;
  double *distweight;
  double *cumweight;
  int *parttogen;
  int *distorder;
  class FixTemplateSphere **templates;

};

}

#endif
#endif
