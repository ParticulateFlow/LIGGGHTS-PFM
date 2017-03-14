
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

#include "fix.h"

enum{RAN_STYLE_CONSTANT_FPD,RAN_STYLE_UNIFORM_FPD,RAN_STYLE_GAUSSIAN_FPD};


namespace LAMMPS_NS {

class FixParticledistributionDiscrete : public Fix {
 public:
  friend class FixPourDev;
  FixParticledistributionDiscrete(class LAMMPS *, int, char **);
  ~FixParticledistributionDiscrete();
  void write_restart(FILE *);
  void restart(char *);

  int setmask();

  double vol_expect();
  double mass_expect();

  int max_type();
  int min_type();

  double min_rad(int);
  double max_rad(int);

  double min_rad() const
  { return minrad; }
  double max_rad() const
  { return maxrad; }
  double max_r_bound() const
  { return maxrbound; }

  int max_nspheres();

  int random_init_single(int);         //NP tells the distribution how many particles are to be generated
  class Region* randomize_single();    //NP generates one random particle

  void random_init_list(int);
  int randomize_list(int,int,int);     //NP generates a list of random particles

  void finalize_insertion();

  class ParticleToInsert *pti;

  class ParticleToInsert **pti_list;
  int n_pti, n_pti_max;
  void pre_insert();
  int insert(int n);

  inline int n_particletemplates()
  { return ntemplates; }
  inline class FixTemplateSphere** particletemplates()
  { return templates; }

 protected:

  //NP called at insertion
  int ninserted;
  int ninsert;

  class RanPark *random;
  int seed;

  int iarg;

  // particle templates
  int ntemplates;
  double *distweight;
  double *cumweight;
  int *parttogen;
  int *distorder;
  class FixTemplateSphere **templates;

  // mass and volume expectancy of this discrete distribution
  double volexpect;
  double massexpect;

  // min/maximum particle type to be inserted
  int maxtype;
  int mintype;

  // maximum number of spheres a template has
  int maxnspheres;

  // maximum radius and bounding sphere radius
  double minrad,maxrad,maxrbound;
};

}

#endif
#endif
