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


// abstract class
#ifdef FIX_CLASS

#else

#ifndef LMP_FIX_PARTICLEDISTRIBUTION_H
#define LMP_FIX_PARTICLEDISTRIBUTION_H

#include <vector>
#include "fix.h"

enum{RAN_STYLE_CONSTANT_FPD,RAN_STYLE_UNIFORM_FPD,RAN_STYLE_GAUSSIAN_FPD};


namespace LAMMPS_NS {

class ParticleToInsert;
class FixPropertyAtom;

class FixParticledistribution : public Fix {
 public:
  FixParticledistribution(class LAMMPS *, int, char **);
  ~FixParticledistribution();
  void write_restart(FILE *);
  virtual void restart(char *);

  virtual int setmask();

  double vol_expect() const
  { return volexpect; }
  double mass_expect() const
  { return massexpect; }

  int max_type() const
  { return maxtype; }
  int min_type() const
  { return mintype; }

  virtual double min_rad(int) const { return 0.0; }
  virtual double max_rad(int) const { return 0.0; }

  double min_rad() const
  { return minrad; }
  double max_rad() const
  { return maxrad; }
  double max_r_bound() const
  { return maxrbound; }

  int max_nspheres() const
  { return maxnspheres; }

  virtual void random_init_list(int) = 0;
  virtual int randomize_list(int,int,int) = 0;   // generate a list of random particles

  typedef std::vector<ParticleToInsert*> pti_list_type;
  pti_list_type pti_list;

  virtual void pre_insert(int n=0, FixPropertyAtom *fp=NULL, double val=0., int idx=-1, int ival=0, int iidx=-1);
  virtual int insert(int n);
  // pti_delete_flag allows to delete pti in pti_list
  // only set if pti were inserted from outside
  virtual void finalize_insertion(bool pti_delete_flag = false);

  virtual int n_particletemplates() const { return -1; }
  virtual class FixTemplateSphere** particletemplates() { return NULL; }

 protected:

  int ninserted;
  int ninsert;

  class RanPark *random;
  int seed;

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
