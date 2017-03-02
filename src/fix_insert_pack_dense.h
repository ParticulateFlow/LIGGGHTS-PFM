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

/* ----------------------------------------------------------------------
   contributing authors
   Philippe Seil (JKU)
   ---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/pack/dense,FixInsertPackDense)

#else

#ifndef LMP_FIX_INSERT_PACK_DENSE_H
#define LMP_FIX_INSERT_PACK_DENSE_H

#include "fix.h"

#include <list>

#include "region_neighbor_list.h"

using namespace LIGGGHTS;

namespace LAMMPS_NS {

class FixInsertPackDense : public Fix {
 public:

  FixInsertPackDense(class LAMMPS *, int, char **);
  ~FixInsertPackDense();

  // virtual void restart(char *);

  virtual int setmask();

  virtual void post_create();
  
  virtual void pre_exchange();
protected:

  double *x_init; // starting point for insertion

  double maxrad; // maximum radius to be expected
  
  // region to be used for insertion
  class Region *ins_region;
  char *idregion;

  //particle distribution
  class FixParticledistributionDiscrete *fix_distribution;

  static const double max_volfrac;
  double target_volfrac;
  typedef std::list<Particle> ParticleList;
  typedef std::list<class ParticleToInsert*> PTIList;

  ParticleList frontSpheres;
  PTIList rejectedSpheres;

  ParticleList candidatePoints;
  
  class RanPark *random;
  int seed;

  RegionNeighborList neighlist;

  BoundingBox ins_bbox;

  bool insertion_done;
  double region_volume, region_volume_local;
  int n_insert_estim, n_insert_estim_local;
  double radius_factor;
  
  void prepare_insertion();
  
  void insert_first_particles();
  bool insert_next_particle(); // returns false if no insertion possible

  void handle_next_front_sphere();
  
  void compute_and_append_candidate_points(Particle const &p1,
                                           Particle const &p2,
                                           Particle const &p3,
                                           double const r_insert);
  void generate_initial_config(class ParticleToInsert *&p1,
                               class ParticleToInsert *&p2,
                               class ParticleToInsert *&p3);

  class ParticleToInsert* get_next_pti();
  
  Particle particle_from_pti(class ParticleToInsert *pti);
  bool is_completely_in_subregion(Particle &p);
  bool candidate_point_is_valid(Particle &p);

  int n_inserted, n_inserted_local;
};

}

#endif /* FIX_INSERT_PACK_DENSE */
#endif /* FIX_CLASS */
