/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2017-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
   Daniel Queteschiner (JKU Linz)
   ---------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(insert/pack/dense,FixInsertPackDense)

#else

#ifndef LMP_FIX_INSERT_PACK_DENSE_H
#define LMP_FIX_INSERT_PACK_DENSE_H

#include "fix.h"

#include <list>

#include "region_neighbor_list.h"
#include "region_distance_field.h"

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
  typedef std::list<LIGGGHTS::Particle> ParticleList;
  typedef std::list<class ParticleToInsert*> PTIList;

  ParticleList frontSpheres;
  PTIList rejectedSpheres;

  ParticleList candidatePoints;

  class RanPark *random;
  int seed;

  LIGGGHTS::RegionNeighborList neighlist;
  LIGGGHTS::RegionDistanceField distfield;

  BoundingBox ins_bbox;

  bool insertion_done;
  double region_volume, region_volume_local;
  int n_insert_estim, n_insert_estim_local;
  double radius_factor;

  void prepare_insertion();

  void insert_first_particles();
  bool insert_next_particle(); // returns false if no insertion possible

  void handle_next_front_sphere();

  void compute_and_append_candidate_points(LIGGGHTS::Particle const &p1,
                                           LIGGGHTS::Particle const &p2,
                                           LIGGGHTS::Particle const &p3,
                                           double const r_insert,
                                           double const type_insert);
  void generate_initial_config(class ParticleToInsert *&p1,
                               class ParticleToInsert *&p2,
                               class ParticleToInsert *&p3);

  class ParticleToInsert* get_next_pti();

  LIGGGHTS::Particle particle_from_pti(class ParticleToInsert *pti);
  bool is_completely_in_subregion(LIGGGHTS::Particle &p);
  bool candidate_point_is_valid(LIGGGHTS::Particle &p);

  int n_inserted, n_inserted_local;
};

}

#endif /* FIX_INSERT_PACK_DENSE */
#endif /* FIX_CLASS */
