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

#include <queue>

#include "region_neighbor_list.h"
#include "region_distance_field.h"

namespace LAMMPS_NS {

class FixPropertyAtom;

class FixInsertPackDense : public Fix {
 public:

  FixInsertPackDense(class LAMMPS *, int, char **);
  ~FixInsertPackDense();

  // virtual void restart(char *);

  virtual int setmask();

  virtual void post_create();

  virtual void init();

  virtual void pre_exchange();
protected:

  double *x_init; // starting point for insertion

  double maxrad; // maximum radius to be expected

  // region to be used for insertion
  class Region *ins_region;
  char *idregion;

  // variable used to identify if insertion should take place
  char *var;
  double var_insvalid;
  double var_threshold;

  //particle distribution
  class FixParticledistributionDiscrete *fix_distribution;

  static const double max_volfrac;
  double target_volfrac;

  typedef std::vector<LIGGGHTS::Particle> ParticleVector;
  typedef std::queue<LIGGGHTS::Particle> ParticleQueue;
  typedef std::queue<class ParticleToInsert*> PTIList;

  ParticleQueue frontSpheres;
  PTIList rejectedSpheres;

  ParticleVector candidatePoints;

  class RanPark *random;
  int seed;

  LIGGGHTS::RegionNeighborList neighlist;
  LIGGGHTS::RegionDistanceField distfield;

  BoundingBox ins_bbox;

  bool insertion_done;
  bool is_inserter; // indicates if proc inserts
  double region_volume, region_volume_local;
  int n_insert_estim, n_insert_estim_local;
  int insert_every, most_recent_ins_step;
  double radius_factor;

  bool prepare_insertion();

  bool insert_first_particles(); // returns false if no insertion possible

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

private:
  char *property_name;
  FixPropertyAtom *fix_property;
  double fix_property_value;
  int fix_property_ivalue;
  int property_index;
  int property_iindex;
};

}

#endif /* FIX_INSERT_PACK_DENSE */
#endif /* FIX_CLASS */
