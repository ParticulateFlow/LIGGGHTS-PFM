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

FixStyle(remove,FixRemove)

#else

#ifndef LMP_FIX_REMOVE_H
#define LMP_FIX_REMOVE_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixRemove : public Fix {

 friend class FixMultisphere;

 public:

  FixRemove(class LAMMPS *, int, char **);
  ~FixRemove();

  int setmask();
  void init();
  void pre_exchange();
  void write_restart(FILE *fp);
  void restart(char *buf);
  double compute_vector(int n);

 private:

  bool count_eligible(double &mass_eligible_me,double &mass_eligible,
              double &mass_shrink_me,double &mass_to_remove_me,
              double &ratio_ms_to_remove_me);
  void delete_all(double mass_eligible_me,double ratio_ms_to_remove_me,
              double &mass_removed_this_me,double &heat_removed_this_me,
              int &nremoved_this_me);
  void shrink(double &mass_to_remove_me,double mass_shrink_me,
              double &mass_removed_this_me,double &heat_removed_this_me,
              int &nremoved_this_me);
  void delete_partial_particles(double &mass_to_remove_me,
              double &mass_removed_this_me,double &heat_removed_this_me,
              int &nremoved_this_me);
  void delete_partial_particles_bodies(double &mass_to_remove_me,
              double &mass_removed_this_me,double &heat_removed_this_me,
              int &nremoved_this_me,double ratio_ms_to_remove_me);

  void delete_particle(int);
  void delete_bodies();

  class RanPark *random_;

  int type_remove_;                 //NP atom type of particles to remove
  int iregion_;
  int style_;
  double delete_below_;
  double rate_remove_;              //NP in kg per sec
  double m_remove_min_;             // minimal amount of mass to be removed at once
  int seed_;
  bool integrated_error_;

  bool variable_rate_;
  char* rate_name_;
  int ivar_;

  double mass_removed_,mass_to_remove_,integrated_rate_;
  int time_origin_,time_last_;
  double dt_;
  bool restart_read_;
  bool restart_write_;
  bool verbose_;
  int compress_flag_;

  bool monitor_heat_;
  double heat_removed_;
  class FixPropertyAtom *fix_temp_;
  class FixPropertyGlobal* fix_capacity_;

  class FixMultisphere *fix_ms_;
  class MultisphereParallel *ms_;

  //NP list of bodies to remove
  std::vector<int> atom_tags_eligible_;
  std::vector<int> body_tags_eligible_;
  std::vector<int> body_tags_delete_;
};

}

#endif
#endif
