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

FixStyle(insert/stream/moving,FixInsertStreamMoving)

#else

#ifndef LMP_FIX_INSERT_STREAM_MOVING_H
#define LMP_FIX_INSERT_STREAM_MOVING_H

#include "fix_insert_stream.h"

namespace LAMMPS_NS {

class FixInsertStreamMoving : public FixInsertStream {
 public:

  FixInsertStreamMoving(class LAMMPS *, int, char **);
  ~FixInsertStreamMoving();
  void post_create();

  virtual void init();
  virtual void end_of_step();

  virtual int release_step_index()
  { return 6; }

 private:

  virtual void finalize_insertion(int);

  class TriMeshPlanar *ins_face_planar;
};

}

#endif
#endif
