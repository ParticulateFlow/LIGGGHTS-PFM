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

#ifndef LMP_MULTISPHERE_PARALLEL_H
#define LMP_MULTISPHERE_PARALLEL_H

#include "multisphere.h"
#include "comm.h"

namespace LAMMPS_NS {

  class MultisphereParallel : public Multisphere {

    friend class FixMultisphere;

    public:

      MultisphereParallel(LAMMPS *lmp);
      ~MultisphereParallel();

      int pack_exchange_rigid(int i, double *buf);
      int unpack_exchange_rigid(double *buf);

      void writeRestart(FILE *fp);
      void restart(double *list);

    private:

      //NP buffers for body communication
      void exchange();
      void grow_send(int, int);
      void grow_recv(int);

      // current size of send/recv buffer
      // send buffer and recv buffer for all comm
      int maxsend_,maxrecv_;
      double *buf_send_;
      double *buf_recv_;
  };

  // *************************************
  #include "multisphere_parallel_I.h"
  // *************************************
}

#endif
