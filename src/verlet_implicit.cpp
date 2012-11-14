/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "verlet_implicit.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
/*NL*/#include "debug_liggghts.h"

using namespace LAMMPS_NS;

/*NL*/ #define DEBUG_VERLET_IMP false //(update->ntimestep>54500) //false  true //

/* ---------------------------------------------------------------------- */

VerletImplicit::VerletImplicit(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   run for N steps iteratively
------------------------------------------------------------------------- */

void VerletImplicit::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_post_force = modify->n_post_force;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"starting time-step " BIGINT_FORMAT "\n",update->ntimestep);__debug__(lmp);}

    do
    {
        ev_set(ntimestep);

        // initial time integration

        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing initial integrate\n");__debug__(lmp);}
        modify->initial_integrate(vflag);
        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing post integrate\n");__debug__(lmp);}
        if (n_post_integrate) modify->post_integrate();

        // regular communication vs neighbor list rebuild

        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing neigh stuff\n");__debug__(lmp);}
        nflag = neighbor->decide();

        if (nflag == 0) {
          timer->stamp();
          comm->forward_comm();
          timer->stamp(TIME_COMM);
        } else {
          /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing pre_exchange\n");__debug__(lmp);}
          if (n_pre_exchange) modify->pre_exchange();
          /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing pbc, reset_box, comm setup\n");__debug__(lmp);}
          if (triclinic) domain->x2lamda(atom->nlocal);
          domain->pbc();
          if (domain->box_change) {
            domain->reset_box();
            /*NL*/ //if(balanceflag) domain->loadbalance_local_boxes();
            comm->setup();
            if (neighbor->style) neighbor->setup_bins();
          }
          timer->stamp();
          /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing exchange\n");__debug__(lmp);}
          comm->exchange();
          /*NL*/ //fprintf(screen,"proc %d has %d owned particles\n",comm->me,atom->nlocal);
          if (sortflag && ntimestep >= atom->nextsort) atom->sort();
          comm->borders();
          if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
          timer->stamp(TIME_COMM);
          /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing pre neigh\n");__debug__(lmp);}
          if (n_pre_neighbor) modify->pre_neighbor();
          /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing neigh\n");__debug__(lmp);}
          neighbor->build();
          timer->stamp(TIME_NEIGHBOR);
        }

        // force computations
        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing pre force\n");__debug__(lmp);}
        force_clear();
        if (n_pre_force) modify->pre_force(vflag);

        timer->stamp();

        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing pair force\n");__debug__(lmp);}
        if (force->pair) {
          force->pair->compute(eflag,vflag);
          timer->stamp(TIME_PAIR);
        }

        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing molecular, kspace\n");__debug__(lmp);}
        if (atom->molecular) {
          if (force->bond) force->bond->compute(eflag,vflag);
          if (force->angle) force->angle->compute(eflag,vflag);
          if (force->dihedral) force->dihedral->compute(eflag,vflag);
          if (force->improper) force->improper->compute(eflag,vflag);
          timer->stamp(TIME_BOND);
        }

        if (kspace_compute_flag) {
          force->kspace->compute(eflag,vflag);
          timer->stamp(TIME_KSPACE);
        }

        // reverse communication of forces
        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing reverse comm\n");__debug__(lmp);}
        if (force->newton) {
          comm->reverse_comm();
          timer->stamp(TIME_COMM);
        }

        // force modifications, final time integration, diagnostics
        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing post force\n");__debug__(lmp);}
        if (n_post_force) modify->post_force(vflag);
        /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing final integrate\n");__debug__(lmp);}
        modify->final_integrate();

    }
    while(modify->iterate_implicitly());


    /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing end of step\n");__debug__(lmp);}
    if (n_end_of_step) modify->end_of_step();

    // all output
    /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    doing output\n");__debug__(lmp);}
    if (ntimestep == output->next) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(TIME_OUTPUT);
    }
    /*NL*/if(DEBUG_VERLET_IMP) {MPI_Barrier(world);if(comm->me==0)fprintf(screen,"    ts finished\n");__debug__(lmp);}
  }
}
