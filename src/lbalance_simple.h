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

/*NP
simple load balancing algorithm

makes three sets independant axis-parallel cuts
where an optimal number of particles is assigned to each slice

uses binary search for this

exact for low number of processes, e.g. 8

does not guarantee optimal load balance for more processes
*/


#ifdef LB_CLASS

LBStyle(nlocal/simple,LbalanceSimple)

#else


#ifndef LMP_LBALANCE_SIMPLE_H
#define LMP_LBALANCE_SIMPLE_H

#include "lbalance.h"

namespace LAMMPS_NS {

class LbalanceSimple : public Lbalance {

   public:
    LbalanceSimple(class LAMMPS *lmp, int narg, char **arg);
    ~LbalanceSimple();

   protected:
    virtual void loadbalance_local_boxes();

    void loadbalance_local_boxes_simple(int *lodim, int *hidim);

    void calc_max_shift(int ,int);
    int count_particles(int,double);
    double calc_border(int,int,int,double *,int *);
    void apply_border(double *bal_sublo, double *bal_subhi,int idim,int idim_proc);

    int ntry_simple;
    int lodim[3],hidim[3]; //NP proc stencil

    int *myloc;
    int *procgrid;
    double boxhi_stencil[3];
    double *boxhi,*boxlo;
    double *subhi,*sublo;

    double minextent;
    double max_shift[2]; //0 is for lowering the border, 1 for raising it
};
}

#endif
#endif
