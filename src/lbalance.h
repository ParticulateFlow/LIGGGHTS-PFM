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

//NP pure virtual class

#ifndef LMP_LBALANCE_H
#define LMP_LBALANCE_H

#include "neighbor.h"
#include "string.h"

#define BIG   1.0e20
#define SMALL 1.0e-4
#define DELTA 1
#define FABS(a) ((a) > 0 ? (a) : (-a))

namespace LAMMPS_NS {

class Lbalance : protected Pointers {

   public:
    Lbalance(class LAMMPS *lmp, int narg, char **arg): Pointers(lmp)
    {
        int n = strlen(arg[0]) + 1;
        style = new char[n];
        strcpy(style,arg[0]);
        iarg = 1;
    }

    virtual ~Lbalance()
    {
        delete []style;
    }

    virtual void loadbalance_local_boxes()=0;

    double cutneighmax() {return neighbor->cutneighmax;}

   protected:
    int iarg;
    char *style;
};
}

#endif
