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

#ifndef LMP_COMM_I_H
#define LMP_COMM_I_H

#include "atom.h"
#include "domain_wedge.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   decide if use comm optimizations for granular systems
   don't use for triclinic (current implementation not valid for triclinic,
   would have to translate radius into triclinic coordinates)
------------------------------------------------------------------------- */

inline bool Comm::use_gran_opt()
{
    return (0 == domain->triclinic && atom->radius);
}

/* ----------------------------------------------------------------------
   decide if border element, optimization for granular
------------------------------------------------------------------------- */

inline bool Comm::decide(int i,int dim,double lo,double hi,int ineed)
{
    double **x = atom->x;
    double *radius = atom->radius;

    /*NL*/ //if (screen) fprintf(screen,"called 'decide'\n");

    //NP only extracted half cutghost before in setup()
    //NP so account for radius here
    //NP going left first line / going right second line
    if( ((ineed % 2 == 0) && x[i][dim] >= lo && x[i][dim] <= (hi + (use_gran_opt()? (radius[i]) : 0.)) ) ||
        ((ineed % 2 == 1) && x[i][dim] >= (lo - (use_gran_opt()? radius[i] : 0.)) && x[i][dim] <= hi )   )
        return true;

    return false;
}

/* ----------------------------------------------------------------------
   decide if border element for wedge case, optimization for granular
------------------------------------------------------------------------- */

inline bool Comm::decide_wedge(int i,int dim,double lo,double hi,int ineed)
{
    double **x = atom->x;
    double *radius = atom->radius;
    double coo[2],d[2];
    coo[0] = x[i][iphi];
    coo[1] = x[i][(iphi+1)%3];

    /*NL*/ //if (screen) printVec3D(screen,"coo",x[i]);
    /*NL*/ //if (screen) fprintf(screen,"coo %f %f\n",coo[0],coo[1]);
    /*NL*/ //error->all(FLERR,"end");

    //NP going left
    if (ineed % 2 == 0)
    {
        vectorSubtract2D(coo,pleft,d);
        if(vectorDot2D(d,nleft) >= -(use_gran_opt()? radius[i] : 0.))
        {
            /*NL*/ //if (screen) fprintf(screen,"particle tag %d: left YES, vectorDot2D(d,nleft) %f ineed %d\n",atom->tag[i],vectorDot2D(d,nleft),ineed);
            /*NL*/ //if (screen) printVec3D(screen," x",atom->x[i]);
            /*NL*/ //if (screen) printVec2D(screen," coo",coo);
            /*NL*/ //if (screen) printVec2D(screen," d",d);
            /*NL*/ //if (screen) printVec2D(screen," pleft",pleft);
            return true;
        }
    }
    //NP going right
    else if (ineed % 2 == 1)
    {
        vectorSubtract2D(coo,pright,d);
        if(vectorDot2D(d,nright) >= -(use_gran_opt()? radius[i] : 0.))
        {
            /*NL*/ //if (screen) fprintf(screen,"particle tag %d: right YES, vectorDot2D(d,nright) %f ineed %d\n",atom->tag[i],vectorDot2D(d,nright),ineed);
            /*NL*/ //if (screen) printVec3D(screen," x",atom->x[i]);
            /*NL*/ //if (screen) printVec2D(screen," d",d);
            /*NL*/ //if (screen) printVec2D(screen," pright",pright);
            return true;
        }
    }
    /*NL*/ //if (screen) fprintf(screen,"particle tag %d: NO, vectorDot2D(d,nl/r) %f %f ineed %d\n",atom->tag[i],vectorDot2D(d,nleft),vectorDot2D(d,nright),ineed);
    /*NL*/ //if (screen) printVec3D(screen," x",atom->x[i]);
    /*NL*/ //if (screen) printVec2D(screen," coo",coo);
    /*NL*/ //if (screen) printVec2D(screen," d",d);
    /*NL*/ //if (screen) printVec2D(screen," pleft",pleft);
    /*NL*/ //if (screen) printVec2D(screen," pright",pright);
    return false;
}

#endif
