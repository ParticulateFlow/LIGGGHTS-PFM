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

#ifndef LMP_DEBUG_LIGGGHTS_H
#define LMP_DEBUG_LIGGGHTS_H

#include "lammps.h"
#include "comm.h"
#include "string.h"
#include "stdlib.h"
#include "style_fix.h"
#include "vector_liggghts.h"
#include "container.h"
#include "atom.h"

namespace LAMMPS_NS {


inline void __trace__()
{
    /*NL*/ //double test[3];
    /*NL*/ //test[4] = 0.;
}

inline void __debug__(LAMMPS* lmp)
{
    /*NL*///fprintf(lmp->screen,"bond_hist %d bond_per_atom %d\n",lmp->atom->bond_hist,lmp->atom->bond_per_atom);
    /*NL*///fprintf(lmp->screen,"step " BIGINT_FORMAT " nparticles %d \n",lmp->update->ntimestep,lmp->atom->nlocal);
    /*NL*///printVec3D(lmp->screen,"pos for tag 206",lmp->atom->x[lmp->atom->map(206)]);
    /*NL*///printVec3D(lmp->screen,"vel for tag 206",lmp->atom->v[lmp->atom->map(206)]);
    /*NL*/if(1==lmp->comm->me)printVec3D(lmp->screen,"f for atom 253",lmp->atom->f[253]);//lmp->atom->f[lmp->atom->map(1)]);
     /*NL*///Atom *atom = lmp->atom;
     /*NL*///int nlocal = atom->nlocal;
     /*NL*///fprintf(lmp->screen,"nlocal %d\n",nlocal);
     /*NL*///int nghost = atom->nghost;
     /*NL*///int nall = atom->nlocal+atom->nghost;

     /*NL*///for(int i = 0; i < nlocal; i++)
     /*NL*///{
     /*NL*///  printVec3D(lmp->screen,"pos",atom->x[i]);
     /*NL*/  //if(atom->tag[i] == 25 || atom->tag[i] == 30 )fprintf(lmp->screen,"tag %d vel %f %f %f \n",atom->tag[i],atom->v[i][0],atom->v[i][1],atom->v[i][2]);
     /*NL*/  //if(atom->map(atom->tag[i]) >= nlocal)
     /*NL*/   //lmp->error->all(FLERR,"catch");
     /*NL*///}

/*NP
    for(int i = 0; i < lmp->modify->nfix; i++)
    {

        if(strcmp(lmp->modify->fix[i]->style,"move/mesh") == 0)
        {
            FixMoveMesh *fmm = static_cast<FixMoveMesh*>(lmp->modify->fix[i]);
            MultiVectorContainer<double,3,3> *v;
            v = fmm->mesh()->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
            if(v) printVec3D(lmp->screen,"VDEBUG",v->begin()[0][0]);

        }
    }
*/
}

}

#endif
