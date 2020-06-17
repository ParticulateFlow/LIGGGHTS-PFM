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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include <math.h>
#include "vector_liggghts.h"
#include "fix_cfd_coupling.h"
#include "fix_multisphere.h"
#include "cfd_datacoupling_one2one.h"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

CfdDatacouplingOne2One::CfdDatacouplingOne2One(LAMMPS *lmp,int iarg, int narg, char **arg,FixCfdCoupling* fc) :
  CfdDatacoupling(lmp, iarg, narg, arg,fc)
{
  liggghts_is_active = false;

  if(!atom->tag_enable) error->one(FLERR,"CFD-DEM coupling via MPI requires particles to have tags");

  this->fc_ = fc;

  if(comm->me == 0) error->message(FLERR,"nevery as specified in LIGGGHTS is overriden by calling external program",1);

  //NP do not make inital grow; this is done at the first ts together with callers arrays
  //NP this is to ensure that caller arrays and allred arrays always have same length
}

CfdDatacouplingOne2One::~CfdDatacouplingOne2One()
{}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingOne2One::exchange()
{
    // does nothing since done by OF
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingOne2One::pull(const char *name,const char *type,void *&from,const char *datatype)
{
    CfdDatacoupling::pull(name,type,from,datatype);
/*
    if(strcmp(datatype,"double") == 0)
        pull_mpi<double>(name, type, from, world);
    else if(strcmp(datatype,"int") == 0)
        pull_mpi<int>(name, type, from, world);
    else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, valid datatypes are 'int' and double'");

    */
  error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, use with twoWayOne2One");

}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingOne2One::push(const char *name,const char *type,void *&to,const char *datatype)
{
    CfdDatacoupling::push(name,type,to,datatype);
/*
    if(strcmp(datatype,"double") == 0)
        push_mpi<double>(name, type, to, world);
    else if(strcmp(datatype,"int") == 0)
        push_mpi<int>(name, type, to, world);
    else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, valid datatypes are 'int' and double'");
    */
  error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, use with twoWayOne2One");

}
