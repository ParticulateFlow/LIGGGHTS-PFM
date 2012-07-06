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

/* ----------------------------------------------------------------------
   CFD-DEM Coupling Stuff
------------------------------------------------------------------------- */

#include "mpi.h"
#include "string.h"
#include "library_cfd_coupling.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_cfd_coupling.h"
/*NL*/#include "fix_rigid_multisphere.h"
#include "cfd_regionmodel.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "cfd_datacoupling.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000
/*NL*/ #define LMP_OF_DEBUG false

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->tag_max();
}

/* ---------------------------------------------------------------------- */

/*NL*/int liggghts_get_maxtag_ms(void *ptr)
/*NL*/{
/*NL*/  // currently no possibility to delete multisphere bodies
/*NL*/  // so just return # of bodies
/*NL*/
/*NL*/  LAMMPS *lmp = (LAMMPS *) ptr;
/*NL*/  FixRigidMultisphere *frm = static_cast<FixRigidMultisphere*>(lmp->modify->find_fix_style_strict("rigid/multisphere",0));
/*NL*/  if(!frm) return 0;
/*NL*/  return frm->tag_max_body();
/*NL*/}

/* ---------------------------------------------------------------------- */

void* locate_coupling_fix(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    int ifix = -1;
    for(int i=0;i<lmp->modify->nfix;i++)
      if(strcmp(lmp->modify->fix[i]->style,"couple/cfd") == 0)
        ifix = i;

    if(ifix ==-1) lmp->error->all(FLERR,"No fix of style 'couple/cfd' found, aborting.");

    return ((void*)lmp->modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void data_liggghts_to_of(char *name,char *type,void *ptr,void *&data,char* datatype)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->push(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

void data_of_to_liggghts(char *name,char *type,void *ptr,void *data,char* datatype)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->pull(name,type,data,datatype);
}

/* ---------------------------------------------------------------------- */

//NP update region model
void update_rm(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //NP call region model
    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,char *keyword,int    initvalue,void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,char* keyword,double initvalue,void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

//NP check if all requested quantities have been communicated
//NP    since last call of this function
void check_datatransfer(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}
