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

#include <mpi.h>
#include <string.h>
#include "library_cfd_coupling_one2one.h"
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix_cfd_coupling.h"
#include "fix_multisphere.h"
#include "cfd_regionmodel.h"
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "cfd_datacoupling.h"
#include "cfd_datacoupling_one2one.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000

double** o2o_liggghts_get_boundingbox(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    double** bbox = new double*[2];
    bbox[0] = lmp->domain->sublo;
    bbox[1] = lmp->domain->subhi;
    return bbox;
}

/* ---------------------------------------------------------------------- */

int o2o_liggghts_get_maxtag(void *ptr)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  return lmp->atom->tag_max();
}

/* ---------------------------------------------------------------------- */

/*
int o2o_liggghts_get_maxtag_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->tag_max_body();
}
*/
/* ---------------------------------------------------------------------- */
/*
int o2o_liggghts_get_ntypes_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->ntypes();
}
*/
/* ---------------------------------------------------------------------- */
/*
double* o2o_liggghts_get_vclump_ms(void *ptr)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  return fix_ms->vclump();
}
*/
/* ---------------------------------------------------------------------- */

void* o2o_locate_coupling_fix(void *ptr)
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

void o2o_data_liggghts_to_of
(
    const char *name,
    const char *type,
    void *ptr,
    void *&data,
    const char *datatype
)
{
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    CfdDatacouplingOne2One* dc = static_cast<CfdDatacouplingOne2One*>(fcfd->get_dc()); 

    if(strcmp(datatype,"double") == 0)
        dc->push_mpi<double>(name, type, data);
    else if(strcmp(datatype,"int") == 0)
        dc->push_mpi<int>(name, type, data);
//    else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull_mpi, valid datatypes are 'int' and double'");
}

/* ---------------------------------------------------------------------- */

void o2o_data_of_to_liggghts
(
    const char *name,
    const char *type,
    void *ptr,
    void *data,
    const char* datatype,
    const int* ids,
    const int ncollected
)
{
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    CfdDatacouplingOne2One* dc = static_cast<CfdDatacouplingOne2One*>(fcfd->get_dc());

    if(strcmp(datatype,"double") == 0)
        dc->pull_mpi<double>(name, type, data, ids, ncollected);
    else if(strcmp(datatype,"int") == 0)
        dc->pull_mpi<int>(name, type, data, ids, ncollected);
 //   else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, valid datatypes are 'int' and double'");

}

/* ---------------------------------------------------------------------- */

//NP update region model
void o2o_update_rm(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    //FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
//    o2o_locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //NP call region model
    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void o2o_allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void o2o_allocate_external_int(int    **&data, int len2,const char *keyword,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void o2o_allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void o2o_allocate_external_double(double **&data, int len2, int len1, const char* keyword,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

//NP check if all requested quantities have been communicated
//NP    since last call of this function
void o2o_check_datatransfer(void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)o2o_locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}
