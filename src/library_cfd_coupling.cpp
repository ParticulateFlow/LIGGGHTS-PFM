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
#include "library_cfd_coupling.h"
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
#include "universe.h"

using namespace LAMMPS_NS;

#define LMP_GROW_DELTA 11000
/*NL*/ #define LMP_OF_DEBUG false

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag(void *ptr, int iworld)
{
  LAMMPS *lmp = (LAMMPS *) ptr;
  int tag_max = lmp->atom->tag_max();
  // bcast from root of iworld to all procs in universe
  if (lmp->universe->existflag == 1)
    MPI_Bcast(&tag_max, 1, MPI_INT, lmp->universe->root_proc[iworld], lmp->universe->uworld);
  return tag_max;
}

/* ---------------------------------------------------------------------- */

int liggghts_get_maxtag_ms(void *ptr, int iworld)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;

  int tag_max_body = fix_ms->tag_max_body();
  // bcast from root of iworld to all procs in universe
  if (lmp->universe->existflag == 1)
    MPI_Bcast(&tag_max_body, 1, MPI_INT, lmp->universe->root_proc[iworld], lmp->universe->uworld);
  return tag_max_body;
}

/* ---------------------------------------------------------------------- */

int liggghts_get_ntypes_ms(void *ptr, int iworld)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;

  int ntypes_ms = fix_ms->ntypes();
  // bcast from root of iworld to all procs in universe
  if (lmp->universe->existflag == 1)
    MPI_Bcast(&ntypes_ms, 1, MPI_INT, lmp->universe->root_proc[iworld], lmp->universe->uworld);
  return ntypes_ms;
}

/* ---------------------------------------------------------------------- */

double* liggghts_get_vclump_ms(void *ptr, int iworld)
{
  // currently no possibility to delete multisphere bodies
  // so just return # of bodies

  LAMMPS *lmp = (LAMMPS *) ptr;
  FixMultisphere *fix_ms = static_cast<FixMultisphere*>(lmp->modify->find_fix_style_strict("multisphere",0));
  if(!fix_ms) return 0;
  if (lmp->universe->existflag == 1) lmp->error->all(FLERR,"Not implemented for universe");
  return fix_ms->vclump();
}

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

void data_liggghts_to_of(const char *name, const char *type, void *ptr, void *&data, const char *datatype, int iworld)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->push(name,type,data,datatype,iworld);
}

/* ---------------------------------------------------------------------- */

void data_of_to_liggghts(const char *name,const char *type,void *ptr,void *data,const char* datatype, int iworld)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->pull(name,type,data,datatype,iworld);
}

/* ---------------------------------------------------------------------- */

//NP update region model
void update_region_model(void *ptr, int iworld)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    //FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    //locate_coupling_fix(ptr);
    //CfdRegionmodel *rm = fcfd->rm;

    //NP call region model
    //if(rm) rm->rm_update();
    lmp->error->all(FLERR,"Region model update not implemented aborting.");
}

/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}
/* ---------------------------------------------------------------------- */

void allocate_external_int(int    **&data, int len2,const char *keyword,int    initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,len1,initvalue);
}

/* ---------------------------------------------------------------------- */

void allocate_external_double(double **&data, int len2,const char* keyword,double initvalue,void *ptr)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->allocate_external(data,len2,keyword,initvalue);
}

/* ---------------------------------------------------------------------- */

//NP check if all requested quantities have been communicated
//NP    since last call of this function
void check_datatransfer(void *ptr, int iworld)
{
    //LAMMPS *lmp = (LAMMPS *) ptr;
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    fcfd->get_dc()->check_datatransfer();
}

/* ---------------------------------------------------------------------- */

double** o2o_liggghts_get_boundingbox(void *ptr)
{
    LAMMPS *lmp = (LAMMPS *) ptr;
    double** bbox = new double*[2];
    bbox[0] = lmp->domain->sublo;
    bbox[1] = lmp->domain->subhi;
    return bbox;
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
    FixCfdCoupling* fcfd = (FixCfdCoupling*)locate_coupling_fix(ptr);
    CfdDatacouplingOne2One* dc = static_cast<CfdDatacouplingOne2One*>(fcfd->get_dc());

    if(strcmp(datatype,"double") == 0)
        dc->pull_mpi<double>(name, type, data, ids, ncollected);
    else if(strcmp(datatype,"int") == 0)
        dc->pull_mpi<int>(name, type, data, ids, ncollected);
 //   else error->one(FLERR,"Illegal call to CfdDatacouplingOne2One::pull, valid datatypes are 'int' and double'");

}

