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
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "set_sph.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "pair.h"
#include "random_park.h"
#include "math_extra.h"
#include "error.h"
#include "modify.h" 
#include "fix_property_atom.h"
#include "sph_kernels.h"
#include "fix_sph.h"

using namespace LAMMPS_NS;

enum{ATOM,GROUP,REGION};
enum{TYPE,X,Y,Z,VX,VY,VZ,DIAMETER,DENSITY,MASS};

/* ---------------------------------------------------------------------- */

SetSph::SetSph(LAMMPS *lmp) : Pointers(lmp)
{
  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

void SetSph::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR, "SetSph command before simulation box is defined");
  if (atom->natoms == 0)
    error->all(FLERR, "SetSph command with no atoms existing");
  if (narg < 3) error->all(FLERR, "Illegal set command");

  // style and ID

  if (strcmp(arg[0],"atom") == 0) style = ATOM;
  else if (strcmp(arg[0],"group") == 0) style = GROUP;
  else if (strcmp(arg[0],"region") == 0) style = REGION;
  else error->all(FLERR, "Illegal set command");

  int n = strlen(arg[1]) + 1;
  id = new char[n];
  strcpy(id,arg[1]);
  select = NULL;

  // loop over keyword/value pairs
  // call appropriate routine to reset attributes

  if (comm->me == 0 && screen) fprintf(screen,"Setting atom values ...\n");

  int allcount,origarg;

  int iarg = 2;
  while (iarg < narg) {
    count = 0;
    origarg = iarg;

    if (strcmp(arg[iarg],"type") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      ivalue = atoi(arg[iarg+1]);
      if (ivalue <= 0 || ivalue > atom->ntypes)
	error->all(FLERR, "Invalid value in set command");
      set(TYPE);
      iarg += 2;
    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(X);
      iarg += 2;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(Y);
      iarg += 2;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(Z);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VX);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      set(VZ);
      iarg += 2;
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->radius_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
      set(DIAMETER);
      iarg += 2;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->density_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
      set(DENSITY);
      iarg += 2;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      dvalue = atof(arg[iarg+1]);
      if (!atom->rmass_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
      set(MASS);
      iarg += 2;
    } else if (strcmp(arg[iarg],"sphkernel") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal set command");
      kernel_style = NULL;
      kernel_style = new char[strlen(arg[iarg+1])+1];
      strcpy(kernel_style,arg[iarg+1]);
/*      if (!atom->rmass_flag)
        error->all(FLERR, "Cannot set this attribute for this atom style");
*/
      // check uniqueness of kernel IDs

      int flag = SPH_KERNEL_NS::sph_kernels_unique_id();
      if(flag < 0) error->all(FLERR, "Cannot proceed, sph kernels need unique IDs, check all sph_kernel_* files");

      // get kernel id

      kernel_id = SPH_KERNEL_NS::sph_kernel_id(kernel_style);
      if(kernel_id < 0) error->all(FLERR, "Illegal pair_style sph command, unknown sph kernel");

      // set kernel_id in all sph fixes

      if (comm->me == 0 && screen) {
        fprintf(screen,"Setting undefined fix_sph kernel IDs ...\n");
        fprintf(screen,"  Sph styles with undefined kernel_id found: \n");
      }
      for (int ifix = 0; ifix < modify->nfix; ifix++)
      {
        if (strstr(modify->fix[ifix]->style,"sph")) {
          if (((FixSph *)(modify->fix[ifix]))->kernel_flag && ((FixSph *)(modify->fix[ifix]))->get_kernel_id() < 0) {
            if (comm->me == 0 && screen) fprintf(screen,"  Fix style = %s\n",modify->fix[ifix]->style);
            ((FixSph *)(modify->fix[ifix]))->set_kernel_id(kernel_id);
            count++;
          }
        }
      }

      iarg += 2;
    } else error->all(FLERR, "Illegal set command");

    // statistics

    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);

    if (comm->me == 0) {
      if (screen) fprintf(screen,"  %d settings on %d procs made for %s\n",
			  allcount,comm->nprocs,arg[origarg]);
      if (logfile) fprintf(logfile,"  %d settings on %d procs made for %s\n",
			   allcount,comm->nprocs,arg[origarg]);
    }
  }

  // free local memory

  delete [] id;
  delete [] select;
  delete [] kernel_style;
}

/* ----------------------------------------------------------------------
   select atoms according to ATOM, GROUP, REGION style
   n = nlocal or nlocal+nghost depending on keyword
------------------------------------------------------------------------- */

void SetSph::selection(int n)
{
  delete [] select;
  select = new int[n];

  if (style == ATOM) {
    if (atom->tag_enable == 0)
      error->all(FLERR, "Cannot use set atom with no atom IDs defined");
    int idatom = atoi(id);

    int *tag = atom->tag;
    for (int i = 0; i < n; i++)
      if (idatom == tag[i]) select[i] = 1;
      else select[i] = 0;

  } else if (style == GROUP) {
    int igroup = group->find(id);
    if (igroup == -1) error->all(FLERR, "Could not find set group ID");
    int groupbit = group->bitmask[igroup];

    int *mask = atom->mask;
    for (int i = 0; i < n; i++)
      if (mask[i] & groupbit) select[i] = 1;
      else select[i] = 0;

  } else {
    int iregion = domain->find_region(id);
    if (iregion == -1) error->all(FLERR, "SetSph region ID does not exist");

    double **x = atom->x;
    for (int i = 0; i < n; i++)
      if (domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2]))
	select[i] = 1;
      else select[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set an owned atom property directly
------------------------------------------------------------------------- */

void SetSph::set(int keyword)
{
  selection(atom->nlocal);

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++)
    if (select[i]) {
      if (keyword == TYPE) atom->type[i] = ivalue;
      else if (keyword == X) atom->x[i][0] = dvalue;
      else if (keyword == Y) atom->x[i][1] = dvalue;
      else if (keyword == Z) atom->x[i][2] = dvalue;
      else if (keyword == VX) atom->v[i][0] = dvalue;
      else if (keyword == VY) atom->v[i][1] = dvalue;
      else if (keyword == VZ) atom->v[i][2] = dvalue;
      else if (keyword == MASS) atom->rmass[i] = dvalue;
      else if (keyword == DIAMETER) atom->radius[i] = 0.5 * dvalue;
      else if (keyword == DENSITY) atom->density[i] = dvalue;
      count++;
    }
}
