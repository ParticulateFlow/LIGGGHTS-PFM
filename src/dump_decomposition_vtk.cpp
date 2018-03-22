/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include "dump_decomposition_vtk.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define BIG      1.0e30

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::DumpDecompositionVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5)
    error->all(FLERR,"Illegal dump decomposition command");

  format_default = NULL;

  //number of properties written out in one line with buff
  size_one=1;  //dont use buff

  lasttimestep=-1;

  //NP length of data to write

  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;

  //NP data to write to

  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

DumpDecompositionVTK::~DumpDecompositionVTK()
{
  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::init_style()
{
  if (domain->triclinic == 1)
    error->all(FLERR,"Can not perform dump decomposition for triclinic box");
  if (binary)
    error->all(FLERR,"Can not perform dump decomposition in binary mode");

  // default format not needed

  delete [] format;
  format = new char[150];

  // setup function ptrs

  header_choice = &DumpDecompositionVTK::header_item;
  pack_choice = &DumpDecompositionVTK::pack_item;
  write_choice = &DumpDecompositionVTK::write_item;

  openfile();

  //NP realloc since domain decomposition may have changed
  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;
  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump decomposition' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpDecompositionVTK::count()
{
  if (multiproc == 0 && comm->me != 0) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::pack(int *ids)
{
   (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::header_item(bigint ndump)
{
  if (multiproc == 0 && comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/gran/VTK export\nASCII\n");
}

void DumpDecompositionVTK::footer_item()
{
  return;

}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::pack_item()
{
  //NP let the processors write the data
  //NP do this stuff here, since called by all processors

  //NP x-direction
  xdata[0] = -BIG;
  if(comm->myloc[0] == 0) xdata[0] = domain->sublo[0];
  for(int i = 0; i < comm->procgrid[0]; i++)
  {
      xdata[i+1] = -BIG;
      if(comm->myloc[0] == i) xdata[i+1] = domain->subhi[0];
  }

  //NP y-direction
  ydata[0] = -BIG;
  if(comm->myloc[1] == 0) ydata[0] = domain->sublo[1];
  for(int i = 0; i < comm->procgrid[1]; i++)
  {
      ydata[i+1] = -BIG;
      if(comm->myloc[1] == i) ydata[i+1] = domain->subhi[1];
  }

  //NP z-direction
  zdata[0] = -BIG;
  if(comm->myloc[2] == 0) zdata[0] = domain->sublo[2];
  for(int i = 0; i < comm->procgrid[2]; i++)
  {
      zdata[i+1] = -BIG;
      if(comm->myloc[2] == i) zdata[i+1] = domain->subhi[2];
  }

  /*NL*///if (screen) fprintf(screen,"proc %d: xdata %f %f %f\n",comm->me,xdata[0],xdata[1],xdata[2]);

  //NP perform allreduce to get correct boundaries

  MPI_Allreduce(xdata,xdata_all,len[0],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(ydata,ydata_all,len[1],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(zdata,zdata_all,len[2],MPI_DOUBLE,MPI_MAX,world);

  /*NL*///if (screen) fprintf(screen,"proc %d: xdata_all %f %f %f\n",comm->me,xdata_all[0],xdata_all[1],xdata_all[2]);

  return;
}

/* ---------------------------------------------------------------------- */

void DumpDecompositionVTK::write_item(int n, double *mybuf)
{
  //NP only proc 0 writes if multiproc=0
  if (multiproc == 0 && comm->me!=0) return;

  //ensure it is only written once in multi-proc (work-around)
  if(lasttimestep==update->ntimestep)return;
  lasttimestep=update->ntimestep;

  //write the data
  fprintf(fp,"DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\n",len[0],len[1],len[2]);

  if(multiproc == 1) {
    fprintf(fp,"X_COORDINATES 2 float\n");
    fprintf(fp,"%f %f\n",domain->sublo[0], domain->subhi[0]);

    fprintf(fp,"Y_COORDINATES 2 float\n");
    fprintf(fp,"%f %f\n",domain->sublo[1], domain->subhi[1]);

    fprintf(fp,"Z_COORDINATES 2 float\n");
    fprintf(fp,"%f %f\n",domain->sublo[2], domain->subhi[2]);
  } else {
    fprintf(fp,"X_COORDINATES %d float\n",len[0]);
    for (int i = 0; i < len[0]; i++)
       fprintf(fp,"%f ",xdata_all[i]);
    fprintf(fp,"\n");

    fprintf(fp,"Y_COORDINATES %d float\n",len[1]);
    for (int i = 0; i < len[1]; i++)
       fprintf(fp,"%f ",ydata_all[i]);
    fprintf(fp,"\n");

    fprintf(fp,"Z_COORDINATES %d float\n",len[2]);
    for (int i = 0; i < len[2]; i++)
       fprintf(fp,"%f ",zdata_all[i]);
    fprintf(fp,"\n");
  }

  //footer not needed
}
