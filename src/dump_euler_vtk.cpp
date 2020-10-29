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
#include "dump_euler_vtk.h"
#include "fix_ave_euler.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"
#include <stdint.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpEulerVTK::DumpEulerVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg),
  fix_euler_(0),
  n_calls_(0),
  n_all_(0),
  n_all_max_(0),
  buf_all_(0)
{
  if (narg < 5)
    error->all(FLERR,"Illegal dump euler/vtk command");

  fix_euler_name_ = NULL;
  cell_center=true;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ave_euler") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR,"missing argument");
      ++iarg;
      int n = strlen(arg[iarg]) + 1;
      fix_euler_name_ = new char[n];
      strcpy(fix_euler_name_,arg[iarg]);
      ++iarg;
    } else if (strcmp(arg[iarg], "cell_center") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR,"missing argument");
      ++iarg;
      if (strcmp(arg[iarg], "yes") == 0) {
        cell_center = true;
      } else if (strcmp(arg[iarg], "no") == 0) {
        cell_center = false;
      } else {
        error->all(FLERR,"expecting 'yes' or 'no' after cell_center");;
      }
      ++iarg;
    } else {
      error->all(FLERR,"illegal argument");;
    }

  }

  // CURRENTLY ONLY PROC 0 writes

  format_default = NULL;
}

/* ---------------------------------------------------------------------- */

DumpEulerVTK::~DumpEulerVTK()
{
  delete fix_euler_name_;
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::init_style()
{
  if (fix_euler_name_) {
    fix_euler_ = static_cast<FixAveEuler*>(modify->find_fix_id(fix_euler_name_));
  }
  if(!fix_euler_) {
    if (fix_euler_name_) error->warning(FLERR, "dump euler/vtk failed to find named fix ave/euler*");
    fix_euler_ = static_cast<FixAveEuler*>(modify->find_fix_style("ave/euler",0));
  }
  if(!fix_euler_)
    error->all(FLERR,"Illegal dump euler/vtk command, need a fix ave/euler");

  // multifile=1;             // 0 = one big file, 1 = one file per timestep
  // multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multifile != 1)
    error->all(FLERR,"You should use a filename like 'dump*.vtk' for the 'dump euler/vtk' command to produce one file per time-step");
  if (multiproc != 0)
    error->all(FLERR,"Your 'dump euler/vtk' command is writing one file per processor, where all the files contain the same data");

//  if (domain->triclinic == 1)
//    error->all(FLERR,"Cannot dump VTK files for triclinic box");
  if (binary)
    error->all(FLERR,"Cannot dump VTK files in binary mode");

  if (cell_center) {
    // node center (3), av vel (3), volume fraction, radius, pressure, normal and shear stresses (6)
    size_one = 15;
  } else {
    if (fix_euler_->is_parallel())
      error->all(FLERR,"Cannot dump cell information for fix ave/euler running parallel");

    if (strcmp(fix_euler_->style,"ave/euler") == 0) {
      // boxlo (3), av vel (3), volume fraction, radius, pressure, normal and shear stresses (6)
      size_one = 15;
    } else if (strcmp(fix_euler_->style,"ave/euler/region") == 0) {
      // node points (24), av vel (3), volume fraction, radius, pressure, normal and shear stresses (6)
      size_one = 36;
    } else {
      error->all(FLERR,"Unable to dump cells to VTK files");
    }
  }

  delete [] format;
}

/* ---------------------------------------------------------------------- */

int DumpEulerVTK::modify_param(int narg, char **arg)
{
  error->warning(FLERR,"dump_modify keyword is not supported by 'dump euler/vtk' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::write_header(bigint ndump)
{
  write_header_ascii(ndump);
}

void DumpEulerVTK::write_header_ascii(bigint ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/VTK export\nASCII\n");
}

/* ---------------------------------------------------------------------- */

int DumpEulerVTK::count()
{
  n_calls_ = 0;
  n_all_ = 0;
  return fix_euler_->ncells_pack();
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::pack(int *ids)
{
  int m = 0;

  // have to stick with this order (all per-element props)
  // as multiple procs pack

  int ncells = fix_euler_->ncells_pack();

  for(int i = 0; i < ncells; i++)
  {
    if (cell_center) {
      buf[m++] = fix_euler_->cell_center(i,0);
      buf[m++] = fix_euler_->cell_center(i,1);
      buf[m++] = fix_euler_->cell_center(i,2);
    } else {
      if (strcmp(fix_euler_->style,"ave/euler") == 0) {
        buf[m++] = fix_euler_->cell_center(i,0)-0.5*fix_euler_->cell_size(0);
        buf[m++] = fix_euler_->cell_center(i,1)-0.5*fix_euler_->cell_size(1);
        buf[m++] = fix_euler_->cell_center(i,2)-0.5*fix_euler_->cell_size(2);
      } else if (strcmp(fix_euler_->style,"ave/euler/region") == 0) {
        double points[24];
        fix_euler_->cell_points(i,points);
        for (int j=0; j < 24; ++j) {
          buf[m++] = points[j];
        }
      } else {
        error->all(FLERR,"Unable to dump cells to VTK files");
      }
    }

    buf[m++] = fix_euler_->cell_v_av(i,0);
    buf[m++] = fix_euler_->cell_v_av(i,1);
    buf[m++] = fix_euler_->cell_v_av(i,2);

    buf[m++] = fix_euler_->cell_vol_fr(i);
    buf[m++] = fix_euler_->cell_radius(i);
    buf[m++] = fix_euler_->cell_pressure(i);

    buf[m++] = fix_euler_->cell_stress(i,0);
    buf[m++] = fix_euler_->cell_stress(i,1);
    buf[m++] = fix_euler_->cell_stress(i,2);
    buf[m++] = fix_euler_->cell_stress(i,3);
    buf[m++] = fix_euler_->cell_stress(i,4);
    buf[m++] = fix_euler_->cell_stress(i,5);
  }
  return ;
}

/* ---------------------------------------------------------------------- */

void DumpEulerVTK::write_data(int n, double *mybuf)
{
    //only proc 0 writes
    if (comm->me != 0) return;

    n_calls_++;

    // grow buffer if necessary
    if(n_all_+n*size_one > n_all_max_)
    {
        n_all_max_ = n_all_ + n*size_one;
        memory->grow(buf_all_,n_all_max_,"DumpEulerVTK:buf_all_");
    }

    // copy to buffer
    vectorCopyN(mybuf,&(buf_all_[n_all_]),n*size_one);
    n_all_ += n*size_one;

    // write on last call
    if(n_calls_ == comm->nprocs)
        write_data_ascii(n_all_/size_one,buf_all_);
}

void DumpEulerVTK::write_data_ascii(int n, double *mybuf)
{

  int m, buf_pos;

  // n is the number of elements

  /*NL*///if (screen) fprintf(screen,"WRITING ITEM at step %d proc %d with n %d\n",update->ntimestep,comm->me,n);
  /*NL*///error->one(FLERR,"end");

  // write point data
  m = 0;
  buf_pos = 0;
  if (cell_center) {
    fprintf(fp,"DATASET POLYDATA\nPOINTS %d float\n",n);
    for (int i = 0; i < n; i++)
    {
      fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += size_one ;
    }
    buf_pos += 3;

    // write polygon data
    fprintf(fp,"VERTICES %d %d\n",n,2*n);
    for (int i = 0; i < n; i++)
    {
        fprintf(fp,"%d %d\n",1,i);
    }

    // write point data header
    fprintf(fp,"POINT_DATA %d\n",n);

  } else {
    if (strcmp(fix_euler_->style,"ave/euler") == 0) {
      // dump a rectilinear grid of cuboidal cell
      int nx = fix_euler_->ncells(0);
      int ny = fix_euler_->ncells(1);
      int nz = fix_euler_->ncells(2);

      if (nx*ny*nz != n)
        error->all(FLERR,"number of cells does not match nx*ny*nz != n");

      fprintf(fp,"DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\n",nx+1,ny+1,nz+1);

      // x coordinates
      m = 0;
      fprintf(fp,"X_COORDINATES %d float\n", nx+1);
      for (int i = 0; i < nx; i++) {
        fprintf(fp,"%f ",mybuf[m]);
        m += size_one ;
      }
      fprintf(fp,"%f\n", mybuf[m-size_one]+fix_euler_->cell_size(0));

      // y coordinates
      m = 1;
      fprintf(fp,"Y_COORDINATES %d float\n", ny+1);
      for (int i = 0; i < ny; i++) {
        fprintf(fp,"%f ",mybuf[m]);
        m += nx*size_one ;
      }
      fprintf(fp,"%f\n", mybuf[m-nx*size_one]+fix_euler_->cell_size(1));

      // z coordinates
      m = 2;
      fprintf(fp,"Z_COORDINATES %d float\n", nz+1);
      for (int i = 0; i < nz; i++) {
        fprintf(fp,"%f ",mybuf[m]);
        m += nx*ny*size_one ;
      }
      fprintf(fp,"%f\n", mybuf[m-nx*ny*size_one]+fix_euler_->cell_size(2));

      buf_pos += 3;
      // write cell data header
      fprintf(fp,"CELL_DATA %d\n",n);

    } else if (strcmp(fix_euler_->style,"ave/euler/region") == 0) {
      // dump an unstructured grid of hexahedral cells
      fprintf(fp,"DATASET UNSTRUCTURED_GRID\nPOINTS %d float\n",n*8);

      for (int i = 0; i < n; i++) {
        fprintf(fp,"%f %f %f %f %f %f\n%f %f %f %f %f %f\n"
                   "%f %f %f %f %f %f\n%f %f %f %f %f %f\n",
                    mybuf[m],   mybuf[m+1], mybuf[m+2],
                    mybuf[m+3], mybuf[m+4], mybuf[m+5],
                    mybuf[m+6], mybuf[m+7], mybuf[m+8],
                    mybuf[m+9], mybuf[m+10],mybuf[m+11],
                    mybuf[m+12],mybuf[m+13],mybuf[m+14],
                    mybuf[m+15],mybuf[m+16],mybuf[m+17],
                    mybuf[m+18],mybuf[m+19],mybuf[m+20],
                    mybuf[m+21],mybuf[m+22],mybuf[m+23]);
        m += size_one ;
      }

      fprintf(fp,"CELLS %d %d\n",n,n*9);
      for (int i = 0; i < n; i++) {
        fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",
                i*8,i*8+1,i*8+2,i*8+3,i*8+4,i*8+5,i*8+6,i*8+7);
      }

      fprintf(fp,"CELL_TYPES %d\n",n);
      for (int i = 0; i < n; i++) {
        fprintf(fp,"12\n"); // hexahedron = 12
      }

      buf_pos += 24;
      // write cell data header
      fprintf(fp,"CELL_DATA %d\n",n);
    } else {
      error->all(FLERR,"Unable to dump cells to VTK files");
    }
  }

  // write cell data

  fprintf(fp,"VECTORS v_avg float\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
     fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
     m += size_one;
  }
  buf_pos += 3;

  fprintf(fp,"SCALARS volumefraction float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  fprintf(fp,"SCALARS radius float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  fprintf(fp,"SCALARS pressure float 1\nLOOKUP_TABLE default\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      fprintf(fp,"%f\n",mybuf[m]);
      m += size_one;
  }
  buf_pos++;

  fprintf(fp,"VECTORS sigma float\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      // normal stresses
      fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += size_one;
  }
  buf_pos += 3;

  fprintf(fp,"VECTORS tau float\n");
  m = buf_pos;
  for (int i = 0; i < n; i++)
  {
      // shear stresses
      // note: internal order of stresses is (xx, yy, zz, xy, xz, yz)
      //       Voigt notation is             (xx, yy, zz, yz, xz, xy)
      fprintf(fp,"%f %f %f\n",mybuf[m],mybuf[m+1],mybuf[m+2]);
      m += size_one;
  }
  buf_pos += 3;
  // footer not needed
  // if would be needed, would do like in dump stl
}
