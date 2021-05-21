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
   Contributing author:
   Anton Gladky(TU Bergakademie Freiberg), gladky.anton@gmail.com
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK
#include <string.h>
#include "dump_atom_vtk.h"
#include "atom.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DumpATOMVTK::DumpATOMVTK(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump command");
  if (binary || multiproc) error->all(FLERR,"Invalid dump filename");

  sort_flag = 1;
  sortcol = 0;

  size_one = 17;

  char *str = (char *) "%d %g %g %g";
  int n = strlen(str) + 1;
  format_default = new char[n];
  strcpy(format_default,str);
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::init_style()
{

}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::write_header(bigint /* n */)
{
}

/* ---------------------------------------------------------------------- */

int DumpATOMVTK::count()
{
  n_calls_ = 0;

  return Dump::count();
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::pack(int *ids)
{
  int n = 0;
  int m = 0;

  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **o = atom->omega;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (ids) ids[n++] = tag[i];

      double massTemp;

      if (rmass) {
        massTemp=rmass[i];
      } else {
        massTemp=mass[type[i]];
      }

      int me = comm->me;

      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];

      buf[m++] = atom->radius[i];
      buf[m++] = massTemp;
      buf[m++] = static_cast<double>(tag[i]);
      buf[m++] = static_cast<double>(type[i]);

      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];

      buf[m++] = o[i][0];
      buf[m++] = o[i][1];
      buf[m++] = o[i][2];

      buf[m++] = f[i][0];
      buf[m++] = f[i][1];
      buf[m++] = f[i][2];

      buf[m++] = static_cast<double>(me);
    }
  }

  setFileCurrent();
  tmpEXP.setFileName(filecurrent);
  return;
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::write_data(int n, double *mybuf)
{
  if (comm->me != 0) return;
  n_calls_++;
  int m = 0;
  for (int i = 0; i < n; i++) {
    DumpATOMVTK::DataVTK tmpVTKDat(
      V3(mybuf[m+0], mybuf[m+1], mybuf[m+2]),
      mybuf[m+3], mybuf[m+4], static_cast<int> (mybuf[m+5]), static_cast<int> (mybuf[m+6]),
      V3(mybuf[m+7], mybuf[m+8], mybuf[m+9]),
      V3(mybuf[m+10], mybuf[m+11], mybuf[m+12]),
      V3(mybuf[m+13], mybuf[m+14], mybuf[m+15]),
      static_cast<int> (mybuf[m+16]));
    tmpEXP.add(tmpVTKDat);
    m += size_one;
  }

  if(n_calls_ == comm->nprocs) {
    tmpEXP.writeSER();
    tmpEXP.clear();
    delete [] filecurrent;
  }
}

/* ---------------------------------------------------------------------- */

DumpATOMVTK::DataVTK::DataVTK(V3 Pos, double Rad, double Mass, int Id, int Type,
  V3 VelL, V3 VelA, V3 Force, int proc) {
  _Pos = Pos;
  _Rad = Rad;
  _Mass = Mass;
  _Id = Id;
  _Type = Type;
  _VelL = VelL;
  _VelA = VelA;
  _Force = Force;
  _proc = proc;
}
/* ---------------------------------------------------------------------- */

std::string  DumpATOMVTK::DataVTK::serialize() {
  std::string tmp;
  std::ostringstream stringStream;

  stringStream <<_Pos[0]<<' '<<_Pos[1]<<' '<<_Pos[2]<<' '<<_Rad<<' '<<_Mass<<' '<<_Id
   <<' '<<_Type<<' '
   <<_VelL[0]<<' '<<_VelL[1]<<' '<<_VelL[2]<<' '
   <<_VelA[0]<<' '<<_VelA[1]<<' '<<_VelA[2]<<' '
   <<_Force[0]<<' '<<_Force[1]<<' '<<_Force[2]<<' '<<_proc<<'\n';

  tmp  = stringStream.str();
  return tmp;
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::vtkExportData::add(DumpATOMVTK::DataVTK & d) {
  vtkData.push_back(d);
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::vtkExportData::setFileName(const char * fileName) {
  _fileName = fileName;
  _setFileName = true;
}

/* ---------------------------------------------------------------------- */

DumpATOMVTK::vtkExportData::vtkExportData() {
  _setFileName=false;
}
/* ---------------------------------------------------------------------- */

int DumpATOMVTK::vtkExportData::size() {
  return vtkData.size();
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::vtkExportData::writeSER() {

  vtkSmartPointer<vtkPoints>  spheresPos = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("radii");

  vtkSmartPointer<vtkDoubleArray> spheresMass = vtkSmartPointer<vtkDoubleArray>::New();
  spheresMass->SetNumberOfComponents(1);
  spheresMass->SetName("mass");

  vtkSmartPointer<vtkIntArray> spheresId = vtkSmartPointer<vtkIntArray>::New();
  spheresId->SetNumberOfComponents(1);
  spheresId->SetName("id");

  vtkSmartPointer<vtkIntArray> spheresType = vtkSmartPointer<vtkIntArray>::New();
  spheresType->SetNumberOfComponents(1);
  spheresType->SetName("type");

  vtkSmartPointer<vtkIntArray> spheresProc = vtkSmartPointer<vtkIntArray>::New();
  spheresProc->SetNumberOfComponents(1);
  spheresProc->SetName("proc");

  vtkSmartPointer<vtkDoubleArray> spheresVelL = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelL->SetNumberOfComponents(3);
  spheresVelL->SetName("velocity_lin");

  vtkSmartPointer<vtkDoubleArray> spheresVelA = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelA->SetNumberOfComponents(3);
  spheresVelA->SetName("velocity_ang");

  vtkSmartPointer<vtkDoubleArray> spheresForce = vtkSmartPointer<vtkDoubleArray>::New();
  spheresForce->SetNumberOfComponents(3);
  spheresForce->SetName("force");

  for (unsigned int i=0; i < vtkData.size(); i++) {
    vtkIdType pid[1];
    pid[0] = spheresPos->InsertNextPoint(vtkData[i]._Pos[0], vtkData[i]._Pos[1], vtkData[i]._Pos[2]);
    radii->InsertNextValue(vtkData[i]._Rad);

    double vv[3] = {vtkData[i]._VelL[0], vtkData[i]._VelL[1], vtkData[i]._VelL[2]};
#if (VTK_MAJOR_VERSION > 7) || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION     >= 1)
    spheresVelL->InsertNextTypedTuple(vv);
#else
    spheresVelL->InsertNextTupleValue(vv);
#endif

    double oo[3] = {vtkData[i]._VelA[0], vtkData[i]._VelA[1], vtkData[i]._VelA[2]};
#if (VTK_MAJOR_VERSION > 7) || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION     >= 1)
    spheresVelA->InsertNextTypedTuple(oo);
#else
    spheresVelA->InsertNextTupleValue(oo);
#endif

    double ff[3] = {vtkData[i]._Force[0], vtkData[i]._Force[1], vtkData[i]._Force[2]};
#if (VTK_MAJOR_VERSION > 7) || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION     >= 1)
    spheresForce->InsertNextTypedTuple(ff);
#else
    spheresForce->InsertNextTupleValue(ff);
#endif

    spheresMass->InsertNextValue(vtkData[i]._Mass);

    spheresId->InsertNextValue(vtkData[i]._Id);
    spheresType->InsertNextValue(vtkData[i]._Type);
    spheresProc->InsertNextValue(vtkData[i]._proc);

    spheresCells->InsertNextCell(1,pid);
  }

  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();

  spheresUg->SetPoints(spheresPos);
  spheresUg->SetCells(VTK_VERTEX, spheresCells);
  spheresUg->GetPointData()->AddArray(radii);
  spheresUg->GetPointData()->AddArray(spheresId);
  spheresUg->GetPointData()->AddArray(spheresType);
  spheresUg->GetPointData()->AddArray(spheresProc);
  spheresUg->GetPointData()->AddArray(spheresMass);
  spheresUg->GetPointData()->AddArray(spheresVelL);
  spheresUg->GetPointData()->AddArray(spheresVelA);
  spheresUg->GetPointData()->AddArray(spheresForce);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION < 6
  writer->SetInput(spheresUg);
#else
  writer->SetInputData(spheresUg);
#endif
  writer->SetFileName(_fileName);
  writer->Write();
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::vtkExportData::show() {
  for (unsigned int i=0; i < vtkData.size(); i++) {
    std::cerr << vtkData[i].serialize();
  }
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::vtkExportData::clear() {
  vtkData.clear();
}

/* ---------------------------------------------------------------------- */

void DumpATOMVTK::setFileCurrent() {
  if (multifile == 0) filecurrent = filename;
  else {
    filecurrent = new char[strlen(filename) + 16];
    char *ptr = strchr(filename,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filename,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filename,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }
}
#endif
