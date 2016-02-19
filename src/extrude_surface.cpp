/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2016- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#if defined(LAMMPS_VTK) //NP do not use #ifdef here (VS C++ bug)
#include "lmptype.h"
#include "mpi.h"
#include "stdlib.h"
#include "extrude_surface.h"
#include "math_extra.h"
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include "error.h"

using namespace LAMMPS_NS;

#define POLYDATA
/* ---------------------------------------------------------------------- */

ExtrudeSurface::ExtrudeSurface(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world, &me);
  binary = false;
}

/* ---------------------------------------------------------------------- */

ExtrudeSurface::~ExtrudeSurface()
{
}

/* ---------------------------------------------------------------------- */

template<class TReader> vtkDataSet* ExtrudeSurface::read_file(const char*filename)
{
  vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
  reader->SetFileName(filename);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

/* ---------------------------------------------------------------------- */

void ExtrudeSurface::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal extrude_surface command");

  if (me == 0) {
    char *filename = arg[0];
    char *suffix = filename + strlen(filename) - strlen(".vtk");
    vtkDataSet* dset = NULL;

    if (suffix > filename && strcmp(suffix,".vtk") == 0) {
      vtkDataReader *chooser = vtkDataReader::New();
      chooser->SetFileName(filename);
      if (chooser->IsFileUnstructuredGrid()) {
        dset = read_file<vtkUnstructuredGridReader>(filename);
      } else if (chooser->IsFilePolyData()) {
        dset = read_file<vtkPolyDataReader>(filename);
      } else {
        chooser->Delete();
        error->all(FLERR,"extrude_surface requires unstructured grid dataset or polydata");
      }
    } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
      vtkXMLUnstructuredGridReader *chooser = vtkXMLUnstructuredGridReader::New();
      if (!chooser->CanReadFile(filename)) {
        chooser->Delete();
        error->all(FLERR,"extrude_surface cannot read input file");
      }
      dset = read_file<vtkXMLUnstructuredGridReader>(filename);
    } else {
      error->all(FLERR,"extrude_surface: invalid input file");
    }

    // check cell types
    int ncells = dset->GetNumberOfCells();
    for (int i = 0; i < ncells; ++i) {
      if(dset->GetCellType(i) != VTK_POLYGON &&
         dset->GetCellType(i) != VTK_PIXEL &&
         dset->GetCellType(i) != VTK_QUAD) {
        error->all(FLERR,"extrude_surface: input file contains unsupported cell type");
      }
    }

    // check cell data
    if(!dset->GetCellData()->HasArray("face_id")) {
      // create cell data
      vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
      cellData->SetName("face_id");
      cellData->SetNumberOfComponents(1);
      for(int i = 0; i < ncells; ++i) {
        cellData->InsertNextValue(i);
      }
      dset->Update(); // force an update so we can set cell data
      dset->GetCellData()->SetScalars(cellData);
    }

    // triangulate surface
    triangulate(narg, arg, dset);

    // extrude faces to get regions for insertion
    if (narg >= 8) {
      extrude(narg, arg, dset);
    }

    dset->Delete();
  }
}

/* ---------------------------------------------------------------------- */

void ExtrudeSurface::triangulate(int narg, char **arg, vtkDataSet* dset)
{
  vtkNew<vtkTriangleFilter> triangleFilter;
#if VTK_MAJOR_VERSION < 6
  triangleFilter->SetInput(dset);
#else
  triangleFilter->SetInputData(dset);
#endif
  triangleFilter->Update();

  vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
#if VTK_MAJOR_VERSION < 6
  appendFilter->AddInput(triangleFilter->GetOutput());
#else
  appendFilter->AddInputData(triangleFilter->GetOutput());
#endif
  appendFilter->Update();

  // write surface data to file
  if (narg > 3 && strcmp(arg[narg-1],"binary") == 0) binary = true;

  char *filename = arg[2];
  char *suffix = filename + strlen(filename) - strlen(".vtu");

  if (suffix > filename && strcmp(suffix,".vtu") == 0) { // xml format
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->SetFileName(filename);
    writer->Write();
  } else { // legacy format
    vtkNew<vtkUnstructuredGridWriter> writer;
    if (binary) writer->SetFileTypeToBinary();
    else        writer->SetFileTypeToASCII();
    writer->SetInputConnection(appendFilter->GetOutputPort());
    writer->SetHeader("Generated by LIGGGHTS");
    writer->SetFileName(filename);
    writer->Write();
  }
}

/* ---------------------------------------------------------------------- */

void ExtrudeSurface::extrude(int /*narg*/, char **arg, vtkDataSet* dset)
{
  char *filename = arg[3];
  if(strcmp(arg[4],"extrude_length"))
    error->all(FLERR,"extrude_surface: expecting keyword 'extrude_length'");
  double ScaleFactor = atof(arg[5]);
  if(strcmp(arg[6],"min_rad"))
    error->all(FLERR,"extrude_surface: expecting keyword 'min_rad'");
  double atom_radius = atof(arg[7]);

  const double scale_out = atom_radius;
  const double scale_in = atom_radius-ScaleFactor;

  vtkSmartPointer<vtkPolyDataNormals> skinNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION < 6
  skinNormals->SetInput(dset);
#else
  skinNormals->SetInputData(dset);
#endif
  skinNormals->SetFeatureAngle(1.0);

  vtkPolyData *input = skinNormals->GetOutput();
  input->Update();

  vtkIdType numPts, numCells;

  numPts = input->GetNumberOfPoints();
  numCells = input->GetNumberOfCells();

  if (numPts < 1 || numCells < 1) {
    error->all(FLERR,"No data to extrude!");
  }

  vtkSmartPointer<vtkUnstructuredGrid> output = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkPointData *pd = input->GetPointData();
  vtkDataArray *inNormals = NULL;
  vtkPolyData *mesh;
  vtkPoints *inPts;
  vtkCellArray *inPolys;
  vtkIdType inCellId;
  vtkIdType ptId, ncells;
  double x[3], xx[3];
  vtkPoints *newPts;
  vtkCellArray *newHexahedra = NULL;
  vtkIdList *cellPts;
  vtkPointData *outputPD = output->GetPointData();
  inNormals = pd->GetNormals();

  // Build cell data structure.
  mesh = vtkPolyData::New();
  inPts = input->GetPoints();
  inPolys = input->GetPolys();
  mesh->SetPoints(inPts);
  mesh->SetPolys(inPolys);

  if (inPolys->GetNumberOfCells()) {
    mesh->BuildLinks();
  }

  // allocate memory for output.
  output->GetCellData()->CopyNormalsOff();
  output->GetCellData()->CopyAllocate(input->GetCellData(),numCells);

  outputPD->CopyNormalsOff();
  outputPD->CopyAllocate(pd,2*numPts);

  newPts = vtkPoints::New();
  newPts->SetNumberOfPoints(2*numPts);

  ncells = inPolys->GetNumberOfCells();
  newHexahedra = vtkCellArray::New();
  newHexahedra->Allocate(newHexahedra->EstimateSize(ncells,8));

  // copy/add points
  for (ptId=0; ptId < numPts; ++ptId) {
    inPts->GetPoint(ptId, x);
    xx[0] = x[0];xx[1] = x[1];xx[2] = x[2];
    extrude_point_via_normal(x, ptId,inNormals,scale_in);
    extrude_point_via_normal(xx,ptId,inNormals,scale_out);
    newPts->SetPoint(ptId,x);
    newPts->SetPoint(ptId+numPts,xx);
    outputPD->CopyData(pd,ptId,ptId);
    outputPD->CopyData(pd,ptId,ptId+numPts);
  }

  vtkGenericCell *cell = vtkGenericCell::New();

  int pointsleft = 0;
  int pointsneeded = 4;
  double a[3], b[3], c[3];

  for (inCellId=0; inCellId < numCells; ++inCellId) {
    mesh->GetCell(inCellId,cell);
    cellPts = cell->GetPointIds();

    if (cell->GetCellDimension() == 2) {

      newHexahedra->InsertNextCell(8);

      int npoints = cell->GetNumberOfPoints();

      if (npoints < 4) {
        error->all(FLERR,"Unexpected number of points in vtk cell\n");
      } else if (npoints > 4) {

        pointsleft = npoints;
        pointsneeded = 4;
        for (int i = 0; i < npoints && pointsneeded; ++i, --pointsleft) {
          if (pointsleft > pointsneeded) {
            inPts->GetPoint(cellPts->GetId((i-1+npoints)%npoints), a);
            inPts->GetPoint(cellPts->GetId(i), c);
            inPts->GetPoint(cellPts->GetId((i+1)%npoints), b);
            if (collinear(a,b,c))
              continue;
          }
          newHexahedra->InsertCellPoint(cellPts->GetId(i));
          --pointsneeded;
        }

        pointsleft = npoints;
        pointsneeded = 4;
        for (int i = 0; i < npoints && pointsneeded; ++i, --pointsleft) {
          if (pointsleft > pointsneeded) {
            inPts->GetPoint(cellPts->GetId((i-1+npoints)%npoints), a);
            inPts->GetPoint(cellPts->GetId(i), b);
            inPts->GetPoint(cellPts->GetId((i+1)%npoints), c);
            if (collinear(a,b,c))
              continue;
          }
          newHexahedra->InsertCellPoint(cellPts->GetId(i)+numPts);
          --pointsneeded;
        }
      } else {
        for (int i = 0; i < 4; ++i)
          newHexahedra->InsertCellPoint(cellPts->GetId(i));
        for (int i = 0; i < 4; ++i)
          newHexahedra->InsertCellPoint(cellPts->GetId(i)+numPts);
      }
    }
  }

  cell->Delete();

  // copy cell data
  for (int i = 0; i < numCells; ++i) {
    output->GetCellData()->CopyData(input->GetCellData(),i,i);
  }

  // send data to output and release memory
  output->SetPoints(newPts);
  newPts->Delete();
  mesh->Delete();

  output->SetCells(VTK_HEXAHEDRON, newHexahedra);

  output->Squeeze();

  vtkNew<vtkUnstructuredGridWriter> writer;
  if (binary) writer->SetFileTypeToBinary();
  else        writer->SetFileTypeToASCII();
#if VTK_MAJOR_VERSION < 6
  writer->SetInput(output);
#else
  writer->SetInputData(output);
#endif
  writer->SetHeader("Generated by LIGGGHTS");
  writer->SetFileName(filename);
  writer->Write();
}

/* ---------------------------------------------------------------------- */

bool ExtrudeSurface::collinear(double *a, double *b, double *c)
{
  // test if a,b,c are colinear (distance of point (c) to line (a,b))
  double dir1[3], dir2[3], dir3[3], cross[3];
  dir1[0] = c[0] - a[0];
  dir1[1] = c[1] - a[1];
  dir1[2] = c[2] - a[2];
  dir2[0] = c[0] - b[0];
  dir2[1] = c[1] - b[1];
  dir2[2] = c[2] - b[2];
  dir3[0] = b[0] - a[0];
  dir3[1] = b[1] - a[1];
  dir3[2] = b[2] - a[2];
  MathExtra::cross3(dir1, dir2, cross);
  double dsq = MathExtra::lensq3(cross)/MathExtra::lensq3(dir3);
  if(dsq > 1e-10) {
    return false;
  }

  return true;
}

/* ---------------------------------------------------------------------- */

void ExtrudeSurface::extrude_point_via_normal(double x[3], vtkIdType id, vtkDataArray *n, double scale)
{
  double normal[3];

  n->GetTuple(id, normal);
  for (vtkIdType i=0; i<3; ++i) {
    x[i] = x[i] + scale*normal[i];
  }
}

#endif
