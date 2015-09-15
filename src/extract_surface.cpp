/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2015- JKU Linz

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
#include "extract_surface.h"
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
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ExtractSurface::ExtractSurface(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world, &me);
}

/* ---------------------------------------------------------------------- */

ExtractSurface::~ExtractSurface()
{
}

/* ---------------------------------------------------------------------- */

template<class TReader> vtkDataSet* ExtractSurface::read_file(const char*filename)
{
  vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
  reader->SetFileName(filename);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

/* ---------------------------------------------------------------------- */

void ExtractSurface::command(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Illegal extract_surface command");

  if (me == 0) {
    char *filename = arg[0];
    char *suffix = filename + strlen(filename) - strlen(".vtk");
    vtkDataSet* ugrid = NULL;

    if (suffix > filename && strcmp(suffix,".vtk") == 0) {
      vtkDataReader *chooser = vtkDataReader::New();
      chooser->SetFileName(filename);
      if (!chooser->IsFileUnstructuredGrid()) {
        chooser->Delete();
        error->all(FLERR,"extract_surface requires unstructured grid dataset");
      }
      ugrid = read_file<vtkUnstructuredGridReader>(filename);
    } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
      vtkXMLUnstructuredGridReader *chooser = vtkXMLUnstructuredGridReader::New();
      if (!chooser->CanReadFile(filename)) {
        chooser->Delete();
        error->all(FLERR,"extract_surface cannot read input file");
      }
      ugrid = read_file<vtkXMLUnstructuredGridReader>(filename);
    } else {
      error->all(FLERR,"extract_surface: invalid input file");
    }

    // check cell types
    int ncells = ugrid->GetNumberOfCells();
    for (int i = 0; i < ncells; ++i) {
      if(ugrid->GetCellType(i) != VTK_HEXAHEDRON)
        error->all(FLERR,"extract_surface: input file contains non-hexahedral cells");
    }

    if(!ugrid->GetCellData()->HasArray("cell_id")) {
      // create cell data
      vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
      cellData->SetName("cell_id");
      cellData->SetNumberOfComponents(1);
      for(int i = 0; i < ncells; ++i) {
        cellData->InsertNextValue(i);
      }
      ugrid->Update(); // force an update so we can set cell data
      ugrid->GetCellData()->SetScalars(cellData);
    }

    // extract surface
    vtkNew<vtkDataSetSurfaceFilter> surfaceFilter;
#if VTK_MAJOR_VERSION < 6
    surfaceFilter->SetInput(ugrid);
#else
    surfaceFilter->SetInputData(ugrid);
#endif
    surfaceFilter->Update();

    vtkPolyData* polyData = surfaceFilter->GetOutput();
    ncells = polyData->GetNumberOfCells();
    // create cell data
    vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
    cellData->SetName("face_id");
    cellData->SetNumberOfComponents(1);
    for (int i = 0; i < ncells; ++i) {
      cellData->InsertNextValue(i);
    }
    polyData->Update(); // force an update so we can set cell data
    polyData->GetCellData()->AddArray(cellData);

    // triangulate surface
    vtkNew<vtkTriangleFilter> triangleFilter;
    triangleFilter->SetInputConnection(surfaceFilter->GetOutputPort());
    triangleFilter->Update();

    vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
#if VTK_MAJOR_VERSION < 6
    appendFilter->AddInput(triangleFilter->GetOutput());
#else
    appendFilter->AddInputData(triangleFilter->GetOutput());
#endif
    appendFilter->Update();

    // write surface data to file
    bool binary = false;
    if (narg > 3 && strcmp(arg[narg-1],"binary") == 0) binary = true;

    filename = arg[2];
    suffix = filename + strlen(filename) - strlen(".vtu");

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

    // extrude faces to get regions for insertion
    if (narg > 6) {
      filename = arg[3];
      if(strcmp(arg[4],"extrude_length"))
        error->all(FLERR,"extract_surface: expecting keyword 'extrude_length'");
      double ScaleFactor = atof(arg[5]);
      if(strcmp(arg[6],"min_rad"))
        error->all(FLERR,"extract_surface: expecting keyword 'min_rad'");
      double atom_radius = atof(arg[7]);

      const double scale_out = atom_radius;
      const double scale_in = atom_radius-ScaleFactor;

      vtkSmartPointer<vtkPolyDataNormals> skinNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
      skinNormals->SetInputConnection(surfaceFilter->GetOutputPort());
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

      for (inCellId=0; inCellId < numCells; ++inCellId) {
        mesh->GetCell(inCellId,cell);
        cellPts = cell->GetPointIds();

        if (cell->GetCellDimension() == 2) {
          if (cell->GetNumberOfPoints() != 4)
            error->all(FLERR,"Unexpected number of points in vtk cell\n");

          newHexahedra->InsertNextCell(8);

          for (int i = 0; i < 4; ++i)
            newHexahedra->InsertCellPoint(cellPts->GetId(i));
          for (int i = 0; i < 4; ++i)
            newHexahedra->InsertCellPoint(cellPts->GetId(i)+numPts);
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

    ugrid->Delete();
  }
}

/* ---------------------------------------------------------------------- */

void ExtractSurface::extrude_point_via_normal(double x[3], vtkIdType id, vtkDataArray *n, double scale)
{
  double normal[3];

  n->GetTuple(id, normal);
  for (vtkIdType i=0; i<3; ++i) {
    x[i] = x[i] + scale*normal[i];
  }
}

#endif
