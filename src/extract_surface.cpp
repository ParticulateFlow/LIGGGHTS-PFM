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
#include <stdlib.h>
#include <string.h>
#include "extract_surface.h"
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
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#if 0 // combination VTK 6.0.0 and c++11 gives compiler error
#include <vtkTriangle.h>
#endif
#include "error.h"

#define SMALL_AREA 5e-7

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ExtractSurface::ExtractSurface(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world, &me);
  binary = false;
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

    // check cell data
    if(!ugrid->GetCellData()->HasArray("cell_id")) {
      // create cell data
      vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
      cellData->SetName("cell_id");
      cellData->SetNumberOfComponents(1);
      for(int i = 0; i < ncells; ++i) {
        cellData->InsertNextValue(i);
      }
#if VTK_MAJOR_VERSION < 6
      ugrid->Update(); // force an update so we can set cell data
#endif
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
#if VTK_MAJOR_VERSION < 6
    polyData->Update(); // force an update so we can set cell data
#endif
    polyData->GetCellData()->AddArray(cellData);

    // triangulate surface
    triangulate(narg, arg, polyData);


    // add vectors for stress direction to grid using surface normals
    vtkSmartPointer<vtkPolyDataNormals> skinNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
    skinNormals->ComputeCellNormalsOn();
    skinNormals->ComputePointNormalsOn();
    skinNormals->SetFeatureAngle(1.0);
    skinNormals->SetInputConnection(surfaceFilter->GetOutputPort());

#if VTK_MAJOR_VERSION < 6
    vtkPolyData *input = skinNormals->GetOutput();
    input->Update();
#else
    skinNormals->Update();
    vtkPolyData *input = skinNormals->GetOutput();
#endif

    vtkIdType numCells = input->GetNumberOfCells();
    if(!ugrid->GetCellData()->HasArray("stress_ctrl_dir")) {
      vtkSmartPointer<vtkDoubleArray> normalsArray = vtkSmartPointer<vtkDoubleArray>::New();
      normalsArray->SetName("stress_ctrl_dir");
      normalsArray->SetNumberOfComponents(3);
      normalsArray->SetNumberOfTuples(ncells);
      double nullvec[3] = {0.0, 0.0, 0.0};
      for(int i = 0; i < ncells; ++i) {
        normalsArray->SetTuple(i, nullvec) ;
      }
#if VTK_MAJOR_VERSION < 6
      ugrid->Update();
#endif
      ugrid->GetCellData()->SetVectors(normalsArray);

      vtkSmartPointer<vtkIntArray> surfcellids = vtkIntArray::SafeDownCast(input->GetCellData()->GetScalars("cell_id"));
      vtkDataArray *normals = input->GetCellData()->GetNormals();
      for (int isurfcell=0; isurfcell < numCells; ++isurfcell) {
        double *snorm = normals->GetTuple(isurfcell);
        int cellid = surfcellids->GetVariantValue(isurfcell).ToInt();

        for (int ivolcell=0; ivolcell<ncells; ++ivolcell) {
          if(ugrid->GetCellData()->GetScalars("cell_id")->GetVariantValue(ivolcell).ToInt() == cellid) {
            double *stress_dir = normalsArray->GetTuple(ivolcell);

            if (snorm[0] < 0.) stress_dir[0] = 1.;
            else if (snorm[0] > 0.) stress_dir[0] = -1.;
            if (snorm[1] < 0.) stress_dir[1] = 1.;
            else if (snorm[1] > 0.) stress_dir[1] = -1.;
            if (snorm[2] < 0.) stress_dir[2] = 1.;
            else if (snorm[2] > 0.) stress_dir[2] = -1.;

            normalsArray->SetTuple(ivolcell, stress_dir);
            break;
          }
        }
      }

      vtkNew<vtkUnstructuredGridWriter> writer;
      if (binary) writer->SetFileTypeToBinary();
      else        writer->SetFileTypeToASCII();
  #if VTK_MAJOR_VERSION < 6
      writer->SetInput(ugrid);
  #else
      writer->SetInputData(ugrid);
  #endif
      const size_t copylen = strlen(filename)-strlen(suffix);
      char *newfilename = new char[copylen+21];
      memcpy(newfilename, filename, copylen);
      newfilename[copylen] = 0;
      strcat(newfilename, "_stress_ctrl_dir.vtk");
      writer->SetHeader("Generated by LIGGGHTS");
      writer->SetFileName(newfilename);
      writer->Write();
      delete newfilename;
    }

    // extrude faces to get regions for insertion
    if (narg >= 8) {
      extrude(narg, arg, input);
    }

    ugrid->Delete();
  }
}

/* ---------------------------------------------------------------------- */

void ExtractSurface::triangulate(int narg, char **arg, vtkDataSet* dset)
{
  vtkNew<vtkTriangleFilter> triangleFilter;
#if VTK_MAJOR_VERSION < 6
  triangleFilter->SetInput(dset);
#else
  triangleFilter->SetInputData(dset);
#endif
  triangleFilter->Update();

  // remove degenerate triangles
  vtkPolyData *pd = triangleFilter->GetOutput();
  vtkIdType numCells = pd->GetNumberOfCells();
  vtkIdType cellId;

  // mark cells for deletion
  for (cellId=0; cellId < numCells; ++cellId) {
#if 0 // combination VTK 6.0.0 and c++11 gives compiler error
    vtkTriangle *tri = dynamic_cast<vtkTriangle*>(pd->GetCell(cellId));
    if (tri->ComputeArea() < SMALL_AREA)
#else
    vtkPoints *pts = pd->GetCell(cellId)->GetPoints();
    double a[3], b[3], c[3], ab[3], ac[3], cross[3];
    pts->GetPoint(0, a);
    pts->GetPoint(1, c);
    pts->GetPoint(2, b);
    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];
    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];
    MathExtra::cross3(ab, ac, cross);
    double area = 0.5 * MathExtra::len3(cross);
    if (area < SMALL_AREA)
#endif
      pd->DeleteCell(cellId);
  }
  // actually remove cells marked for deletion
  pd->RemoveDeletedCells();

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

void ExtractSurface::extrude(int /*narg*/, char **arg, vtkDataSet* dset)
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

  vtkPolyData *input = vtkPolyData::SafeDownCast(dset);

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

bool ExtractSurface::collinear(double *a, double *b, double *c)
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

void ExtractSurface::extrude_point_via_normal(double x[3], vtkIdType id, vtkDataArray *n, double scale)
{
  double normal[3];

  n->GetTuple(id, normal);
  for (vtkIdType i=0; i<3; ++i) {
    x[i] = x[i] + scale*normal[i];
  }
}

#endif
