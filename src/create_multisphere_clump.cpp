/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2019- JKU Linz

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

#if defined(LAMMPS_VTK) // do not use #ifdef here (VS C++ bug)
#include "lmptype.h"
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include "create_multisphere_clump.h"
#include "math_extra.h"
#include "random_park.h"
#include "domain.h"
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSet.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSTLReader.h>
#include <vtkPLYReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkButterflySubdivisionFilter.h>
#if 0 // combination VTK 6.0.0 and c++11 gives compiler error
#include <vtkTriangle.h>
#endif
#include "error.h"

#define SMALL_DIST_SQ 1e-12

using namespace LAMMPS_NS;

int CreateMultisphereClump::tag = 0;
/* ---------------------------------------------------------------------- */

CreateMultisphereClump::CreateMultisphereClump(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world, &me);
  atom_type = 1;
  sign_normals = -1.0;
  binary = false;
  absolute_dmin = true;
  density = 1000.0;
  iarg = 0;
  seed = 1;
  random = new RanPark(lmp,seed);
}

/* ---------------------------------------------------------------------- */

CreateMultisphereClump::~CreateMultisphereClump()
{
}

/* ---------------------------------------------------------------------- */

template<class TReader> vtkDataSet* CreateMultisphereClump::read_file(const char*filename)
{
  vtkSmartPointer<TReader> reader = vtkSmartPointer<TReader>::New();
  reader->SetFileName(filename);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::command(int narg, char **arg)
{
  if (narg < 12) error->all(FLERR,"Illegal create_multisphere_clump command");

  iarg = 0;
  radii.clear();
  sx.clear();
  sy.clear();
  sz.clear();

  if (me == 0) {
    if (strcmp(arg[iarg++],"dmin") != 0)
      error->one(FLERR,"create_multisphere_clump command expects keyword 'dmin'");

    if(strcmp(arg[iarg],"radius_ratio") == 0)
        absolute_dmin = false;
    else if(strcmp(arg[iarg],"absolute") == 0)
        absolute_dmin = true;
    else
        error->one(FLERR,"expecting keyword 'absolute' or 'radius_ratio'");

    ++iarg;

    dmin = atof(arg[iarg++]);

    if (strcmp(arg[iarg++],"rmin") != 0)
      error->one(FLERR,"create_multisphere_clump command expects keyword 'rmin'");
    rmin = atof(arg[iarg++]);

    if (strcmp(arg[iarg++],"pmax") != 0)
      error->one(FLERR,"create_multisphere_clump command expects keyword 'pmax'");
    pmax = atof(arg[iarg++]);

    if (strcmp(arg[iarg++],"seed") != 0)
      error->one(FLERR,"create_multisphere_clump command expects keyword 'seed'");
    seed = atoi(arg[iarg++]);
    random->reset(seed);

    if (strcmp(arg[iarg++],"surfacefile") != 0)
      error->one(FLERR,"create_multisphere_clump command expects keyword 'surfacefile'");

    char *filename = arg[iarg++];
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
        error->one(FLERR,"create_multisphere_clump requires unstructured grid dataset or polydata");
      }
    } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
      vtkXMLUnstructuredGridReader *chooser = vtkXMLUnstructuredGridReader::New();
      if (!chooser->CanReadFile(filename)) {
        chooser->Delete();
        error->one(FLERR,"create_multisphere_clump cannot read input file");
      }
      dset = read_file<vtkXMLUnstructuredGridReader>(filename);
    } else if (suffix > filename && strcmp(suffix,".ply") == 0) {
      vtkPLYReader *chooser = vtkPLYReader::New();
      if (!chooser->CanReadFile(filename)) {
        chooser->Delete();
        error->one(FLERR,"create_multisphere_clump cannot read PLY input file");
      }
      dset = read_file<vtkPLYReader>(filename);
    } else if (suffix > filename && strcmp(suffix,".stl") == 0) {
      dset = read_file<vtkSTLReader>(filename);
    } else {
      error->one(FLERR,"create_multisphere_clump: invalid input file");
    }

    // check cell types
    int ncells = dset->GetNumberOfCells();
    for (int i = 0; i < ncells; ++i) {
      if (dset->GetCell(i)->GetCellDimension() != 2)
        error->one(FLERR,"create_multisphere_clump: input file contains unsupported cell type with dim != 2");
    }

    vtkNew<vtkGeometryFilter> geometryFilter;
#if VTK_MAJOR_VERSION < 6
    geometryFilter->SetInput(dset);
#else
    geometryFilter->SetInputData(dset);
#endif
    geometryFilter->Update();
    vtkSmartPointer<vtkPolyData> pd = geometryFilter->GetOutput();

    // triangulate surface
    pd = triangulate(pd);

    if (strcmp(arg[iarg],"invert_normals") == 0) {
      ++iarg;
      if (strcmp(arg[iarg++],"yes") == 0)
        sign_normals = -1.0;
      else
        sign_normals = 1.0;
    }

    // subdivide faces to get more points for better fit
    if (strcmp(arg[iarg],"subdivide") == 0) {
        ++iarg;
        pd = subdivide(narg, arg, pd);
        write_subdiv_file(narg, arg, pd);
    }

    generate_spheres(pd);

    if (narg <= iarg+1)
      error->one(FLERR, "create_multisphere_clump: not enough arguments");

    if (strcmp(arg[iarg],"clumpfile") == 0) {
      ++iarg;
      write_clump_file(arg[iarg++]);
      if (iarg < narg) {
        write_clump_file_debug(arg[iarg]);
      }
    } else if (strcmp(arg[iarg],"datafile") == 0) {
      ++iarg;
      int farg = iarg;
      atom_type = 1;
      density = 1000.0;
      ++iarg;
      if (iarg < narg) {
        atom_type = atoi(arg[iarg++]);
      }
      if (iarg < narg) {
        density = atof(arg[iarg++]);
      }
      write_data_file(arg[farg]);
    } else {
      error->one(FLERR,"create_multisphere_clump command expects keyword 'clumpfile' or 'datafile'");
    }

    dset->Delete();
  }
}

/* ---------------------------------------------------------------------- */

vtkSmartPointer<vtkPolyData> CreateMultisphereClump::triangulate(vtkPolyData* dset)
{
  vtkNew<vtkTriangleFilter> triangleFilter;
#if VTK_MAJOR_VERSION < 6
  triangleFilter->SetInput(dset);
#else
  triangleFilter->SetInputData(dset);
#endif
  triangleFilter->Update();

  vtkSmartPointer<vtkPolyData> triPD = triangleFilter->GetOutput();

  return triPD;
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::write_data_file(const char* filename)
{
  if (me == 0) {
    FILE* datafile = fopen(filename,"wt");
    if (datafile) {
      fprintf(datafile,"# LIGGGHTS multi-sphere clump data\n"
                        "%lu atoms\n%d atom types\n", radii.size(),(atom_type>atom->ntypes)?atom_type:atom->ntypes);
      fprintf(datafile,"%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nAtoms\n\n",
              domain->boxlo[0], domain->boxhi[0],
              domain->boxlo[1], domain->boxhi[1],
              domain->boxlo[2], domain->boxhi[2]);
      for (unsigned int isphere=0; isphere<radii.size(); ++isphere,++tag) {
        // atom-ID atom-type diameter density x y z
        fprintf(datafile,"%u %d %f %f %f %f %f\n",
                atom->tag_max()+1+tag, atom_type, 2.*radii[isphere], density,
                sx[isphere], sy[isphere], sz[isphere]);
      }
      fclose(datafile);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::write_clump_file(const char* filename)
{
  if (me == 0) {
    FILE* clumpfile = fopen(filename,"wt");
    if (clumpfile) {
      fprintf(clumpfile,"# LIGGGHTS multi-sphere clump data\n# %lu atoms\n# x y z radius\n", radii.size());
      for (unsigned int isphere=0; isphere<radii.size(); ++isphere) {
        fprintf(clumpfile,"%f %f %f %f\n", sx[isphere], sy[isphere], sz[isphere], radii[isphere]); // x y z r
      }
      fclose(clumpfile);
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::write_clump_file_debug(const char* filename)
{
  if (me == 0) {
    FILE* clumpfile = fopen(filename,"wt");
    if (clumpfile) {
      const unsigned int nspheres = radii.size();
      fprintf(clumpfile,"# vtk DataFile Version 3.0\n");
      fprintf(clumpfile,"LIGGGHTS multisphere/clump\n");
      fprintf(clumpfile,"ASCII\n");
      fprintf(clumpfile,"DATASET POLYDATA\n");
      fprintf(clumpfile,"POINTS %d float\n", nspheres);
      for (unsigned int isphere=0; isphere<nspheres; ++isphere) {
        fprintf(clumpfile,"%f %f %f\n", sx[isphere], sy[isphere], sz[isphere]); // x y z
      }
      fprintf(clumpfile,"VERTICES %u %u\n", nspheres, 2*nspheres);
      for (unsigned int isphere=0; isphere<nspheres; ++isphere) {
        fprintf(clumpfile,"1 %d\n", isphere);
      }
      fprintf(clumpfile,"POINT_DATA %u\n", nspheres);
      fprintf(clumpfile,"SCALARS radius float 1\n");
      fprintf(clumpfile,"LOOKUP_TABLE default\n");
      for (unsigned int isphere=0; isphere<nspheres; ++isphere) {
        fprintf(clumpfile,"%f\n", radii[isphere]); // r
      }
      fclose(clumpfile);
    }
  }
}

/* ---------------------------------------------------------------------- */

vtkSmartPointer<vtkPolyData> CreateMultisphereClump::subdivide(int narg, char **arg, vtkPolyData* dset)
{
  int numberOfSubdivisions = 2;
  vtkPolyData *originalMesh = dset;

  vtkSmartPointer<vtkPolyDataAlgorithm> subdivisionFilter;

  int sdtype = -1; // linear=0, loop=1, butterfly=2

  if (strcmp(arg[iarg],"linear") == 0) {
    sdtype = 0;
  } else if (strcmp(arg[iarg],"loop") == 0) {
    sdtype = 1;
  } else if (strcmp(arg[iarg],"butterfly") == 0) {
    sdtype = 2;
  } else {
    error->one(FLERR, "Unknown subdivision type");
  }
  ++iarg;

  if (narg <= iarg)
    error->one(FLERR, "create_multisphere_clump: not enough arguments");

  numberOfSubdivisions = atoi(arg[iarg++]);

  switch (sdtype)
  {
  case 0:
    subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
    dynamic_cast<vtkLinearSubdivisionFilter*>(subdivisionFilter.GetPointer())
        ->SetNumberOfSubdivisions(numberOfSubdivisions);
    break;
  case 1:
    subdivisionFilter = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
    dynamic_cast<vtkLoopSubdivisionFilter*>(subdivisionFilter.GetPointer())
        ->SetNumberOfSubdivisions(numberOfSubdivisions);
    break;
  case 2:
    subdivisionFilter = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
    dynamic_cast<vtkButterflySubdivisionFilter*>(subdivisionFilter.GetPointer())
        ->SetNumberOfSubdivisions(numberOfSubdivisions);
    break;
  default:
    error->one(FLERR, "Unknown subdivision type");
    break;
  }

#if VTK_MAJOR_VERSION < 6
  subdivisionFilter->SetInput(originalMesh);
#else
  subdivisionFilter->SetInputData(originalMesh);
#endif
  subdivisionFilter->Update();

  vtkSmartPointer<vtkPolyData> subdivPD = subdivisionFilter->GetOutput();

  return subdivPD;
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::write_subdiv_file(int narg, char **arg, vtkPolyData* pd)
{
  if (narg <= iarg)
    error->one(FLERR, "create_multisphere_clump: not enough arguments");

  vtkPolyData *subdivMesh = pd;

  if (strcmp(arg[iarg],"subdivisionfile") == 0) {
      ++iarg;
      if (narg <= iarg)
        error->one(FLERR, "create_multisphere_clump: not enough arguments");
      char *filename = arg[iarg++];
      char *suffix = filename + strlen(filename) - strlen(".vtp");

      if (suffix > filename && strcmp(suffix,".vtp") == 0) { // xml format
        vtkNew<vtkXMLPolyDataWriter> writer;
        if (binary) writer->SetDataModeToBinary();
        else        writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION < 6
        writer->SetInput(subdivMesh);
#else
        writer->SetInputData(subdivMesh);
#endif
        writer->SetFileName(filename);
        writer->Write();
      } else { // legacy format
        vtkNew<vtkPolyDataWriter> writer;
        if (binary) writer->SetFileTypeToBinary();
        else        writer->SetFileTypeToASCII();
#if VTK_MAJOR_VERSION < 6
        writer->SetInput(subdivMesh);
#else
        writer->SetInputData(subdivMesh);
#endif
        writer->SetHeader("Generated by LIGGGHTS");
        writer->SetFileName(filename);
        writer->Write();
    }
  }
}

/* ---------------------------------------------------------------------- */

void CreateMultisphereClump::generate_spheres(vtkPolyData* dset)
{
  vtkSmartPointer<vtkPolyDataNormals> skinNormals = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION < 6
  skinNormals->SetInput(dset);
#else
  skinNormals->SetInputData(dset);
#endif
  skinNormals->SplittingOff();

#if VTK_MAJOR_VERSION < 6
  vtkPolyData *input = skinNormals->GetOutput();
  input->Update();
#else
  skinNormals->Update();
  vtkPolyData *input = skinNormals->GetOutput();
#endif

  vtkIdType numPts = input->GetNumberOfPoints();

  if (numPts < 1) {
    error->one(FLERR,"No data to generate clump!");
  }

  vtkPointData *pd = input->GetPointData();
  vtkDataArray *inNormals = pd->GetNormals();
  vtkPoints *inPts = input->GetPoints();
  vtkIdType ptId, ptId1, tmpId;

  double pd_bounds[6];
  input->GetBounds(pd_bounds);
  double sizex = pd_bounds[1] - pd_bounds[0];
  double sizey = pd_bounds[3] - pd_bounds[2];
  double sizez = pd_bounds[5] - pd_bounds[4];
  double r_bound = 0.5*sqrt(sizex*sizex + sizey*sizey + sizez*sizez);
  double px[3], px1[3];
  double normal[3];
  double sxyz[3];
  double srad;

  std::set<vtkIdType> eligiblePts;
  for (ptId=0; ptId < numPts; ++ptId) {
    eligiblePts.insert(ptId);
  }

  for (ptId=0; ptId < numPts; ++ptId) {

    if (eligiblePts.find(ptId) == eligiblePts.end())
      continue;

    if (random->uniform() > pmax)
      continue;

    inPts->GetPoint(ptId, px);
    if (!check_sphere_distance(px))
      continue;

    inNormals->GetTuple(ptId, normal);

    srad = r_bound; // start with max radius a sphere can have in this clump
    MathExtra::addscaled3(px, normal, sign_normals*srad, sxyz); // initial center of sphere

    tmpId = -1;

    // shrink sphere until there are not more points inside
    for (ptId1=0; ptId1 < numPts; ++ptId1) {
      // exclude self
      if (ptId1 == ptId)
        continue;

      inPts->GetPoint(ptId1, px1);

      if (is_point_in_sphere(sxyz, srad, px1)) {
        double dir1[3], cross[3];
        MathExtra::sub3(px1, px, dir1);
        MathExtra::cross3(dir1, normal, cross);
        double shortestdistsq = MathExtra::lensq3(cross); // sq. dist. between line along px-normal and px1
        double distsq = MathExtra::lensq3(dir1); // sq. dist. between px and px1

        double lambdasq = (distsq - shortestdistsq);
        if (lambdasq > SMALL_DIST_SQ) {
          srad = 0.5 * distsq / sqrt(lambdasq); // update sphere radius
          MathExtra::addscaled3(px, normal, sign_normals*srad, sxyz); // update sphere center
          tmpId = ptId1;
        }
      }
    }

    if (tmpId >= 0 && srad >= rmin && is_sphere_in_bounds(pd_bounds, sxyz, srad)) {
      // save sphere data for later
      radii.push_back(srad);
      sx.push_back(sxyz[0]);
      sy.push_back(sxyz[1]);
      sz.push_back(sxyz[2]);
      eligiblePts.erase(tmpId); // don't add another sphere at this sphere's second point
    }
  }
}

/* ---------------------------------------------------------------------- */

bool CreateMultisphereClump::check_sphere_distance(const double *x)
{
  double center[3];
  double radius;

  for (unsigned int isphere=0; isphere<radii.size(); ++isphere) {
    center[0] = sx[isphere];
    center[1] = sy[isphere];
    center[2] = sz[isphere];
    radius = absolute_dmin?(radii[isphere]+dmin):(radii[isphere]*(1.+dmin));

    if (is_point_in_sphere(center, radius, x))
      return false;
  }

  return true;
}

/* ---------------------------------------------------------------------- */

bool CreateMultisphereClump::is_point_in_sphere(const double *center,
                                                double radius,
                                                const double *x)
{
  double delta[3];
  MathExtra::sub3(center, x, delta);

  return (MathExtra::lensq3(delta) < radius*radius);
}

/* ---------------------------------------------------------------------- */

bool CreateMultisphereClump::is_sphere_in_bounds(const double* bounds,
                                                 const double *center,
                                                 double radius)
{
  if (center[0]-radius >= bounds[0] && center[0]+radius <= bounds[1] &&
      center[1]-radius >= bounds[2] && center[1]+radius <= bounds[3] &&
      center[2]-radius >= bounds[4] && center[2]+radius <= bounds[5])
    return true;

  return false;
}

#endif
