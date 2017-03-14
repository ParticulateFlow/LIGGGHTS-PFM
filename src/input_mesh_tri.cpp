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

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "memory.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include <math.h>
#include "vector_liggghts.h"
#include "input_mesh_tri.h"
#include "tri_mesh.h"

using namespace LAMMPS_NS;
enum{UNSUPPORTED_DATA_TYPE,INT,DOUBLE};

InputMeshTri::InputMeshTri(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
  verbose_(false),
  i_exclusion_list_(0),
  size_exclusion_list_(0),
  exclusion_list_(0),
  read_cell_data_(false),
  restart_(false)
{}

InputMeshTri::~InputMeshTri()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile(const char *filename, class TriMesh *mesh, bool verbose,
                               const int size_exclusion_list, int *exclusion_list,
                               bool read_cell_data, bool restart)
{
  verbose_ = verbose;
  size_exclusion_list_ = size_exclusion_list;
  exclusion_list_ = exclusion_list;
  read_cell_data_ = read_cell_data;
  restart_ = restart;

  if(strlen(filename) < 5)
    error->all(FLERR,"Illegal command, file name too short for input of triangular mesh");
  const char *ext = &(filename[strlen(filename)-3]);

  // read file
  // case STL file or VTK file

  bool is_stl = (strcmp(ext,"stl") == 0) || (strcmp(ext,"STL") == 0);
  bool is_vtk = (strcmp(ext,"vtk") == 0) || (strcmp(ext,"VTK") == 0);

  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set stl___file

  if (me == 0)
  {
    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[128];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(FLERR,str);
    }
  } else nonlammps_file = NULL;

  if(is_stl)
  {
      if (!read_cell_data_ && !restart_)
      {
          if (comm->me == 0) fprintf(screen,"\nReading STL file '%s' \n",filename);
          meshtrifile_stl(mesh);
      }
  }
  else if(is_vtk)
  {
      if (comm->me == 0 && screen) fprintf(screen,"\nReading VTK file '%s' \n",filename);
      meshtrifile_vtk(mesh);
  }
  else error->all(FLERR,"Illegal command, need either an STL file or a VTK file as input for triangular mesh.");

  if(nonlammps_file) fclose(nonlammps_file);
}

/* ----------------------------------------------------------------------
   process VTK file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_vtk(class TriMesh *mesh)
{
  int n,m;

  double **points = NULL;
  int ipoint = 0,npoints = 0;

  int **cells = NULL, *lines = NULL;
  int icell = 0,ncells = 0;

  int ntris = 0;
  int iLine = 0;

  int nLines = 0;
  int nPointLines = 0;

  bool parse_scalars = false;
  bool parse_field = false;
  bool parse_vectors = false;
  int dataType = UNSUPPORTED_DATA_TYPE;
  int nFieldArrays = 0;
  int itri = 0;
  char propertyname[256]={};
  const char communication[]="comm_exchange_borders";

  while (1)
  {
    // read a line from input script
    // n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line

    if (me == 0) {
      m = 0;
      while (1) {
        if (maxline-m < 2) reallocate(line,maxline,0);
        if (fgets(&line[m],maxline-m,nonlammps_file) == NULL) {
          if (m) n = strlen(line) + 1;
          else n = 0;
          break;
        }
        m = strlen(line);
        if (line[m-1] != '\n') continue;

        m--;
        while (m >= 0 && isspace(line[m])) m--;
        if (m < 0 || line[m] != '&') {
          line[m+1] = '\0';
          n = m+2;
          break;
        }
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      break;
    }

    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // lines start with 1 (not 0)
    nLines++;

    // parse one line from the file
    parse_nonlammps();
    // skip empty lines
    if(narg == 0){
         if (me == 0 && screen && verbose_)
            fprintf(screen,"Note: Skipping empty line in VTK mesh file\n");
      continue;
    }

    // Note that the first line of an ASCII VTK files start with '#' and will be skipped

    //increase line counter
    iLine++;

    if(iLine < 2) continue;

    if(iLine == 2)
    {
        if(strcmp(arg[0],"ASCII"))
            error->all(FLERR,"Expecting ASCII VTK mesh file, cannot continue");
        continue;
    }

    if(iLine == 3)
    {
        if(strcmp(arg[0],"DATASET") || strcmp(arg[1],"UNSTRUCTURED_GRID"))
            error->all(FLERR,"Expecting ASCII VTK unstructured grid mesh file, cannot continue");
        continue;
    }

    if(iLine == 4)
    {
        if(strcmp(arg[0],"POINTS"))
            error->all(FLERR,"Expecting 'POINTS' section in ASCII VTK mesh file, cannot continue");
        npoints = atoi(arg[1]);
        if(!read_cell_data_ && !restart_) memory->create<double>(points,npoints,3,"input_mesh:points");
        continue;
    }

    if(ipoint < npoints)
    {
        if(narg % 3)
            error->all(FLERR,"Expecting multiple of 3 values of point data in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        for(int i = 0; i < narg/3; ++i)
        {
            if(!read_cell_data_ && !restart_)
            {
                for(int j = 0; j < 3; ++j)
                {
                    points[ipoint][j] = atof(arg[3*i+j]);
                }
            }
            ++ipoint;
        }
        ++nPointLines;
        continue;
    }

    if(iLine == 5+nPointLines)
    {
        if(strcmp(arg[0],"CELLS"))
            error->all(FLERR,"Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);
        memory->create<int>(cells,ncells,3,"input_mesh:cells");
        if(!read_cell_data_ && !restart_) memory->create<int>(lines,ncells,"input_mesh:lines");
        continue;
    }

    //copy data of all which have 3 values - can be tri, polygon etc
    if(iLine <= 5+nPointLines+ncells)
    {
        if(narg == 4)
            for (int j=0;j<3;j++) cells[icell][j] = atoi(arg[1+j]);
        else
            cells[icell][0] = -1;

        if(!read_cell_data_ && !restart_) lines[icell] = nLines;

        icell++;
        continue;
    }

    if(iLine == 6+nPointLines+ncells)
    {
        if(strcmp(arg[0],"CELL_TYPES"))
            error->all(FLERR,"Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))
            error->all(FLERR,"Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        icell = 0;
        continue;
    }

    //only take triangles (cell type 5 according to VTK standard) - count them
    if(iLine <= 6+nPointLines+2*ncells)
    {
        if(strcmp(arg[0],"5")) cells[icell][0] = -1; //remove if not a tri
        else ntris++;
        icell++;
        continue;
    }

    // parsing cell data
    if(read_cell_data_)
    {
        if(iLine == 7+nPointLines+2*ncells)
        {
            if(strcmp(arg[0],"CELL_DATA"))
                error->all(FLERR,"Expecting 'CELL_DATA' section in ASCII VTK mesh file, cannot continue");
            if(ncells != atoi(arg[1]))
                error->all(FLERR,"Inconsistency in 'CELL_DATA' section in ASCII VTK mesh file, cannot continue");
            continue;
        }

        if(strcmp(arg[0],"SCALARS") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'SCALARS' section in ASCII VTK mesh file, cannot continue");
            if(narg > 3 && atoi(arg[3]) > 1)
                error->all(FLERR,"Support for 'SCALARS' with more than 1 component in ASCII VTK mesh file not implemented, cannot continue");
            strcpy(propertyname, arg[1]);
            if(strcmp(arg[2],"int") == 0 || strcmp(arg[2],"long") == 0)
            {
                mesh->prop().addElementProperty< ScalarContainer<int> >(propertyname,communication,"frame_invariant","restart_yes");
                dataType = INT;
            }
            else if(strcmp(arg[2],"float") == 0 || strcmp(arg[2],"double") == 0)
            {
                mesh->prop().addElementProperty< ScalarContainer<double> >(propertyname,communication,"frame_invariant","restart_yes");
                dataType = DOUBLE;
            }
            else
            {
                dataType = UNSUPPORTED_DATA_TYPE;
                error->warning(FLERR,"Data type for 'SCALARS' in ASCII VTK mesh file not implemented, skipping data");
            }
            parse_scalars = true;
            icell = 0;
            itri = 0;
            continue;
        }

        if(strcmp(arg[0],"LOOKUP_TABLE") == 0)
        {
            if(!parse_scalars)
                error->all(FLERR,"Unexpected 'LOOKUP_TABLE' section in ASCII VTK mesh file, cannot continue");
            if(strcmp(arg[1],"default"))
                error->all(FLERR,"Parsing of custom lookup table in ASCII VTK mesh file not implemented, cannot continue");
            continue;
        }

        if(strcmp(arg[0],"FIELD") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'FIELD' section in ASCII VTK mesh file, cannot continue");
            nFieldArrays = atoi(arg[2]);
            parse_field = true;
            icell = -1;
            continue;
        }

        if(strcmp(arg[0],"COLOR_SCALARS") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'COLOR_SCALARS' section in ASCII VTK mesh file, cannot continue");
            error->warning(FLERR,"Support for 'COLOR_SCALARS' section in ASCII VTK mesh file not implemented, skipping data");
            continue;
        }

        if(strcmp(arg[0],"VECTORS") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'VECTORS' section in ASCII VTK mesh file, cannot continue");
            strcpy(propertyname, arg[1]);
            if(strcmp(arg[2],"int") == 0 || strcmp(arg[2],"long") == 0)
            {
                mesh->prop().addElementProperty< VectorContainer<int,3> >(propertyname,communication,"frame_invariant","restart_yes");

                dataType = INT;
            }
            else if(strcmp(arg[2],"float") == 0 || strcmp(arg[2],"double") == 0)
            {
                mesh->prop().addElementProperty< VectorContainer<double,3> >(propertyname,communication,"frame_invariant","restart_yes");
                dataType = DOUBLE;
            }
            else
            {
                dataType = UNSUPPORTED_DATA_TYPE;
                error->warning(FLERR,"Data type for 'VECTORS' in ASCII VTK mesh file not implemented, skipping data");
            }
            parse_vectors = true;
            icell = 0;
            itri = 0;
            continue;
        }

        if(strcmp(arg[0],"NORMALS") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'NORMALS' section in ASCII VTK mesh file, cannot continue");
            error->warning(FLERR,"Support for 'NORMALS' section in ASCII VTK mesh file not implemented, skipping data");
            continue;
        }

        if(strcmp(arg[0],"TEXTURE_COORDINATES") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'TEXTURE_COORDINATES' section in ASCII VTK mesh file, cannot continue");
            error->warning(FLERR,"Support for 'TEXTURE_COORDINATES' section in ASCII VTK mesh file not implemented, skipping data");
            continue;
        }

        if(strcmp(arg[0],"TENSORS") == 0)
        {
            if(parse_scalars || parse_field || parse_vectors)
                error->all(FLERR,"Unexpected 'TENSORS' section in ASCII VTK mesh file, cannot continue");
            error->warning(FLERR,"Support for 'TENSORS' section in ASCII VTK mesh file not implemented, skipping data");
            continue;
        }

        if(parse_scalars)
        {
            if(restart_)
            {
                icell += narg;
            }
            else
            {
                switch(dataType)
                {
                case INT:
                    {
                        ScalarContainer<int> *ep = mesh->prop().getElementProperty<ScalarContainer<int> >(propertyname);
                        for(int i = 0; i < narg; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                ep->set(itri, atoi(arg[i]));
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case DOUBLE:
                    {
                        ScalarContainer<double> *ep = mesh->prop().getElementProperty<ScalarContainer<double> >(propertyname);
                        for(int i = 0; i < narg; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                ep->set(itri, atof(arg[i]));
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case UNSUPPORTED_DATA_TYPE:
                    icell += narg;
                    break;
                }
            }

            if(icell == ncells)
                parse_scalars = false;
            continue;
        }

        if(parse_field)
        {
            if(icell < 0)
            {
                strcpy(propertyname, arg[0]);
                if(atoi(arg[1]) > 1)
                    error->all(FLERR,"Support for field arrays with more than 1 component in ASCII VTK mesh file not implemented, cannot continue");
                if(ncells != atoi(arg[2]))
                    error->all(FLERR,"Inconsistency in 'FIELD' section in ASCII VTK mesh file, cannot continue");
                if(strcmp(arg[3],"int") == 0 || strcmp(arg[3],"long") == 0)
                {
                    mesh->prop().addElementProperty< ScalarContainer<int> >(propertyname,communication,"frame_invariant","restart_yes");
                    dataType = INT;
                }
                else if(strcmp(arg[3],"float") == 0 || strcmp(arg[3],"double") == 0)
                {
                    mesh->prop().addElementProperty< ScalarContainer<double> >(propertyname,communication,"frame_invariant","restart_yes");
                    dataType = DOUBLE;
                }
                else
                {
                    dataType = UNSUPPORTED_DATA_TYPE;
                    error->warning(FLERR,"Data type for 'FIELD' in ASCII VTK mesh file not implemented, skipping data");
                }
                icell = 0;
                itri = 0;
                continue;
            }

            if(restart_)
            {
                icell += narg;
            }
            else
            {
                switch(dataType)
                {
                case INT:
                    {
                        ScalarContainer<int> *ep = mesh->prop().getElementProperty<ScalarContainer<int> >(propertyname);
                        for(int i = 0; i < narg; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                ep->set(itri, atoi(arg[i]));
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case DOUBLE:
                    {
                        ScalarContainer<double> *ep = mesh->prop().getElementProperty<ScalarContainer<double> >(propertyname);
                        for(int i = 0; i < narg; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                ep->set(itri, atof(arg[i]));
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case UNSUPPORTED_DATA_TYPE:
                    icell += narg;
                    break;
                }
            }

            if(icell == ncells)
            {
                --nFieldArrays;
                if(nFieldArrays > 0)
                    icell = -1;
                else
                    parse_field = false;
            }
            continue;
        }

        if(parse_vectors)
        {
            if(narg % 3)
                error->all(FLERR,"Expecting multiple of 3 values of cell data in 'VECTORS' section of ASCII VTK mesh file, cannot continue");

            if(restart_)
            {
                icell += narg/3;
            }
            else
            {
                switch(dataType)
                {
                case INT:
                    {
                        VectorContainer<int,3> *ep = mesh->prop().getElementProperty<VectorContainer<int,3> >(propertyname);
                        for(int i = 0; i < narg/3; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                int vector[3];
                                vector[0] = atoi(arg[3*i + 0]);
                                vector[1] = atoi(arg[3*i + 1]);
                                vector[2] = atoi(arg[3*i + 2]);
                                ep->set(itri, vector);
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case DOUBLE:
                    {
                        VectorContainer<double,3> *ep = mesh->prop().getElementProperty<VectorContainer<double,3> >(propertyname);
                        for(int i = 0; i < narg/3; ++i)
                        {
                            if(cells[icell][0] != -1)
                            {
                                double vector[3];
                                vector[0] = atof(arg[3*i + 0]);
                                vector[1] = atof(arg[3*i + 1]);
                                vector[2] = atof(arg[3*i + 2]);
                                ep->set(itri, vector);
                                ++itri;
                            }
                            ++icell;
                        }
                        break;
                    }
                case UNSUPPORTED_DATA_TYPE:
                    icell += narg/3;
                    break;
                }
            }

            if(icell == ncells)
                parse_vectors = false;
            continue;
        }
    }
  }

  if(!read_cell_data_ && !restart_)
  {
      //now that everything is parsed, write the data into the mesh
      for(int i = 0; i < ncells; i++)
      {
          if(cells[i][0] == -1) continue;
          if(size_exclusion_list_ > 0 && lines[i] == exclusion_list_[i_exclusion_list_])
          {
             if(i_exclusion_list_ < size_exclusion_list_-1)
                i_exclusion_list_++;
             continue;
          }
          addTriangle(mesh,points[cells[i][0]],points[cells[i][1]],points[cells[i][2]],lines[i]);
      }
  }

  if(!read_cell_data_ && !restart_) memory->destroy<double>(points);
  if(!read_cell_data_ && !restart_) memory->destroy<int>(lines);
  memory->destroy<int>(cells);
}

/* ----------------------------------------------------------------------
   process STL file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile_stl(class TriMesh *mesh)
{
  int n,m;
  int iVertex = 0;
  double vertices[3][3];
  bool insideSolidObject = false;
  bool insideFacet = false;
  bool insideOuterLoop = false;

  int nLines = 0, nLinesTri = 0;

  while (1)
  {
    // read a line from input script
    // n = length of line including str terminator, 0 if end of file
    // if line ends in continuation char '&', concatenate next line

    if (me == 0) {
      m = 0;
      while (1) {
        if (maxline-m < 2) reallocate(line,maxline,0);
        if (fgets(&line[m],maxline-m,nonlammps_file) == NULL) {
          if (m) n = strlen(line) + 1;
          else n = 0;
          break;
        }
        m = strlen(line);
        if (line[m-1] != '\n') continue;

        m--;
        while (m >= 0 && isspace(line[m])) m--;
        if (m < 0 || line[m] != '&') {
          line[m+1] = '\0';
          n = m+2;
          break;
        }
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file

    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      break;
    }

    if (n > maxline) reallocate(line,maxline,n);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // lines start with 1 (not 0)
    nLines++;

    // parse one line from the stl file
    parse_nonlammps();

    // skip empty lines
    if(narg==0){
      if (me == 0 && screen && verbose_)
        fprintf(screen,"Note: Skipping empty line in STL file\n");
      continue;
    }

    // detect begin and end of a solid object, facet and vertices
    if (strcmp(arg[0],"solid") == 0)
    {
      if (insideSolidObject)
        error->all(FLERR,"Corrupt or unknown STL file: New solid object begins without closing prior solid object.");
      insideSolidObject=true;
      if (me == 0 && screen && verbose_){
        fprintf(screen,"Solid body detected in STL file\n");
      }
    }
    else if (strcmp(arg[0],"endsolid") == 0)
    {
      if (!insideSolidObject)
        error->all(FLERR,"Corrupt or unknown STL file: End of solid object found, but no begin.");
      insideSolidObject=false;
      if (me == 0 && screen && verbose_) {
        fprintf(screen,"End of solid body detected in STL file.\n");
      }
    }

    // detect begin and end of a facet within a solids object
    else if (strcmp(arg[0],"facet") == 0)
    {
      if (insideFacet)
        error->all(FLERR,"Corrupt or unknown STL file: New facet begins without closing prior facet.");
      if (!insideSolidObject)
        error->all(FLERR,"Corrupt or unknown STL file: New facet begins outside solid object.");
      insideFacet = true;

      nLinesTri = nLines;

      // check for keyword normal belonging to facet
      if (strcmp(arg[1],"normal") != 0)
        error->all(FLERR,"Corrupt or unknown STL file: Facet normal not defined.");

      // do not import facet normal (is calculated later)
    }
    else if (strcmp(arg[0],"endfacet") == 0)
    {
       if (!insideFacet)
         error->all(FLERR,"Corrupt or unknown STL file: End of facet found, but no begin.");
       insideFacet = false;
       if (iVertex != 3)
         error->all(FLERR,"Corrupt or unknown STL file: Number of vertices not equal to three (no triangle).");

      // add triangle to mesh
      //if (screen) printVec3D(screen,"vertex",vertices[0]);
      //if (screen) printVec3D(screen,"vertex",vertices[1]);
      //if (screen) printVec3D(screen,"vertex",vertices[2]);
      if(size_exclusion_list_ > 0 && nLinesTri == exclusion_list_[i_exclusion_list_])
      {
         if(i_exclusion_list_ < size_exclusion_list_-1)
            i_exclusion_list_++;
      }
      else
         addTriangle(mesh,vertices[0],vertices[1],vertices[2],nLinesTri);

      //if (me == 0 && screen) {
        //fprintf(screen,"  End of facet detected in in solid body.\n");
      //}
    }

    //detect begin and end of an outer loop within a facet
    else if (strcmp(arg[0],"outer") == 0)
    {
      if (insideOuterLoop)
        error->all(FLERR,"Corrupt or unknown STL file: New outer loop begins without closing prior outer loop.");
      if (!insideFacet)
        error->all(FLERR,"Corrupt or unknown STL file: New outer loop begins outside facet.");
      insideOuterLoop = true;
      iVertex = 0;

      //if (me == 0 && screen){
        //fprintf(screen,"    Outer loop detected in facet.\n");
      //}
    }
    else if (strcmp(arg[0],"endloop") == 0)
    {
      if (!insideOuterLoop)
        error->all(FLERR,"Corrupt or unknown STL file: End of outer loop found, but no begin.");
      insideOuterLoop=false;
      //if (me == 0 && screen) {
        //fprintf(screen,"    End of outer loop detected in facet.\n");
      //}
    }

    else if (strcmp(arg[0],"vertex") == 0)
    {
      if (!insideOuterLoop)
        error->all(FLERR,"Corrupt or unknown STL file: Vertex found outside a loop.");

      //if (me == 0 && screen) {
        //fprintf(screen,"      Vertex found.\n");
      //}

      // read the vertex
      for (int j=0;j<3;j++)
        vertices[iVertex][j]=atof(arg[1+j]);

      iVertex++;
      if (iVertex > 3)
        error->all(FLERR,"Corrupt or unknown STL file: Can not have more than 3 vertices "
                          "in a facet (only triangular meshes supported).");
    }
  }
}

/* ----------------------------------------------------------------------
   add a triangle to the mesh
------------------------------------------------------------------------- */

void InputMeshTri::addTriangle(TriMesh *mesh,double *a, double *b, double *c,int lineNumber)
{
    double **nodeTmp = create<double>(nodeTmp,3,3);
    for(int i=0;i<3;i++){
      nodeTmp[0][i] = a[i];
      nodeTmp[1][i] = b[i];
      nodeTmp[2][i] = c[i];
    }
    mesh->addElement(nodeTmp,lineNumber);
    destroy<double>(nodeTmp);
}
