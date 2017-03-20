/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "math.h"
#include "vector_liggghts.h"
#include "input_mesh_hex.h"
#include "region_mesh_hex.h"

using namespace LAMMPS_NS;
enum{UNSUPPORTED_DATA_TYPE,INT,DOUBLE};

#define MAXLINE 2048

InputMeshHex::InputMeshHex(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
  verbose_(false),
  read_cell_data_(false)
{}

InputMeshHex::~InputMeshHex()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshHex::meshhexfile(const char *filename, class RegHexMesh *mesh, bool verbose, bool read_cell_data)
{
  verbose_ = verbose;
  read_cell_data_ = read_cell_data;
  if(strlen(filename) < 5)
    error->all(FLERR,"Illegal command, file name too short for input of tet mesh");

  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set stl___file

  if (me == 0) {

    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[128];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(FLERR,str);
    }
  } else nonlammps_file = NULL;

  meshhexfile_vtk(mesh);

  if(nonlammps_file) fclose(nonlammps_file);

}

/* ----------------------------------------------------------------------
   process VTK file
------------------------------------------------------------------------- */

void InputMeshHex::meshhexfile_vtk(class RegHexMesh *mesh)
{
  int n,m;

  double phix = (mesh->rot_angle[0])*M_PI/180.;
  double phiy = (mesh->rot_angle[1])*M_PI/180.;
  double phiz = (mesh->rot_angle[2])*M_PI/180.;

  int flag_outside = 0;

  double **points = NULL;
  int ipoint = 0, npoints = 0;
  double vert_before_rot[3], vert_after_rot[3];

  int **cells = NULL;
  int icell = 0, ncells = 0;

  int nhexs = 0;
  int iLine = 0;

  int nLines = 0;
  int nPointLines = 0;

  bool parse_scalars = false;
  bool parse_field = false;
  bool parse_vectors = false;
  int dataType = UNSUPPORTED_DATA_TYPE;
  int nFieldArrays = 0;
  int ihex = 0;
  char propertyname[256]={};
  const char communication[]="comm_exchange_borders";
  int flag_other_than_hex = 0;

  while (1) {

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

    //parse one line from the file
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
        if(strcmp(arg[0],"ASCII")) error->all(FLERR,"Expecting ASCII VTK mesh file, cannot continue");
        continue;
    }

    if(iLine == 3)
    {
        if(strcmp(arg[0],"DATASET") || strcmp(arg[1],"UNSTRUCTURED_GRID")) error->all(FLERR,"Expecting ASCII VTK unstructured grid mesh file, cannot continue");
        continue;
    }

    if(iLine == 4)
    {
        if(strcmp(arg[0],"POINTS")) error->all(FLERR,"Expecting 'POINTS' section in ASCII VTK mesh file, cannot continue");
        npoints = atoi(arg[1]);
        memory->create(points,npoints,3,"input_mesh:points");
        continue;
    }

    if(ipoint < npoints)
    {
        if(narg % 3)
            error->all(FLERR,"Expecting multiple of 3 values of point data in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        //read the vertex, translate and scale it
        for(int i = 0; i < narg/3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                vert_before_rot[j] = (atof(arg[3*i+j])+(mesh->off_fact[j]))*(mesh->scale_fact);
            }

            //rotate the vertex
            vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
            vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
            vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

            if (!domain->is_in_domain(vert_after_rot))
                flag_outside = 1;

            //store the vertex

            vectorCopy3D(vert_after_rot,points[ipoint]);
            ++ipoint;
        }
        ++nPointLines;
        continue;
    }

    if(iLine == 5+nPointLines)
    {
        if(strcmp(arg[0],"CELLS")) error->all(FLERR,"Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);

        memory->create(cells,ncells,8,"input_mesh_hex:cells");
        continue;
    }

    //copy data of all which have 8 values - can be tet, quad, poly_line or triangle_strip
    if(iLine <= 5+nPointLines+ncells)
    {
        if(narg == 9)
            for (int j=0;j<8;j++) cells[icell][j] = atoi(arg[1+j]);
        else
            cells[icell][0] = -1;

        icell++;
        continue;
    }

    if(iLine == 6+nPointLines+ncells)
    {
        if(strcmp(arg[0],"CELL_TYPES")) error->all(FLERR,"Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))  error->all(FLERR,"Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        icell = 0;
        continue;
    }

    //only take hexahedra (cell type 12 according to VTK standard) - count them
    if(iLine <= 6+nPointLines+2*ncells)
    {
        if(strcmp(arg[0],"12"))
        {
            cells[icell][0] = -1; //remove if not a hex
            flag_other_than_hex = 1;
        }
        else nhexs++;
        icell++;
        continue;
    }

    if(strcmp(arg[0],"POINT_DATA") == 0)
    {
        error->warning(FLERR,"Support for 'POINT_DATA' in ASCII VTK mesh file not implemented, skipping data");
        break;
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
                mesh->prop().addElementProperty< ScalarContainer<int> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);
                dataType = INT;
            }
            else if(strcmp(arg[2],"float") == 0 || strcmp(arg[2],"double") == 0)
            {
                mesh->prop().addElementProperty< ScalarContainer<double> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);
                dataType = DOUBLE;
            }
            else
            {
                dataType = UNSUPPORTED_DATA_TYPE;
                error->warning(FLERR,"Data type for 'SCALARS' in ASCII VTK mesh file not implemented, skipping data");
            }
            parse_scalars = true;
            icell = 0;
            ihex = 0;
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
                mesh->prop().addElementProperty< VectorContainer<int,3> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);

                dataType = INT;
            }
            else if(strcmp(arg[2],"float") == 0 || strcmp(arg[2],"double") == 0)
            {
                mesh->prop().addElementProperty< VectorContainer<double,3> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);
                dataType = DOUBLE;
            }
            else
            {
                dataType = UNSUPPORTED_DATA_TYPE;
                error->warning(FLERR,"Data type for 'VECTORS' in ASCII VTK mesh file not implemented, skipping data");
            }
            parse_vectors = true;
            icell = 0;
            ihex = 0;
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
            switch(dataType)
            {
            case INT:
                {
                    ScalarContainer<int> *ep = mesh->prop().getElementProperty<ScalarContainer<int> >(propertyname);
                    for(int i = 0; i < narg; ++i)
                    {
                        if(cells[icell][0] != -1)
                        {
                            ep->set(ihex, atoi(arg[i]));
                            ++ihex;
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
                            ep->set(ihex, atof(arg[i]));
                            ++ihex;
                        }
                        ++icell;
                    }
                    break;
                }
            case UNSUPPORTED_DATA_TYPE:
                icell += narg;
                break;
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
                    mesh->prop().addElementProperty< ScalarContainer<int> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);
                    dataType = INT;
                }
                else if(strcmp(arg[3],"float") == 0 || strcmp(arg[3],"double") == 0)
                {
                    mesh->prop().addElementProperty< ScalarContainer<double> >(propertyname,communication,"frame_invariant","restart_yes",1,nhexs);
                    dataType = DOUBLE;
                }
                else
                {
                    dataType = UNSUPPORTED_DATA_TYPE;
                    error->warning(FLERR,"Data type for 'FIELD' in ASCII VTK mesh file not implemented, skipping data");
                }
                icell = 0;
                ihex = 0;
                continue;
            }

            switch(dataType)
            {
            case INT:
                {
                    ScalarContainer<int> *ep = mesh->prop().getElementProperty<ScalarContainer<int> >(propertyname);
                    for(int i = 0; i < narg; ++i)
                    {
                        if(cells[icell][0] != -1)
                        {
                            ep->set(ihex, atoi(arg[i]));
                            ++ihex;
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
                            ep->set(ihex, atof(arg[i]));
                            ++ihex;
                        }
                        ++icell;
                    }
                    break;
                }
            case UNSUPPORTED_DATA_TYPE:
                icell += narg;
                break;
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
                            ep->set(ihex, vector);
                            ++ihex;
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
                            ep->set(ihex, vector);
                            ++ihex;
                        }
                        ++icell;
                    }
                    break;
                }
            case UNSUPPORTED_DATA_TYPE:
                icell += narg/3;
                break;
            }

            if(icell == ncells)
                parse_vectors = false;
            continue;
        }
    }
  }

  //must throw an error here since regions may not extend outside box
  if(flag_outside)
    error->all(FLERR,"VTK mesh file is incompatible with simulation box: One or more vertices outside simulation box");

  //now that everything is parsed, write the data into the mesh
  while(mesh->nHexMax < nhexs) mesh->grow_arrays();

  double **hexnodes;
  memory->create(hexnodes,8,3,"input_mesh_hex:hexnodes");

  for(int i = 0; i < ncells; i++)
  {
      if(cells[i][0] == -1) continue;
      for (int j = 0; j < 8; j++)
        vectorCopy3D(points[cells[i][j]],hexnodes[j]);

      mesh->add_hex(hexnodes);
  }

  if(flag_other_than_hex && comm->me == 0)
    error->warning(FLERR,"VTK file contains other types than hexahedra - only hexs are currently supported in LIGGGHTS, other cells are discarded");

  if(ncells == 0)
    error->all(FLERR,"VTK mesh file containing no hex cell - cannnot continue");

  memory->destroy(hexnodes);
  memory->destroy(points);
  memory->destroy(cells);
}

#endif
