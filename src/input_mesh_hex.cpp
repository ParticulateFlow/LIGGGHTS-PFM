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

#define MAXLINE 2048

InputMeshHex::InputMeshHex(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv),
verbose_(false)
{}

InputMeshHex::~InputMeshHex()
{}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshHex::meshhexfile(const char *filename, class RegHexMesh *mesh, bool verbose)
{
  verbose_ = verbose;

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
         if (me == 0 && verbose_)
            fprintf(screen,"Note: Skipping empty line in VTK mesh file\n");
      continue;
    }

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

    if(iLine <= 4+npoints)
    {
        if(narg != 3) error->all(FLERR,"Expecting 3 values for each point in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        //read the vertex, translate and scale it
        for (int j=0;j<3;j++) vert_before_rot[j]=(atof(arg[j])+(mesh->off_fact[j]))*(mesh->scale_fact);

        //rotate the vertex
        vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
        vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
        vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

        if (!domain->is_in_domain(vert_after_rot))
            flag_outside = 1;

        //store the vertex

        vectorCopy3D(vert_after_rot,points[ipoint]);
        ipoint++;
        continue;
    }

    if(iLine == 5+npoints)
    {
        if(strcmp(arg[0],"CELLS")) error->all(FLERR,"Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);

        memory->create(cells,ncells,8,"input_mesh_hex:cells");
        continue;
    }

    //copy data of all which have 8 values - can be tet, quad, poly_line or triangle_strip
    if(iLine <= 5+npoints+ncells)
    {
        if(narg == 9)
        {
            for (int j=0;j<8;j++) cells[icell][j] = atoi(arg[1+j]);

        }
        else
        {
            cells[icell][0] = -1;
        }

        icell++;
        continue;
    }

    if(iLine == 6+npoints+ncells)
    {
        if(strcmp(arg[0],"CELL_TYPES")) error->all(FLERR,"Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))  error->all(FLERR,"Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        icell = 0;

        continue;
    }

    //only take hexahedra (cell type 12 according to VTK standard) - count them
    if(iLine <= 6+npoints+2*ncells)
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
