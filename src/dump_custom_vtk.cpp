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
   Daniel Queteschiner, daniel.queteschiner@dcs-computing.com
   Richard Berger <richard.berger@jku.at>
------------------------------------------------------------------------- */

#ifdef LAMMPS_VTK
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "dump_custom_vtk.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "group.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include <vector>
#include <sstream>
#include <vtkVersion.h>
#ifndef VTK_MAJOR_VERSION
#include <vtkConfigure.h>
#endif
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h>

using namespace LAMMPS_NS;

// customize by
// * adding an enum constant (add vector components in consecutive order)
// * adding a pack_*(int) function for the value
// * adjusting parse_fields function to add the pack_* function to pack_choice
//   (in case of vectors, adjust identify_vectors as well)
// * adjusting thresh part in modify_param and count functions

enum{X,Y,Z, // required for vtk, must come first
     ID,MOL,TYPE,ELEMENT,MASS,
     XS,YS,ZS,XSTRI,YSTRI,ZSTRI,XU,YU,ZU,XUTRI,YUTRI,ZUTRI,
     XSU,YSU,ZSU,XSUTRI,YSUTRI,ZSUTRI,
     IX,IY,IZ,
     VX,VY,VZ,FX,FY,FZ,
     Q, MUX,MUY,MUZ,MU,RADIUS,DIAMETER,
     OMEGAX,OMEGAY,OMEGAZ,ANGMOMX,ANGMOMY,ANGMOMZ,
     TQX,TQY,TQZ,SPIN,ERADIUS,ERVEL,ERFORCE,
     DENSITY, RHO, P,
// superquadric start
     SHAPEX, SHAPEY, SHAPEZ,
     QUAT1, QUAT2, QUAT3, QUAT4,
     BLOCKINESS1, BLOCKINESS2,
     INERTIAX, INERTIAY, INERTIAZ,
// superquadric end
     THREAD, //NP modified C.K. included DENSITY .. A.A. included RHO and P .. R.B. included THREAD
     VARIABLE,COMPUTE,FIX,INAME,DNAME,
     ATTRIBUTES}; // must come last
enum{LT,LE,GT,GE,EQ,NEQ};
enum{INT,DOUBLE,STRING};    // same as in DumpCFG
enum{VTK,VTP,VTU,PVTP,PVTU}; // file formats

/* ---------------------------------------------------------------------- */

DumpCustomVTK::DumpCustomVTK(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg == 5) error->all(FLERR,"No dump custom/vtk arguments specified");

  clearstep = 1;

  nevery = force->inumeric(FLERR,arg[3]);

  // size_one may be shrunk below if additional optional args exist

  size_one = nfield = narg - 5;
  pack_choice.clear();
  vtype.clear();
  name.clear();

  iregion = -1;
  idregion = NULL;
  nthresh = 0;
  thresh_array = NULL;
  thresh_op = NULL;
  thresh_value = NULL;

  // computes, fixes, variables which the dump accesses

  ncompute = 0;
  id_compute = NULL;
  compute = NULL;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;

  nvariable = 0;
  id_variable = NULL;
  variable = NULL;
  vbuf = NULL;

  ncustom = 0;
  id_custom = NULL;
  flag_custom = NULL;

  myarrays.clear();
  n_calls_ = 0;

  // process attributes
  // ioptional = start of additional optional args
  // only dump image style processes optional args

  ioptional = parse_fields(narg,arg);

  if (ioptional < narg)
    error->all(FLERR,"Invalid attribute in dump custom/vtk command");
  size_one = pack_choice.size();
  current_pack_choice_key = -1;

  if (filewriter) reset_vtk_data_containers();

  // atom selection arrays

  maxlocal = 0;
  choose = NULL;
  dchoose = NULL;
  clist = NULL;

  // element names

  ntypes = atom->ntypes;
  typenames = NULL;

  label = NULL; //NP modified C.K.

  write_domain = 1;

  {
    // parallel vtp/vtu requires proc number to be preceded by underscore '_'
    multiname_ex = NULL;
    char *ptr = strchr(filename,'%');
    if (ptr) {
      multiname_ex = new char[strlen(filename) + 16];
      *ptr = '\0';
      sprintf(multiname_ex,"%s_%d%s",filename,me,ptr+1);
      *ptr = '%';
    }
  }

  vtk_file_format = VTK;

  char *suffix = filename + strlen(filename) - strlen(".vtp");
  if (suffix > filename && strcmp(suffix,".vtp") == 0) {
    if (multiproc) vtk_file_format = PVTP;
    else           vtk_file_format = VTP;
  } else if (suffix > filename && strcmp(suffix,".vtu") == 0) {
    if (multiproc) vtk_file_format = PVTU;
    else           vtk_file_format = VTU;
  }

  if (vtk_file_format == VTK) { // no multiproc support for legacy vtk format
    if (me != 0) filewriter = 0;
    fileproc = 0;
    multiproc = 0;
    nclusterprocs = nprocs;
  }

  filecurrent = NULL;
  domainfilecurrent = NULL;
  parallelfilecurrent = NULL;
}

/* ---------------------------------------------------------------------- */

DumpCustomVTK::~DumpCustomVTK()
{
  delete [] filecurrent;
  delete [] domainfilecurrent;
  delete [] parallelfilecurrent;
  delete [] multiname_ex;

  delete [] idregion;
  memory->destroy(thresh_array);
  memory->destroy(thresh_op);
  memory->destroy(thresh_value);

  for (int i = 0; i < ncompute; i++) delete [] id_compute[i];
  memory->sfree(id_compute);
  delete [] compute;

  for (int i = 0; i < nfix; i++) delete [] id_fix[i];
  memory->sfree(id_fix);
  delete [] fix;

  for (int i = 0; i < nvariable; i++) delete [] id_variable[i];
  memory->sfree(id_variable);
  delete [] variable;
  for (int i = 0; i < nvariable; i++) memory->destroy(vbuf[i]);
  delete [] vbuf;

  for (int i = 0; i < ncustom; i++) delete [] id_custom[i];
  memory->sfree(id_custom);
  delete [] flag_custom;

  memory->destroy(choose);
  memory->destroy(dchoose);
  memory->destroy(clist);

  if (typenames) {
    for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
    delete [] typenames;
  }

  delete [] label; //NP modified C.K.
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::init_style()
{
  // default for element names = C

  if (typenames == NULL) {
    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      typenames[itype] = new char[2];
      strcpy(typenames[itype],"C");
    }
  }

  // setup boundary string

  domain->boundary_string(boundstr);

  // setup function ptrs

  header_choice = &DumpCustomVTK::header_vtk;

  if (vtk_file_format == VTP || vtk_file_format == PVTP)
    write_choice = &DumpCustomVTK::write_vtp;
  else if (vtk_file_format == VTU || vtk_file_format == PVTU)
    write_choice = &DumpCustomVTK::write_vtu;
  else
    write_choice = &DumpCustomVTK::write_vtk;

  // find current ptr for each compute,fix,variable
  // check that fix frequency is acceptable

  int icompute;
  for (int i = 0; i < ncompute; i++) {
    icompute = modify->find_compute(id_compute[i]);
    if (icompute < 0) error->all(FLERR,"Could not find dump custom/vtk compute ID");
    compute[i] = modify->compute[icompute];
  }

  int ifix;
  for (int i = 0; i < nfix; i++) {
    ifix = modify->find_fix(id_fix[i]);
    if (ifix < 0) error->all(FLERR,"Could not find dump custom/vtk fix ID");
    fix[i] = modify->fix[ifix];
    if (nevery % modify->fix[ifix]->peratom_freq)
      error->all(FLERR,"Dump custom/vtk and fix not computed at compatible times");
  }

  int ivariable;
  for (int i = 0; i < nvariable; i++) {
    ivariable = input->variable->find(id_variable[i]);
    if (ivariable < 0)
      error->all(FLERR,"Could not find dump custom/vtk variable name");
    variable[i] = ivariable;
  }

  int icustom;
  for (int i = 0; i < ncustom; i++) {
    icustom = atom->find_custom(id_custom[i],flag_custom[i]);
    if (icustom < 0)
      error->all(FLERR,"Could not find custom per-atom property ID");
  }

  // set index and check validity of region

  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for dump custom/vtk does not exist");
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_header(bigint /*ndump*/)
{
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::header_vtk(bigint)
{
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTK::count()
{
  n_calls_ = 0;

  int i;

  // grow choose and variable vbuf arrays if needed

  int nlocal = atom->nlocal;
  if (nlocal > maxlocal) {
    maxlocal = atom->nmax;

    memory->destroy(choose);
    memory->destroy(dchoose);
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(dchoose,maxlocal,"dump:dchoose");
    memory->create(clist,maxlocal,"dump:clist");

    for (i = 0; i < nvariable; i++) {
      memory->destroy(vbuf[i]);
      memory->create(vbuf[i],maxlocal,"dump:vbuf");
    }
  }

  // invoke Computes for per-atom quantities

  if (ncompute) {
    for (i = 0; i < ncompute; i++)
      if (!(compute[i]->invoked_flag & INVOKED_PERATOM)) {
        compute[i]->compute_peratom();
        compute[i]->invoked_flag |= INVOKED_PERATOM;
      }
  }

  // evaluate atom-style Variables for per-atom quantities

  if (nvariable)
    for (i = 0; i < nvariable; i++)
      input->variable->compute_atom(variable[i],igroup,vbuf[i],1,0);

  // choose all local atoms for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in group

  if (igroup) {
    int *mask = atom->mask;
    for (i = 0; i < nlocal; i++)
      if (!(mask[i] & groupbit))
        choose[i] = 0;
  }

  // un-choose if not in region

  if (iregion >= 0) {
    Region *region = domain->regions[iregion];
    double **x = atom->x;
    for (i = 0; i < nlocal; i++)
      if (choose[i] && region->match(x[i][0],x[i][1],x[i][2]) == 0)
        choose[i] = 0;
  }

  // un-choose if any threshold criterion isn't met

  if (nthresh) {
    double *ptr;
    double value;
    int nstride;
    int nlocal = atom->nlocal;

    for (int ithresh = 0; ithresh < nthresh; ithresh++) {

      // customize by adding to if statement

      if (thresh_array[ithresh] == ID) {
        int *tag = atom->tag;
        for (i = 0; i < nlocal; i++) dchoose[i] = tag[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MOL) {
        if (!atom->molecule_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        int *molecule = atom->molecule;
        for (i = 0; i < nlocal; i++) dchoose[i] = molecule[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == TYPE) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ELEMENT) {
        int *type = atom->type;
        for (i = 0; i < nlocal; i++) dchoose[i] = type[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == MASS) {
        if (atom->rmass) {
          ptr = atom->rmass;
          nstride = 1;
        } else {
          double *mass = atom->mass;
          int *type = atom->type;
          for (i = 0; i < nlocal; i++) dchoose[i] = mass[type[i]];
          ptr = dchoose;
          nstride = 1;
        }

      } else if (thresh_array[ithresh] == X) {
        ptr = &atom->x[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == Y) {
        ptr = &atom->x[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == Z) {
        ptr = &atom->x[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == XS) {
        double **x = atom->x;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YS) {
        double **x = atom->x;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][1] - boxylo) * invyprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZS) {
        double **x = atom->x;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) + h_inv[4]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSTRI) {
        double **x = atom->x;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]);
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double xprd = domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][0] + ((image[i] & IMGMASK) - IMGMAX) * xprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double yprd = domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][1] +
            ((image[i] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double zprd = domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = x[i][2] + ((image[i] >> IMG2BITS) - IMGMAX) * zprd;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *h = domain->h;
        int xbox,ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          xbox = (image[i] & IMGMASK) - IMGMAX;
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *h = domain->h;
        int ybox,zbox;
        for (i = 0; i < nlocal; i++) {
          ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][1] + h[1]*ybox + h[3]*zbox;
        }
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *h = domain->h;
        int zbox;
        for (i = 0; i < nlocal; i++) {
          zbox = (image[i] >> IMG2BITS) - IMGMAX;
          dchoose[i] = x[i][2] + h[2]*zbox;
        }
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double boxxlo = domain->boxlo[0];
        double invxprd = 1.0/domain->xprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][0] - boxxlo) * invxprd +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == YSU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double boxylo = domain->boxlo[1];
        double invyprd = 1.0/domain->yprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] =
            (x[i][1] - boxylo) * invyprd +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == ZSU) {
        double **x = atom->x;
        tagint *image = atom->image;
        double boxzlo = domain->boxlo[2];
        double invzprd = 1.0/domain->zprd;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (x[i][2] - boxzlo) * invzprd +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == XSUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[0]*(x[i][0]-boxlo[0]) +
            h_inv[5]*(x[i][1]-boxlo[1]) +
            h_inv[4]*(x[i][2]-boxlo[2]) +
            (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == YSUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[1]*(x[i][1]-boxlo[1]) +
            h_inv[3]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ZSUTRI) {
        double **x = atom->x;
        tagint *image = atom->image;
        double *boxlo = domain->boxlo;
        double *h_inv = domain->h_inv;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = h_inv[2]*(x[i][2]-boxlo[2]) +
            (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == IX) {
        tagint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IY) {
        tagint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == IZ) {
        tagint *image = atom->image;
        for (i = 0; i < nlocal; i++)
          dchoose[i] = (image[i] >> IMG2BITS) - IMGMAX;
        ptr = dchoose;
        nstride = 1;

      } else if (thresh_array[ithresh] == VX) {
        ptr = &atom->v[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == VY) {
        ptr = &atom->v[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == VZ) {
        ptr = &atom->v[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == FX) {
        ptr = &atom->f[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == FY) {
        ptr = &atom->f[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == FZ) {
        ptr = &atom->f[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == Q) {
        if (!atom->q_flag)
          error->all(FLERR,"Threshhold for an atom property that isn't allocated");
        ptr = atom->q;
        nstride = 1;
      } else if (thresh_array[ithresh] == P) { //NP modified C.K.
        if (!atom->p_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = atom->p;
        nstride = 1;
      } else if (thresh_array[ithresh] == RHO) { //NP modified C.K.
        if (!atom->rho_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->rho;
        nstride = 1;
      } else if (thresh_array[ithresh] == DENSITY) { //NP modified C.K.
        if (!atom->density_flag)
          error->all(FLERR,"Threshold for an atom property that isn't allocated");
        ptr = atom->density;
        nstride = 1;
      } else if (thresh_array[ithresh] == MUX) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][0];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUY) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][1];
        nstride = 4;
      } else if (thresh_array[ithresh] == MUZ) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][2];
        nstride = 4;
      } else if (thresh_array[ithresh] == MU) {
        if (!atom->mu_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->mu[0][3];
        nstride = 4;

      } else if (thresh_array[ithresh] == RADIUS) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->radius;
        nstride = 1;
      } else if (thresh_array[ithresh] == DIAMETER) {
        if (!atom->radius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        double *radius = atom->radius;
        for (i = 0; i < nlocal; i++) dchoose[i] = 2.0*radius[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == OMEGAX) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAY) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == OMEGAZ) {
        if (!atom->omega_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->omega[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMX) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMY) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == ANGMOMZ) {
        if (!atom->angmom_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->angmom[0][2];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQX) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][0];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQY) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][1];
        nstride = 3;
      } else if (thresh_array[ithresh] == TQZ) {
        if (!atom->torque_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = &atom->torque[0][2];
        nstride = 3;

      } else if (thresh_array[ithresh] == SPIN) {
        if (!atom->spin_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        int *spin = atom->spin;
        for (i = 0; i < nlocal; i++) dchoose[i] = spin[i];
        ptr = dchoose;
        nstride = 1;
      } else if (thresh_array[ithresh] == ERADIUS) {
        if (!atom->eradius_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->eradius;
        nstride = 1;
      } else if (thresh_array[ithresh] == ERVEL) {
        if (!atom->ervel_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->ervel;
        nstride = 1;
      } else if (thresh_array[ithresh] == ERFORCE) {
        if (!atom->erforce_flag)
          error->all(FLERR,
                     "Threshold for an atom property that isn't allocated");
        ptr = atom->erforce;
        nstride = 1;

      } else if (thresh_array[ithresh] == COMPUTE) {
        i = ATTRIBUTES + nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = compute[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &compute[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = compute[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == FIX) {
        i = ATTRIBUTES + nfield + ithresh;
        if (argindex[i] == 0) {
          ptr = fix[field2index[i]]->vector_atom;
          nstride = 1;
        } else {
          ptr = &fix[field2index[i]]->array_atom[0][argindex[i]-1];
          nstride = fix[field2index[i]]->size_peratom_cols;
        }

      } else if (thresh_array[ithresh] == VARIABLE) {
        i = ATTRIBUTES + nfield + ithresh;
        ptr = vbuf[field2index[i]];
        nstride = 1;

      } else if (thresh_array[ithresh] == DNAME) {
        int iwhich,tmp;
        i = ATTRIBUTES + nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);
        ptr = atom->dvector[iwhich];
        nstride = 1;

      } else if (thresh_array[ithresh] == INAME) {
        int iwhich,tmp;
        i = ATTRIBUTES + nfield + ithresh;
        iwhich = atom->find_custom(id_custom[field2index[i]],tmp);

        int *ivector = atom->ivector[iwhich];
        for (i = 0; i < nlocal; i++)
          dchoose[i] = ivector[i];
        ptr = dchoose;
        nstride = 1;
      }
      else if (thresh_array[ithresh] == THREAD) {
        int *thread = atom->thread;
        for (i = 0; i < nlocal; i++) dchoose[i] = comm->me * comm->nprocs + thread[i];
        ptr = dchoose;
        nstride = 1;
      }

      // unselect atoms that don't meet threshold criterion

      value = thresh_value[ithresh];

      switch (thresh_op[ithresh]) {
      case LT:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr >= value) choose[i] = 0;
        break;
      case LE:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr > value) choose[i] = 0;
        break;
      case GT:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr <= value) choose[i] = 0;
        break;
      case GE:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr < value) choose[i] = 0;
        break;
      case EQ:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr != value) choose[i] = 0;
        break;
      case NEQ:
        for (i = 0; i < nlocal; i++, ptr += nstride)
          if (choose[i] && *ptr == value) choose[i] = 0;
        break;
      }
    }
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  return nchoose;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write()
{
  // simulation box bounds

  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    domain->box_corners();
    boxcorners = domain->corners;
  }

  // nme = # of dump lines this proc contributes to dump

  nme = count();

  // ntotal = total # of dump lines in snapshot
  // nmax = max # of dump lines on any proc

  bigint bnme = nme;
  MPI_Allreduce(&bnme,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  int nmax;
  if (multiproc != nprocs) MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  else nmax = nme;

  // write timestep header
  // for multiproc,
  //   nheader = # of lines in this file via Allreduce on clustercomm

  bigint nheader = ntotal;
  if (multiproc)
    MPI_Allreduce(&bnme,&nheader,1,MPI_LMP_BIGINT,MPI_SUM,clustercomm);

  if (filewriter) write_header(nheader);

  // insure buf is sized for packing and communicating
  // use nmax to insure filewriter proc can receive info from others
  // limit nmax*size_one to int since used as arg in MPI calls

  if (nmax > maxbuf) {
    if ((bigint) nmax * size_one > MAXSMALLINT)
      error->all(FLERR,"Too much per-proc info for dump");
    maxbuf = nmax;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // insure ids buffer is sized for sorting

  if (sort_flag && sortcol == 0 && nmax > maxids) {
    maxids = nmax;
    memory->destroy(ids);
    memory->create(ids,maxids,"dump:ids");
  }

  // pack my data into buf
  // if sorting on IDs also request ID list from pack()
  // sort buf as needed

  if (sort_flag && sortcol == 0) pack(ids);
  else pack(NULL);
  if (sort_flag) sort();

  // filewriter = 1 = this proc writes to file
  //   ping each proc in my cluster, receive its data, write data to file
  // else wait for ping from fileproc, send my data to fileproc

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;

  // comm and output buf of doubles

  if (filewriter) {
    for (int iproc = 0; iproc < nclusterprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(buf,maxbuf*size_one,MPI_DOUBLE,me+iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,me+iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&nlines);
        nlines /= size_one;
      } else nlines = nme;

      write_data(nlines,buf);
    }
  } else {
    MPI_Recv(&tmp,0,MPI_INT,fileproc,0,world,&status);
    MPI_Rsend(buf,nme*size_one,MPI_DOUBLE,fileproc,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack(int *ids)
{
  int n = 0;
  for (std::map<int,FnPtrPack>::iterator it=pack_choice.begin(); it!=pack_choice.end(); ++it, ++n) {
      current_pack_choice_key = it->first; // work-around for pack_compute, pack_fix, pack_variable
      (this->*(it->second))(n);
  }

  if (ids) {
    int *tag = atom->tag;
    for (int i = 0; i < nchoose; i++)
      ids[i] = tag[clist[i]];
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::setFileCurrent() {
  delete [] filecurrent;
  filecurrent = NULL;

  char *filestar = filename;
  if (multiproc) {
    if (multiproc > 1) { // if dump_modify fileper or nfile was used
      delete [] multiname_ex;
      multiname_ex = NULL;
      char *ptr = strchr(filename,'%');
      if (ptr) {
        int id;
        if (me + nclusterprocs == nprocs) // last filewriter
          id = multiproc -1;
        else
          id = me/nclusterprocs;
        multiname_ex = new char[strlen(filename) + 16];
        *ptr = '\0';
        sprintf(multiname_ex,"%s_%d%s",filename,id,ptr+1);
        *ptr = '%';
      }
    } // else multiname_ex built in constructor is OK
    filestar = multiname_ex;
  }

  if (multifile == 0) {
    filecurrent = new char[strlen(filestar) + 1];
    strcpy(filecurrent, filestar);
  } else {
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0) {
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    } else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
  }

  // filename of domain box data file
  delete [] domainfilecurrent;
  domainfilecurrent = NULL;
  if (multiproc) {
    // remove '%' character
    char *ptr = strchr(filename,'%');
    domainfilecurrent = new char[strlen(filename)];
    *ptr = '\0';
    sprintf(domainfilecurrent,"%s%s",filename,ptr+1);
    *ptr = '%';
    // insert "_boundingBox" string
    ptr = strrchr(domainfilecurrent,'.');
    filestar = new char[strlen(domainfilecurrent)+16];
    *ptr = '\0';
    sprintf(filestar,"%s_boundingBox.%s",domainfilecurrent,ptr+1);
    delete [] domainfilecurrent;
    domainfilecurrent = NULL;

    if (multifile == 0) {
      domainfilecurrent = new char[strlen(filestar) + 1];
      strcpy(domainfilecurrent, filestar);
    } else {
      domainfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(domainfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(domainfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  } else {
    domainfilecurrent = new char[strlen(filecurrent) + 16];
    char *ptr = strrchr(filecurrent,'.');
    *ptr = '\0';
    sprintf(domainfilecurrent,"%s_boundingBox.%s",filecurrent,ptr+1);
    *ptr = '.';
  }

  // filename of parallel file
  if (multiproc && me == 0) {
    delete [] parallelfilecurrent;
    parallelfilecurrent = NULL;

    // remove '%' character and add 'p' to file extension
    // -> string length stays the same
    char *ptr = strchr(filename,'%');
    filestar = new char[strlen(filename) + 1];
    *ptr = '\0';
    sprintf(filestar,"%s%s",filename,ptr+1);
    *ptr = '%';
    ptr = strrchr(filestar,'.');
    ptr++;
    *ptr++='p';
    *ptr++='v';
    *ptr++='t';
    *ptr++= (vtk_file_format == PVTP)?'p':'u';
    *ptr++= 0;

    if (multifile == 0) {
      parallelfilecurrent = new char[strlen(filestar) + 1];
      strcpy(parallelfilecurrent, filestar);
    } else {
      parallelfilecurrent = new char[strlen(filestar) + 16];
      char *ptr = strchr(filestar,'*');
      *ptr = '\0';
      if (padflag == 0) {
        sprintf(parallelfilecurrent,"%s" BIGINT_FORMAT "%s",
                filestar,update->ntimestep,ptr+1);
      } else {
        char bif[8],pad[16];
        strcpy(bif,BIGINT_FORMAT);
        sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
        sprintf(parallelfilecurrent,pad,filestar,update->ntimestep,ptr+1);
      }
      *ptr = '*';
    }
    delete [] filestar;
    filestar = NULL;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::buf2arrays(int n, double *mybuf)
{
  for (int iatom=0; iatom < n; ++iatom) {
    vtkIdType pid[1];
    pid[0] = points->InsertNextPoint(mybuf[iatom*size_one],mybuf[iatom*size_one+1],mybuf[iatom*size_one+2]);

    int j=3; // 0,1,2 = x,y,z handled just above
    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      vtkAbstractArray *paa = it->second;
      if (it->second->GetNumberOfComponents() == 3) {
        switch (vtype[it->first]) {
          case INT:
            {
              int iv3[3] = { static_cast<int>(mybuf[iatom*size_one+j  ]),
                             static_cast<int>(mybuf[iatom*size_one+j+1]),
                             static_cast<int>(mybuf[iatom*size_one+j+2]) };
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);

#if (VTK_MAJOR_VERSION > 7) || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION     >= 1)
              pia->InsertNextTypedTuple(iv3);
#else
              pia->InsertNextTupleValue(iv3);
#endif
              break;
            }
          case DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
#if (VTK_MAJOR_VERSION > 7) || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION     >= 1)
              pda->InsertNextTypedTuple(&mybuf[iatom*size_one+j]);
#else
              pda->InsertNextTupleValue(&mybuf[iatom*size_one+j]);
#endif
              break;
            }
        }
        j+=3;
      } else {
        switch (vtype[it->first]) {
          case INT:
            {
              vtkIntArray *pia = static_cast<vtkIntArray*>(paa);
              pia->InsertNextValue(mybuf[iatom*size_one+j]);
              break;
            }
          case DOUBLE:
            {
              vtkDoubleArray *pda = static_cast<vtkDoubleArray*>(paa);
              pda->InsertNextValue(mybuf[iatom*size_one+j]);
              break;
            }
          case STRING:
            {
              vtkStringArray *psa = static_cast<vtkStringArray*>(paa);
              psa->InsertNextValue(typenames[static_cast<int>(mybuf[iatom*size_one+j])]);
              break;
            }
        }
        ++j;
      }
    }

    pointsCells->InsertNextCell(1,pid);
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::prepare_domain_data(vtkRectilinearGrid *rgrid)
{
  vtkSmartPointer<vtkDoubleArray> xCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  xCoords->InsertNextValue(boxxlo);
  xCoords->InsertNextValue(boxxhi);
  vtkSmartPointer<vtkDoubleArray> yCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  yCoords->InsertNextValue(boxylo);
  yCoords->InsertNextValue(boxyhi);
  vtkSmartPointer<vtkDoubleArray> zCoords =  vtkSmartPointer<vtkDoubleArray>::New();
  zCoords->InsertNextValue(boxzlo);
  zCoords->InsertNextValue(boxzhi);

  rgrid->SetDimensions(2,2,2);
  rgrid->SetXCoordinates(xCoords);
  rgrid->SetYCoordinates(yCoords);
  rgrid->SetZCoordinates(zCoords);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::prepare_domain_data_triclinic(vtkUnstructuredGrid *hexahedronGrid)
{
  vtkSmartPointer<vtkPoints> hexahedronPoints = vtkSmartPointer<vtkPoints>::New();
  hexahedronPoints->SetNumberOfPoints(8);
  hexahedronPoints->InsertPoint(0, boxcorners[0][0], boxcorners[0][1], boxcorners[0][2]);
  hexahedronPoints->InsertPoint(1, boxcorners[1][0], boxcorners[1][1], boxcorners[1][2]);
  hexahedronPoints->InsertPoint(2, boxcorners[3][0], boxcorners[3][1], boxcorners[3][2]);
  hexahedronPoints->InsertPoint(3, boxcorners[2][0], boxcorners[2][1], boxcorners[2][2]);
  hexahedronPoints->InsertPoint(4, boxcorners[4][0], boxcorners[4][1], boxcorners[4][2]);
  hexahedronPoints->InsertPoint(5, boxcorners[5][0], boxcorners[5][1], boxcorners[5][2]);
  hexahedronPoints->InsertPoint(6, boxcorners[7][0], boxcorners[7][1], boxcorners[7][2]);
  hexahedronPoints->InsertPoint(7, boxcorners[6][0], boxcorners[6][1], boxcorners[6][2]);
  vtkSmartPointer<vtkHexahedron> hexahedron = vtkSmartPointer<vtkHexahedron>::New();
  hexahedron->GetPointIds()->SetId(0, 0);
  hexahedron->GetPointIds()->SetId(1, 1);
  hexahedron->GetPointIds()->SetId(2, 2);
  hexahedron->GetPointIds()->SetId(3, 3);
  hexahedron->GetPointIds()->SetId(4, 4);
  hexahedron->GetPointIds()->SetId(5, 5);
  hexahedron->GetPointIds()->SetId(6, 6);
  hexahedron->GetPointIds()->SetId(7, 7);

  hexahedronGrid->Allocate(1, 1);
  hexahedronGrid->InsertNextCell(hexahedron->GetCellType(),
                                  hexahedron->GetPointIds());
  hexahedronGrid->SetPoints(hexahedronPoints);
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtk()
{
  vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  prepare_domain_data(rgrid.GetPointer());

  vtkSmartPointer<vtkRectilinearGridWriter> gwriter = vtkSmartPointer<vtkRectilinearGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LIGGGHTS");

  if (binary) gwriter->SetFileTypeToBinary();
  else        gwriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(rgrid);
#else
  gwriter->SetInputData(rgrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtk_triclinic()
{
  vtkSmartPointer<vtkUnstructuredGrid> hexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  prepare_domain_data_triclinic(hexahedronGrid.GetPointer());

  vtkSmartPointer<vtkUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();

  if(label) gwriter->SetHeader(label);
  else      gwriter->SetHeader("Generated by LIGGGHTS");

  if (binary) gwriter->SetFileTypeToBinary();
  else        gwriter->SetFileTypeToASCII();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtr()
{
  vtkSmartPointer<vtkRectilinearGrid> rgrid = vtkSmartPointer<vtkRectilinearGrid>::New();
  prepare_domain_data(rgrid.GetPointer());

  vtkSmartPointer<vtkXMLRectilinearGridWriter> gwriter = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();

  if (binary) gwriter->SetDataModeToBinary();
  else        gwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(rgrid);
#else
  gwriter->SetInputData(rgrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_domain_vtu_triclinic()
{
  vtkSmartPointer<vtkUnstructuredGrid> hexahedronGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  prepare_domain_data_triclinic(hexahedronGrid.GetPointer());

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> gwriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

  if (binary) gwriter->SetDataModeToBinary();
  else        gwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
  gwriter->SetInput(hexahedronGrid);
#else
  gwriter->SetInputData(hexahedronGrid);
#endif
  gwriter->SetFileName(domainfilecurrent);
  gwriter->Write();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtk(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but only proc 0 is a filewriter (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
#ifdef UNSTRUCTURED_GRID_VTK
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VERTEX, pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
#else
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetVerts(pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
#endif

    if(label) writer->SetHeader(label);
    else      writer->SetHeader("Generated by LIGGGHTS");

    if (binary) writer->SetFileTypeToBinary();
    else        writer->SetFileTypeToASCII();

#ifdef UNSTRUCTURED_GRID_VTK
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
  #else
    writer->SetInputData(unstructuredGrid);
  #endif
#else
  #if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
  #else
    writer->SetInputData(polyData);
  #endif
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (write_domain) {
      if (domain->triclinic == 0)
        write_domain_vtk();
      else
        write_domain_vtk_triclinic();
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtp(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

    polyData->SetPoints(points);
    polyData->SetVerts(pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      polyData->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(polyData);
#else
    writer->SetInputData(polyData);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPPolyDataWriter> pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(polyData);
#else
        pwriter->SetInputData(polyData);
#endif
        pwriter->Write();
      }

      if (write_domain) {
        if (domain->triclinic == 0) {
          domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
          write_domain_vtr();
        } else {
          domainfilecurrent[strlen(domainfilecurrent)-1] = 'u'; // adjust filename extension
          write_domain_vtu_triclinic();
        }
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::write_vtu(int n, double *mybuf)
{
  ++n_calls_;

  buf2arrays(n, mybuf);

  if (n_calls_ < nclusterprocs)
    return; // multiple processors but not all are filewriters (-> nclusterprocs procs contribute to the filewriter's output data)

  setFileCurrent();

  {
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_VERTEX, pointsCells);

    for (std::map<int, vtkSmartPointer<vtkAbstractArray> >::iterator it=myarrays.begin(); it!=myarrays.end(); ++it) {
      unstructuredGrid->GetPointData()->AddArray(it->second);
    }

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    if (binary) writer->SetDataModeToBinary();
    else        writer->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->SetFileName(filecurrent);
    writer->Write();

    if (me == 0) {
      if (multiproc) {
        vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
        pwriter->SetFileName(parallelfilecurrent);
        pwriter->SetNumberOfPieces((multiproc > 1)?multiproc:nprocs);
        if (binary) pwriter->SetDataModeToBinary();
        else        pwriter->SetDataModeToAscii();

#if VTK_MAJOR_VERSION < 6
        pwriter->SetInput(unstructuredGrid);
#else
        pwriter->SetInputData(unstructuredGrid);
#endif
        pwriter->Write();
      }

      if (write_domain) {
        if (domain->triclinic == 0) {
          domainfilecurrent[strlen(domainfilecurrent)-1] = 'r'; // adjust filename extension
          write_domain_vtr();
        } else {
          write_domain_vtu_triclinic();
        }
      }
    }
  }

  reset_vtk_data_containers();
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::reset_vtk_data_containers()
{
  points = vtkSmartPointer<vtkPoints>::New();
  pointsCells = vtkSmartPointer<vtkCellArray>::New();

  std::map<int,int>::iterator it=vtype.begin();
  ++it; ++it; ++it;
  for (; it!=vtype.end(); ++it) {
    switch(vtype[it->first]) {
      case INT:
        myarrays[it->first] = vtkSmartPointer<vtkIntArray>::New();
        break;
      case DOUBLE:
        myarrays[it->first] = vtkSmartPointer<vtkDoubleArray>::New();
        break;
      case STRING:
        myarrays[it->first] = vtkSmartPointer<vtkStringArray>::New();
        break;
    }

    if (vector_set.find(it->first) != vector_set.end()) {
      myarrays[it->first]->SetNumberOfComponents(3);
      myarrays[it->first]->SetName(name[it->first].c_str());
      ++it; ++it;
    } else {
      myarrays[it->first]->SetName(name[it->first].c_str());
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTK::parse_fields(int narg, char **arg)
{

  pack_choice[X] = &DumpCustomVTK::pack_x;
  vtype[X] = DOUBLE;
  name[X] = "x";
  pack_choice[Y] = &DumpCustomVTK::pack_y;
  vtype[Y] = DOUBLE;
  name[Y] = "y";
  pack_choice[Z] = &DumpCustomVTK::pack_z;
  vtype[Z] = DOUBLE;
  name[Z] = "z";

  // customize by adding to if statement
  int i;
  for (int iarg = 5; iarg < narg; iarg++) {
    i = iarg-5;

    if (strcmp(arg[iarg],"id") == 0) {
      pack_choice[ID] = &DumpCustomVTK::pack_id;
      vtype[ID] = INT;
      name[ID] = arg[iarg];
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (!atom->molecule_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MOL] = &DumpCustomVTK::pack_molecule;
      vtype[MOL] = INT;
      name[MOL] = arg[iarg];
    } else if (strcmp(arg[iarg],"type") == 0) {
      pack_choice[TYPE] = &DumpCustomVTK::pack_type;
      vtype[TYPE] = INT;
      name[TYPE] =arg[iarg];
    } else if (strcmp(arg[iarg],"element") == 0) {
      pack_choice[ELEMENT] = &DumpCustomVTK::pack_type;
      vtype[ELEMENT] = STRING;
      name[ELEMENT] = arg[iarg];
    } else if (strcmp(arg[iarg],"mass") == 0) {
      pack_choice[MASS] = &DumpCustomVTK::pack_mass;
      vtype[MASS] = DOUBLE;
      name[MASS] = arg[iarg];

    } else if (strcmp(arg[iarg],"x") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"y") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"z") == 0) {
      // required property
    } else if (strcmp(arg[iarg],"xs") == 0) {
      if (domain->triclinic) pack_choice[XS] = &DumpCustomVTK::pack_xs_triclinic;
      else pack_choice[XS] = &DumpCustomVTK::pack_xs;
      vtype[XS] = DOUBLE;
      name[XS] = arg[iarg];
    } else if (strcmp(arg[iarg],"ys") == 0) {
      if (domain->triclinic) pack_choice[YS] = &DumpCustomVTK::pack_ys_triclinic;
      else pack_choice[YS] = &DumpCustomVTK::pack_ys;
      vtype[YS] = DOUBLE;
      name[YS] = arg[iarg];
    } else if (strcmp(arg[iarg],"zs") == 0) {
      if (domain->triclinic) pack_choice[ZS] = &DumpCustomVTK::pack_zs_triclinic;
      else pack_choice[ZS] = &DumpCustomVTK::pack_zs;
      vtype[ZS] = DOUBLE;
      name[ZS] = arg[iarg];
    } else if (strcmp(arg[iarg],"xu") == 0) {
      if (domain->triclinic) pack_choice[XU] = &DumpCustomVTK::pack_xu_triclinic;
      else pack_choice[XU] = &DumpCustomVTK::pack_xu;
      vtype[XU] = DOUBLE;
      name[XU] = arg[iarg];
    } else if (strcmp(arg[iarg],"yu") == 0) {
      if (domain->triclinic) pack_choice[YU] = &DumpCustomVTK::pack_yu_triclinic;
      else pack_choice[YU] = &DumpCustomVTK::pack_yu;
      vtype[YU] = DOUBLE;
      name[YU] = arg[iarg];
    } else if (strcmp(arg[iarg],"zu") == 0) {
      if (domain->triclinic) pack_choice[ZU] = &DumpCustomVTK::pack_zu_triclinic;
      else pack_choice[ZU] = &DumpCustomVTK::pack_zu;
      vtype[ZU] = DOUBLE;
      name[ZU] = arg[iarg];
    } else if (strcmp(arg[iarg],"xsu") == 0) {
      if (domain->triclinic) pack_choice[XSU] = &DumpCustomVTK::pack_xsu_triclinic;
      else pack_choice[XSU] = &DumpCustomVTK::pack_xsu;
      vtype[XSU] = DOUBLE;
      name[XSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"ysu") == 0) {
      if (domain->triclinic) pack_choice[YSU] = &DumpCustomVTK::pack_ysu_triclinic;
      else pack_choice[YSU] = &DumpCustomVTK::pack_ysu;
      vtype[YSU] = DOUBLE;
      name[YSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"zsu") == 0) {
      if (domain->triclinic) pack_choice[ZSU] = &DumpCustomVTK::pack_zsu_triclinic;
      else pack_choice[ZSU] = &DumpCustomVTK::pack_zsu;
      vtype[ZSU] = DOUBLE;
      name[ZSU] = arg[iarg];
    } else if (strcmp(arg[iarg],"ix") == 0) {
      pack_choice[IX] = &DumpCustomVTK::pack_ix;
      vtype[IX] = INT;
      name[IX] = arg[iarg];
    } else if (strcmp(arg[iarg],"iy") == 0) {
      pack_choice[IY] = &DumpCustomVTK::pack_iy;
      vtype[IY] = INT;
      name[IY] = arg[iarg];
    } else if (strcmp(arg[iarg],"iz") == 0) {
      pack_choice[IZ] = &DumpCustomVTK::pack_iz;
      vtype[IZ] = INT;
      name[IZ] = arg[iarg];

    } else if (strcmp(arg[iarg],"vx") == 0) {
      pack_choice[VX] = &DumpCustomVTK::pack_vx;
      vtype[VX] = DOUBLE;
      name[VX] = arg[iarg];
    } else if (strcmp(arg[iarg],"vy") == 0) {
      pack_choice[VY] = &DumpCustomVTK::pack_vy;
      vtype[VY] = DOUBLE;
      name[VY] = arg[iarg];
    } else if (strcmp(arg[iarg],"vz") == 0) {
      pack_choice[VZ] = &DumpCustomVTK::pack_vz;
      vtype[VZ] = DOUBLE;
      name[VZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"fx") == 0) {
      pack_choice[FX] = &DumpCustomVTK::pack_fx;
      vtype[FX] = DOUBLE;
      name[FX] = arg[iarg];
    } else if (strcmp(arg[iarg],"fy") == 0) {
      pack_choice[FY] = &DumpCustomVTK::pack_fy;
      vtype[FY] = DOUBLE;
      name[FY] = arg[iarg];
    } else if (strcmp(arg[iarg],"fz") == 0) {
      pack_choice[FZ] = &DumpCustomVTK::pack_fz;
      vtype[FZ] = DOUBLE;
      name[FZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"q") == 0) {
      if (!atom->q_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[Q] = &DumpCustomVTK::pack_q;
      vtype[Q] = DOUBLE;
      name[Q] = arg[iarg];
   } else if (strcmp(arg[iarg],"density") == 0) {
      if (!atom->density_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[DENSITY] = &DumpCustomVTK::pack_density;
      vtype[DENSITY] = DOUBLE;
      name[DENSITY] = arg[iarg];
   } else if (strcmp(arg[iarg],"p") == 0) {
      if (!atom->p_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[P] = &DumpCustomVTK::pack_p;
      vtype[P] = DOUBLE;
      name[P] = arg[iarg];
   } else if (strcmp(arg[iarg],"rho") == 0) {
      if (!atom->rho_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[RHO] = &DumpCustomVTK::pack_rho;
      vtype[RHO] = DOUBLE;
      name[RHO] = arg[iarg];
    } else if (strcmp(arg[iarg],"mux") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUX] = &DumpCustomVTK::pack_mux;
      vtype[MUX] = DOUBLE;
      name[MUX] = arg[iarg];
    } else if (strcmp(arg[iarg],"muy") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUY] = &DumpCustomVTK::pack_muy;
      vtype[MUY] = DOUBLE;
      name[MUY] = arg[iarg];
    } else if (strcmp(arg[iarg],"muz") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MUZ] = &DumpCustomVTK::pack_muz;
      vtype[MUZ] = DOUBLE;
      name[MUZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"mu") == 0) {
      if (!atom->mu_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[MU] = &DumpCustomVTK::pack_mu;
      vtype[MU] = DOUBLE;
      name[MU] = arg[iarg];

    } else if (strcmp(arg[iarg],"radius") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[RADIUS] = &DumpCustomVTK::pack_radius;
      vtype[RADIUS] = DOUBLE;
      name[RADIUS] = arg[iarg];
    } else if (strcmp(arg[iarg],"diameter") == 0) {
      if (!atom->radius_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[DIAMETER] = &DumpCustomVTK::pack_diameter;
      vtype[DIAMETER] = DOUBLE;
      name[DIAMETER] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegax") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAX] = &DumpCustomVTK::pack_omegax;
      vtype[OMEGAX] = DOUBLE;
      name[OMEGAX] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegay") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAY] = &DumpCustomVTK::pack_omegay;
      vtype[OMEGAY] = DOUBLE;
      name[OMEGAY] = arg[iarg];
    } else if (strcmp(arg[iarg],"omegaz") == 0) {
      if (!atom->omega_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[OMEGAZ] = &DumpCustomVTK::pack_omegaz;
      vtype[OMEGAZ] = DOUBLE;
      name[OMEGAZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomx") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMX] = &DumpCustomVTK::pack_angmomx;
      vtype[ANGMOMX] = DOUBLE;
      name[ANGMOMX] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomy") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMY] = &DumpCustomVTK::pack_angmomy;
      vtype[ANGMOMY] = DOUBLE;
      name[ANGMOMY] = arg[iarg];
    } else if (strcmp(arg[iarg],"angmomz") == 0) {
      if (!atom->angmom_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[ANGMOMZ] = &DumpCustomVTK::pack_angmomz;
      vtype[ANGMOMZ] = DOUBLE;
      name[ANGMOMZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqx") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQX] = &DumpCustomVTK::pack_tqx;
      vtype[TQX] = DOUBLE;
      name[TQX] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqy") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQY] = &DumpCustomVTK::pack_tqy;
      vtype[TQY] = DOUBLE;
      name[TQY] = arg[iarg];
    } else if (strcmp(arg[iarg],"tqz") == 0) {
      if (!atom->torque_flag)
        error->all(FLERR,"Dumping an atom property that isn't allocated");
      pack_choice[TQZ] = &DumpCustomVTK::pack_tqz;
      vtype[TQZ] = DOUBLE;
      name[TQZ] = arg[iarg];

    } else if (strcmp(arg[iarg],"spin") == 0) {
      if (!atom->spin_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[SPIN] = &DumpCustomVTK::pack_spin;
      vtype[SPIN] = INT;
      name[SPIN] = arg[iarg];
    } else if (strcmp(arg[iarg],"eradius") == 0) {
      if (!atom->eradius_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[ERADIUS] = &DumpCustomVTK::pack_eradius;
      vtype[ERADIUS] = DOUBLE;
      name[ERADIUS] = arg[iarg];
    } else if (strcmp(arg[iarg],"ervel") == 0) {
      if (!atom->ervel_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[ERVEL] = &DumpCustomVTK::pack_ervel;
      vtype[ERVEL] = DOUBLE;
      name[ERVEL] = arg[iarg];
    } else if (strcmp(arg[iarg],"erforce") == 0) {
      if (!atom->erforce_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[ERFORCE] = &DumpCustomVTK::pack_erforce;
      vtype[ERFORCE] = DOUBLE;
      name[ERFORCE] = arg[iarg];
    } else if (strcmp(arg[iarg],"thread") == 0) {
      pack_choice[THREAD] = &DumpCustomVTK::pack_thread;
      vtype[THREAD] = INT;
      name[THREAD] = arg[iarg];

// superquadric start
    } else if (strcmp(arg[iarg],"shapex") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[SHAPEX] = &DumpCustomVTK::pack_shapex;
      vtype[SHAPEX] = DOUBLE;
      name[SHAPEX] = arg[iarg];
    } else if (strcmp(arg[iarg],"shapey") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[SHAPEY] = &DumpCustomVTK::pack_shapey;
      vtype[SHAPEY] = DOUBLE;
      name[SHAPEY] = arg[iarg];
    } else if (strcmp(arg[iarg],"shapez") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[SHAPEZ] = &DumpCustomVTK::pack_shapez;
      vtype[SHAPEZ] = DOUBLE;
      name[SHAPEZ] = arg[iarg];
    } else if (strcmp(arg[iarg],"quat1") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[QUAT1] = &DumpCustomVTK::pack_quat1;
      vtype[QUAT1] = DOUBLE;
      name[QUAT1] = arg[iarg];
    } else if (strcmp(arg[iarg],"quat2") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[QUAT2] = &DumpCustomVTK::pack_quat2;
      vtype[QUAT2] = DOUBLE;
      name[QUAT2] = arg[iarg];
    } else if (strcmp(arg[iarg],"quat3") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[QUAT3] = &DumpCustomVTK::pack_quat3;
      vtype[QUAT3] = DOUBLE;
      name[QUAT3] = arg[iarg];
    } else if (strcmp(arg[iarg],"quat4") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[QUAT4] = &DumpCustomVTK::pack_quat4;
      vtype[QUAT4] = DOUBLE;
      name[QUAT4] = arg[iarg];
    } else if (strcmp(arg[iarg],"blockiness1") == 0 || strcmp(arg[iarg],"roundness1") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      if(strcmp(arg[iarg],"roundness1") == 0)
        error->warning(FLERR,"Keyword 'roundness1' will be deprecated in future, please use 'blockiness1' istead");
      pack_choice[BLOCKINESS1] = &DumpCustomVTK::pack_blockiness1;
      vtype[BLOCKINESS1] = DOUBLE;
      name[BLOCKINESS1] = arg[iarg];
    } else if (strcmp(arg[iarg],"blockiness2") == 0 || strcmp(arg[iarg],"roundness2") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      if(strcmp(arg[iarg],"roundness2") == 0)
        error->warning(FLERR,"Keyword 'roundness2' will be deprecated in future, please use 'blockiness2' istead");
      pack_choice[BLOCKINESS2] = &DumpCustomVTK::pack_blockiness2;
      vtype[BLOCKINESS2] = DOUBLE;
      name[BLOCKINESS2] = arg[iarg];
    } else if (strcmp(arg[iarg],"inertiax") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[INERTIAX] = &DumpCustomVTK::pack_inertiax;
      vtype[INERTIAX] = DOUBLE;
      name[INERTIAX] = arg[iarg];
    } else if (strcmp(arg[iarg],"inertiay") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[INERTIAY] = &DumpCustomVTK::pack_inertiay;
      vtype[INERTIAY] = DOUBLE;
      name[INERTIAY] = arg[iarg];
    } else if (strcmp(arg[iarg],"inertiaz") == 0) {
      if (!atom->superquadric_flag)
        error->all(FLERR,"Dumping an atom quantity that isn't allocated");
      pack_choice[INERTIAZ] = &DumpCustomVTK::pack_inertiaz;
      vtype[INERTIAZ] = DOUBLE;
      name[INERTIAZ] = arg[iarg];
// superquadric end

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is int between []

    } else if (strncmp(arg[iarg],"c_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpCustomVTK::pack_compute;
      vtype[ATTRIBUTES+i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump custom/vtk command");
        argindex[ATTRIBUTES+i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+i] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom/vtk compute ID");
      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,"Dump custom/vtk compute does not compute per-atom info");
      if (argindex[ATTRIBUTES+i] == 0 && modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump custom/vtk compute does not calculate per-atom vector");
      if (argindex[ATTRIBUTES+i] > 0 && modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,\
                   "Dump custom/vtk compute does not calculate per-atom array");
      if (argindex[ATTRIBUTES+i] > 0 &&
          argindex[ATTRIBUTES+i] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump custom/vtk compute vector is accessed out-of-range");

      field2index[ATTRIBUTES+i] = add_compute(suffix);
      name[ATTRIBUTES+i] = arg[iarg];
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[iarg],"f_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpCustomVTK::pack_fix;
      vtype[ATTRIBUTES+i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump custom/vtk command");
        argindex[ATTRIBUTES+i] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+i] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom/vtk fix ID");
      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump custom/vtk fix does not compute per-atom info");
      if (argindex[ATTRIBUTES+i] == 0 && modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump custom/vtk fix does not compute per-atom vector");
      if (argindex[ATTRIBUTES+i] > 0 && modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump custom/vtk fix does not compute per-atom array");
      if (argindex[ATTRIBUTES+i] > 0 &&
          argindex[ATTRIBUTES+i] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump custom/vtk fix vector is accessed out-of-range");

      field2index[ATTRIBUTES+i] = add_fix(suffix);
      name[ATTRIBUTES+i] = arg[iarg];
      delete [] suffix;

    // variable value = v_name

    } else if (strncmp(arg[iarg],"v_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpCustomVTK::pack_variable;
      vtype[ATTRIBUTES+i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      argindex[ATTRIBUTES+i] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump custom/vtk variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump custom/vtk variable is not atom-style variable");

      field2index[ATTRIBUTES+i] = add_variable(suffix);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    // custom per-atom floating point value = d_ID

    } else if (strncmp(arg[iarg],"d_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpCustomVTK::pack_custom;
      vtype[ATTRIBUTES+i] = DOUBLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[ATTRIBUTES+i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 1)
        error->all(FLERR,"Custom per-atom property ID is not floating point");

      field2index[ATTRIBUTES+i] = add_custom(suffix,1);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    // custom per-atom integer value = i_ID

    } else if (strncmp(arg[iarg],"i_",2) == 0) {
      pack_choice[ATTRIBUTES+i] = &DumpCustomVTK::pack_custom;
      vtype[ATTRIBUTES+i] = INT;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);
      argindex[ATTRIBUTES+i] = 0;

      int tmp = -1;
      n = atom->find_custom(suffix,tmp);
      if (n < 0)
        error->all(FLERR,"Could not find custom per-atom property ID");

      if (tmp != 0)
        error->all(FLERR,"Custom per-atom property ID is not integer");

      field2index[ATTRIBUTES+i] = add_custom(suffix,0);
      name[ATTRIBUTES+i] = suffix;
      delete [] suffix;

    } else return iarg;
  }

  identify_vectors();

  return narg;
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::identify_vectors()
{
  // detect vectors
  vector_set.insert(X); // required

  int vector3_starts[] = {XS, XU, XSU, IX, VX, FX, MUX, OMEGAX, ANGMOMX, TQX, SHAPEX, INERTIAX};
  int num_vector3_starts = sizeof(vector3_starts) / sizeof(int);

  for (int v3s = 0; v3s < num_vector3_starts; v3s++) {
    if(name.count(vector3_starts[v3s]  ) &&
       name.count(vector3_starts[v3s]+1) &&
       name.count(vector3_starts[v3s]+2) )
    {
      std::string vectorName = name[vector3_starts[v3s]];
      vectorName.erase(vectorName.find_first_of('x'));
      name[vector3_starts[v3s]] = vectorName;
      vector_set.insert(vector3_starts[v3s]);
    }
  }

  // compute and fix vectors
  for (std::map<int,std::string>::iterator it=name.begin(); it!=name.end(); ++it) {
    if (it->first < ATTRIBUTES) // neither fix nor compute
      continue;

    if(argindex[it->first] == 0) // single value
      continue;

    // assume components are grouped together and in correct order
    if(name.count(it->first + 1) && name.count(it->first + 2) ) { // more attributes?
      if(it->second.compare(0,it->second.length()-3,name[it->first + 1],0,it->second.length()-3) == 0  && // same attributes?
         it->second.compare(0,it->second.length()-3,name[it->first + 2],0,it->second.length()-3) == 0 )
      {
        it->second.erase(it->second.length()-1);
        std::ostringstream oss;
        oss << "-" << argindex[it->first+2] << "]";
        it->second += oss.str();
        vector_set.insert(it->first);
        ++it; ++it;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   add Compute to list of Compute objects used by dump
   return index of where this Compute is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustomVTK::add_compute(char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,id_compute[icompute]) == 0) break;
  if (icompute < ncompute) return icompute;

  id_compute = (char **)
    memory->srealloc(id_compute,(ncompute+1)*sizeof(char *),"dump:id_compute");
  delete [] compute;
  compute = new Compute*[ncompute+1];

  int n = strlen(id) + 1;
  id_compute[ncompute] = new char[n];
  strcpy(id_compute[ncompute],id);
  ncompute++;
  return ncompute-1;
}

/* ----------------------------------------------------------------------
   add Fix to list of Fix objects used by dump
   return index of where this Fix is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustomVTK::add_fix(char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,id_fix[ifix]) == 0) break;
  if (ifix < nfix) return ifix;

  id_fix = (char **)
    memory->srealloc(id_fix,(nfix+1)*sizeof(char *),"dump:id_fix");
  delete [] fix;
  fix = new Fix*[nfix+1];

  int n = strlen(id) + 1;
  id_fix[nfix] = new char[n];
  strcpy(id_fix[nfix],id);
  nfix++;
  return nfix-1;
}

/* ----------------------------------------------------------------------
   add Variable to list of Variables used by dump
   return index of where this Variable is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustomVTK::add_variable(char *id)
{
  int ivariable;
  for (ivariable = 0; ivariable < nvariable; ivariable++)
    if (strcmp(id,id_variable[ivariable]) == 0) break;
  if (ivariable < nvariable) return ivariable;

  id_variable = (char **)
    memory->srealloc(id_variable,(nvariable+1)*sizeof(char *),
                     "dump:id_variable");
  delete [] variable;
  variable = new int[nvariable+1];
  delete [] vbuf;
  vbuf = new double*[nvariable+1];
  for (int i = 0; i <= nvariable; i++) vbuf[i] = NULL;

  int n = strlen(id) + 1;
  id_variable[nvariable] = new char[n];
  strcpy(id_variable[nvariable],id);
  nvariable++;
  return nvariable-1;
}

/* ----------------------------------------------------------------------
   add custom atom property to list used by dump
   return index of where this property is in list
   if already in list, do not add, just return index, else add to list
------------------------------------------------------------------------- */

int DumpCustomVTK::add_custom(char *id, int flag)
{
  int icustom;
  for (icustom = 0; icustom < ncustom; icustom++)
    if ((strcmp(id,id_custom[icustom]) == 0)
        && (flag == flag_custom[icustom])) break;
  if (icustom < ncustom) return icustom;

  id_custom = (char **)
    memory->srealloc(id_custom,(ncustom+1)*sizeof(char *),"dump:id_custom");
  flag_custom = (int *)
    memory->srealloc(flag_custom,(ncustom+1)*sizeof(int),"dump:flag_custom");

  int n = strlen(id) + 1;
  id_custom[ncustom] = new char[n];
  strcpy(id_custom[ncustom],id);
  flag_custom[ncustom] = flag;

  ncustom++;
  return ncustom-1;
}

/* ---------------------------------------------------------------------- */

int DumpCustomVTK::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"region") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) iregion = -1;
    else {
      iregion = domain->find_region(arg[1]);
      if (iregion == -1)
        error->all(FLERR,"Dump_modify region ID does not exist");
      delete [] idregion;
      int n = strlen(arg[1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[1]);
    }
    return 2;
  }

  if (strcmp(arg[0],"label") == 0) { //NP modified C.K.
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [label]");
     delete [] label;
     int n = strlen(arg[1]) + 1;
     label = new char[n];
     strcpy(label,arg[1]);
     return 2;
   }

  if (strcmp(arg[0],"binary") == 0) {
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [binary]");
     if (strcmp(arg[1],"yes") == 0) binary = 1;
     else if (strcmp(arg[1],"no") == 0) binary = 0;
     else error->all(FLERR,"Illegal dump_modify command [binary]");
     return 2;
  }

  if (strcmp(arg[0],"domainfile") == 0) {
     if (narg < 2) error->all(FLERR,"Illegal dump_modify command [domainfile]");
     if (strcmp(arg[1],"yes") == 0) write_domain = 1;
     else if (strcmp(arg[1],"no") == 0) write_domain = 0;
     else error->all(FLERR,"Illegal dump_modify command [domainfile]");
     return 2;
  }

  if (strcmp(arg[0],"element") == 0) {
    if (narg < ntypes+1)
      error->all(FLERR,"Dump modify: number of element names do not match atom types");

    if (typenames) {
      for (int i = 1; i <= ntypes; i++) delete [] typenames[i];
      delete [] typenames;
      typenames = NULL;
    }

    typenames = new char*[ntypes+1];
    for (int itype = 1; itype <= ntypes; itype++) {
      int n = strlen(arg[itype]) + 1;
      typenames[itype] = new char[n];
      strcpy(typenames[itype],arg[itype]);
    }
    return ntypes+1;
  }

  if (strcmp(arg[0],"thresh") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"none") == 0) {
      if (nthresh) {
        memory->destroy(thresh_array);
        memory->destroy(thresh_op);
        memory->destroy(thresh_value);
        thresh_array = NULL;
        thresh_op = NULL;
        thresh_value = NULL;
      }
      nthresh = 0;
      return 2;
    }

    if (narg < 4) error->all(FLERR,"Illegal dump_modify command");

    // grow threshold arrays

    memory->grow(thresh_array,nthresh+1,"dump:thresh_array");
    memory->grow(thresh_op,(nthresh+1),"dump:thresh_op");
    memory->grow(thresh_value,(nthresh+1),"dump:thresh_value");

    // set attribute type of threshold
    // customize by adding to if statement

    if (strcmp(arg[1],"id") == 0) thresh_array[nthresh] = ID;
    else if (strcmp(arg[1],"mol") == 0) thresh_array[nthresh] = MOL;
    else if (strcmp(arg[1],"type") == 0) thresh_array[nthresh] = TYPE;
    else if (strcmp(arg[1],"mass") == 0) thresh_array[nthresh] = MASS;

    else if (strcmp(arg[1],"x") == 0) thresh_array[nthresh] = X;
    else if (strcmp(arg[1],"y") == 0) thresh_array[nthresh] = Y;
    else if (strcmp(arg[1],"z") == 0) thresh_array[nthresh] = Z;

    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XS;
    else if (strcmp(arg[1],"xs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSTRI;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YS;
    else if (strcmp(arg[1],"ys") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSTRI;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZS;
    else if (strcmp(arg[1],"zs") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSTRI;

    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XU;
    else if (strcmp(arg[1],"xu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XUTRI;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YU;
    else if (strcmp(arg[1],"yu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YUTRI;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZU;
    else if (strcmp(arg[1],"zu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZUTRI;

    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = XSU;
    else if (strcmp(arg[1],"xsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = XSUTRI;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = YSU;
    else if (strcmp(arg[1],"ysu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = YSUTRI;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 0)
      thresh_array[nthresh] = ZSU;
    else if (strcmp(arg[1],"zsu") == 0 && domain->triclinic == 1)
      thresh_array[nthresh] = ZSUTRI;

    else if (strcmp(arg[1],"ix") == 0) thresh_array[nthresh] = IX;
    else if (strcmp(arg[1],"iy") == 0) thresh_array[nthresh] = IY;
    else if (strcmp(arg[1],"iz") == 0) thresh_array[nthresh] = IZ;
    else if (strcmp(arg[1],"vx") == 0) thresh_array[nthresh] = VX;
    else if (strcmp(arg[1],"vy") == 0) thresh_array[nthresh] = VY;
    else if (strcmp(arg[1],"vz") == 0) thresh_array[nthresh] = VZ;
    else if (strcmp(arg[1],"fx") == 0) thresh_array[nthresh] = FX;
    else if (strcmp(arg[1],"fy") == 0) thresh_array[nthresh] = FY;
    else if (strcmp(arg[1],"fz") == 0) thresh_array[nthresh] = FZ;

    else if (strcmp(arg[1],"q") == 0) thresh_array[nthresh] = Q;
    else if (strcmp(arg[1],"density") == 0) thresh_array[nthresh] = DENSITY; //NP modified C.K.
    else if (strcmp(arg[1],"p") == 0) thresh_array[nthresh] = P; //NP modified C.K.
    else if (strcmp(arg[1],"rho") == 0) thresh_array[nthresh] = RHO; //NP modified C.K.
    else if (strcmp(arg[1],"thread") == 0) thresh_array[nthresh] = THREAD; //NP modified R.B.
    else if (strcmp(arg[1],"mux") == 0) thresh_array[nthresh] = MUX;
    else if (strcmp(arg[1],"muy") == 0) thresh_array[nthresh] = MUY;
    else if (strcmp(arg[1],"muz") == 0) thresh_array[nthresh] = MUZ;
    else if (strcmp(arg[1],"mu") == 0) thresh_array[nthresh] = MU;

    else if (strcmp(arg[1],"radius") == 0) thresh_array[nthresh] = RADIUS;
    else if (strcmp(arg[1],"diameter") == 0) thresh_array[nthresh] = DIAMETER;
    else if (strcmp(arg[1],"omegax") == 0) thresh_array[nthresh] = OMEGAX;
    else if (strcmp(arg[1],"omegay") == 0) thresh_array[nthresh] = OMEGAY;
    else if (strcmp(arg[1],"omegaz") == 0) thresh_array[nthresh] = OMEGAZ;
    else if (strcmp(arg[1],"angmomx") == 0) thresh_array[nthresh] = ANGMOMX;
    else if (strcmp(arg[1],"angmomy") == 0) thresh_array[nthresh] = ANGMOMY;
    else if (strcmp(arg[1],"angmomz") == 0) thresh_array[nthresh] = ANGMOMZ;
    else if (strcmp(arg[1],"tqx") == 0) thresh_array[nthresh] = TQX;
    else if (strcmp(arg[1],"tqy") == 0) thresh_array[nthresh] = TQY;
    else if (strcmp(arg[1],"tqz") == 0) thresh_array[nthresh] = TQZ;

    else if (strcmp(arg[1],"spin") == 0) thresh_array[nthresh] = SPIN;
    else if (strcmp(arg[1],"eradius") == 0) thresh_array[nthresh] = ERADIUS;
    else if (strcmp(arg[1],"ervel") == 0) thresh_array[nthresh] = ERVEL;
    else if (strcmp(arg[1],"erforce") == 0) thresh_array[nthresh] = ERFORCE;
// superquadric start
    else if (strcmp(arg[1],"shapex") == 0) thresh_array[nthresh] = SHAPEX;
    else if (strcmp(arg[1],"shapey") == 0) thresh_array[nthresh] = SHAPEY;
    else if (strcmp(arg[1],"shapez") == 0) thresh_array[nthresh] = SHAPEZ;
    else if (strcmp(arg[1],"quat1") == 0) thresh_array[nthresh] = QUAT1;
    else if (strcmp(arg[1],"quat2") == 0) thresh_array[nthresh] = QUAT2;
    else if (strcmp(arg[1],"quat3") == 0) thresh_array[nthresh] = QUAT3;
    else if (strcmp(arg[1],"quat4") == 0) thresh_array[nthresh] = QUAT4;
    else if (strcmp(arg[1],"blockiness1") == 0) thresh_array[nthresh] = BLOCKINESS1;
    else if (strcmp(arg[1],"blockiness2") == 0) thresh_array[nthresh] = BLOCKINESS2;
    else if (strcmp(arg[1],"inertiax") == 0) thresh_array[nthresh] = INERTIAX;
    else if (strcmp(arg[1],"inertiay") == 0) thresh_array[nthresh] = INERTIAY;
    else if (strcmp(arg[1],"inertiaz") == 0) thresh_array[nthresh] = INERTIAZ;
// superquadric end

    // compute value = c_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    else if (strncmp(arg[1],"c_",2) == 0) {
      thresh_array[nthresh] = COMPUTE;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = modify->find_compute(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify compute ID");

      if (modify->compute[n]->peratom_flag == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom info");
      if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
          modify->compute[n]->size_peratom_cols > 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom vector");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          modify->compute[n]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Dump modify compute ID does not compute per-atom array");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          argindex[ATTRIBUTES+nfield+nthresh] > modify->compute[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify compute ID vector is not large enough");

      field2index[ATTRIBUTES+nfield+nthresh] = add_compute(suffix);
      delete [] suffix;

    // fix value = f_ID
    // if no trailing [], then arg is set to 0, else arg is between []

    } else if (strncmp(arg[1],"f_",2) == 0) {
      thresh_array[nthresh] = FIX;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Invalid attribute in dump modify command");
        argindex[ATTRIBUTES+nfield+nthresh] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = modify->find_fix(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify fix ID");

      if (modify->fix[n]->peratom_flag == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom info");
      if (argindex[ATTRIBUTES+nfield+nthresh] == 0 &&
          modify->fix[n]->size_peratom_cols > 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom vector");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          modify->fix[n]->size_peratom_cols == 0)
        error->all(FLERR,"Dump modify fix ID does not compute per-atom array");
      if (argindex[ATTRIBUTES+nfield+nthresh] > 0 &&
          argindex[ATTRIBUTES+nfield+nthresh] > modify->fix[n]->size_peratom_cols)
        error->all(FLERR,"Dump modify fix ID vector is not large enough");

      field2index[ATTRIBUTES+nfield+nthresh] = add_fix(suffix);
      delete [] suffix;

    // variable value = v_ID

    } else if (strncmp(arg[1],"v_",2) == 0) {
      thresh_array[nthresh] = VARIABLE;
      int n = strlen(arg[1]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[1][2]);

      argindex[ATTRIBUTES+nfield+nthresh] = 0;

      n = input->variable->find(suffix);
      if (n < 0) error->all(FLERR,"Could not find dump modify variable name");
      if (input->variable->atomstyle(n) == 0)
        error->all(FLERR,"Dump modify variable is not atom-style variable");

      field2index[ATTRIBUTES+nfield+nthresh] = add_variable(suffix);
      delete [] suffix;

    } else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set operation type of threshold

    if (strcmp(arg[2],"<") == 0) thresh_op[nthresh] = LT;
    else if (strcmp(arg[2],"<=") == 0) thresh_op[nthresh] = LE;
    else if (strcmp(arg[2],">") == 0) thresh_op[nthresh] = GT;
    else if (strcmp(arg[2],">=") == 0) thresh_op[nthresh] = GE;
    else if (strcmp(arg[2],"==") == 0) thresh_op[nthresh] = EQ;
    else if (strcmp(arg[2],"!=") == 0) thresh_op[nthresh] = NEQ;
    else error->all(FLERR,"Invalid dump_modify threshold operator");

    // set threshold value

    thresh_value[nthresh] = force->numeric(FLERR,arg[3]);

    nthresh++;
    return 4;
  }

  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory in buf, choose, variable arrays
------------------------------------------------------------------------- */

bigint DumpCustomVTK::memory_usage()
{
  bigint bytes = Dump::memory_usage();
  bytes += memory->usage(choose,maxlocal);
  bytes += memory->usage(dchoose,maxlocal);
  bytes += memory->usage(clist,maxlocal);
  bytes += memory->usage(vbuf,nvariable,maxlocal);
  return bytes;
}

/* ----------------------------------------------------------------------
   extraction of Compute, Fix, Variable results
------------------------------------------------------------------------- */

void DumpCustomVTK::pack_compute(int n)
{
  double *vector = compute[field2index[current_pack_choice_key]]->vector_atom;
  double **array = compute[field2index[current_pack_choice_key]]->array_atom;
  int index = argindex[current_pack_choice_key];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_fix(int n)
{
  double *vector = fix[field2index[current_pack_choice_key]]->vector_atom;
  double **array = fix[field2index[current_pack_choice_key]]->array_atom;
  int index = argindex[current_pack_choice_key];

  if (index == 0) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = vector[clist[i]];
      n += size_one;
    }
  } else {
    index--;
    for (int i = 0; i < nchoose; i++) {
      buf[n] = array[clist[i]][index];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_variable(int n)
{
  double *vector = vbuf[field2index[current_pack_choice_key]];

  for (int i = 0; i < nchoose; i++) {
    buf[n] = vector[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_custom(int n)
{

  int index = field2index[n];

  if (flag_custom[index] == 0) { // integer
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    int *ivector = atom->ivector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = ivector[clist[i]];
      n += size_one;
    }
  } else if (flag_custom[index] == 1) { // double
    int iwhich,tmp;
    iwhich = atom->find_custom(id_custom[index],tmp);

    double *dvector = atom->dvector[iwhich];
    for (int i = 0; i < nchoose; i++) {
      buf[n] = dvector[clist[i]];
      n += size_one;
    }
  }
}

/* ----------------------------------------------------------------------
   one method for every attribute dump custom/vtk can output
   the atom property is packed into buf starting at n with stride size_one
   customize a new attribute by adding a method
------------------------------------------------------------------------- */

void DumpCustomVTK::pack_id(int n)
{
  int *tag = atom->tag;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = tag[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_molecule(int n)
{
  int *molecule = atom->molecule;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = molecule[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_type(int n)
{
  int *type = atom->type;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = type[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_mass(int n)
{
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;

  if (rmass) {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = rmass[clist[i]];
      n += size_one;
    }
  } else {
    for (int i = 0; i < nchoose; i++) {
      buf[n] = mass[type[clist[i]]];
      n += size_one;
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_x(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_y(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_z(int n)
{
  double **x = atom->x;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = x[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xs(int n)
{
  double **x = atom->x;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][0] - boxxlo) * invxprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ys(int n)
{
  double **x = atom->x;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][1] - boxylo) * invyprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zs(int n)
{
  double **x = atom->x;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (x[clist[i]][2] - boxzlo) * invzprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xs_triclinic(int n)
{
  int j;
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
      h_inv[4]*(x[j][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ys_triclinic(int n)
{
  int j;
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zs_triclinic(int n)
{
  double **x = atom->x;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = h_inv[2]*(x[clist[i]][2]-boxlo[2]);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double xprd = domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][0] + ((image[j] & IMGMASK) - IMGMAX) * xprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_yu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double yprd = domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][1] + ((image[j] >> IMGBITS & IMGMASK) - IMGMAX) * yprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double zprd = domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = x[j][2] + ((image[j] >> IMG2BITS) - IMGMAX) * zprd;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *h = domain->h;
  int xbox,ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    xbox = (image[j] & IMGMASK) - IMGMAX;
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_yu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *h = domain->h;
  int ybox,zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    ybox = (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][1] + h[1]*ybox + h[3]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *h = domain->h;
  int zbox;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    zbox = (image[j] >> IMG2BITS) - IMGMAX;
    buf[n] = x[j][2] + h[2]*zbox;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xsu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double boxxlo = domain->boxlo[0];
  double invxprd = 1.0/domain->xprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][0] - boxxlo) * invxprd + (image[j] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ysu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double boxylo = domain->boxlo[1];
  double invyprd = 1.0/domain->yprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][1] - boxylo) * invyprd + (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zsu(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double boxzlo = domain->boxlo[2];
  double invzprd = 1.0/domain->zprd;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = (x[j][2] - boxzlo) * invzprd + (image[j] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_xsu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[0]*(x[j][0]-boxlo[0]) + h_inv[5]*(x[j][1]-boxlo[1]) +
      h_inv[4]*(x[j][2]-boxlo[2]) + (image[j] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ysu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[1]*(x[j][1]-boxlo[1]) + h_inv[3]*(x[j][2]-boxlo[2]) +
      (image[j] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_zsu_triclinic(int n)
{
  int j;
  double **x = atom->x;
  tagint *image = atom->image;

  double *boxlo = domain->boxlo;
  double *h_inv = domain->h_inv;

  for (int i = 0; i < nchoose; i++) {
    j = clist[i];
    buf[n] = h_inv[2]*(x[j][2]-boxlo[2]) + (image[j] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ix(int n)
{
  tagint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_iy(int n)
{
  tagint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] >> IMGBITS & IMGMASK) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_iz(int n)
{
  tagint *image = atom->image;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = (image[clist[i]] >> IMG2BITS) - IMGMAX;
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_vx(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_vy(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_vz(int n)
{
  double **v = atom->v;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = v[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_fx(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_fy(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_fz(int n)
{
  double **f = atom->f;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = f[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_q(int n)
{
  double *q = atom->q;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = q[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_density(int n) //NP modified C.K.
{
  double *density = atom->density;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = density[clist[i]];
    n += size_one;
  }

}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_p(int n) //NP modified C.K.
{
  double *p = atom->p;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = p[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_rho(int n) //NP modified C.K.
{
  double *rho = atom->rho;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = rho[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_mux(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_muy(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_muz(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_mu(int n)
{
  double **mu = atom->mu;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = mu[clist[i]][3];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_radius(int n)
{
  double *radius = atom->radius;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = radius[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_diameter(int n)
{
  double *radius = atom->radius;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = 2.0*radius[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_omegax(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_omegay(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_omegaz(int n)
{
  double **omega = atom->omega;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = omega[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_angmomx(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_angmomy(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_angmomz(int n)
{
  double **angmom = atom->angmom;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = angmom[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_tqx(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_tqy(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_tqz(int n)
{
  double **torque = atom->torque;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = torque[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_spin(int n)
{
  int *spin = atom->spin;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = spin[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_eradius(int n)
{
  double *eradius = atom->eradius;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = eradius[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_ervel(int n)
{
  double *ervel = atom->ervel;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = ervel[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_erforce(int n)
{
  double *erforce = atom->erforce;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = erforce[clist[i]];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_shapex(int n)
{
  double **shape = atom->shape;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = shape[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_shapey(int n)
{
  double **shape = atom->shape;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = shape[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_shapez(int n)
{
  double **shape = atom->shape;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = shape[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_quat1(int n)
{
  double **quaternion = atom->quaternion;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = quaternion[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_quat2(int n)
{
  double **quaternion = atom->quaternion;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = quaternion[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_quat3(int n)
{
  double **quaternion = atom->quaternion;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = quaternion[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_quat4(int n)
{
  double **quaternion = atom->quaternion;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = quaternion[clist[i]][3];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_blockiness1(int n)
{
  double **blockiness = atom->blockiness;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = blockiness[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_blockiness2(int n)
{
  double **blockiness = atom->blockiness;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = blockiness[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_inertiax(int n)
{
  double **inertia = atom->inertia;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = inertia[clist[i]][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_inertiay(int n)
{
  double **inertia = atom->inertia;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = inertia[clist[i]][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_inertiaz(int n)
{
  double **inertia = atom->inertia;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = inertia[clist[i]][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpCustomVTK::pack_thread(int n)
{
  int *thread = atom->thread;

  for (int i = 0; i < nchoose; i++) {
    buf[n] = comm->me * comm->nprocs + thread[clist[i]];
    n += size_one;
  }
}
#endif
