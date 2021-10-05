/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <stdlib.h>
#include "atom_vec_bond_gran.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include <string.h>
#include "memory.h"
#include "error.h"
#include "comm.h"
#include "update.h"

using namespace LAMMPS_NS;

#define DELTA 10000
#define BOND_GRAN_USE_MOLECULE_ID

/* ---------------------------------------------------------------------- */

AtomVecBondGran::AtomVecBondGran(LAMMPS *lmp) :
  AtomVec(lmp)
{
  molecular = 1;
  bonds_allow = 1;
  mass_type = 1;

  comm_x_only = comm_f_only = 1;
  size_forward = 3;
  size_reverse = 3;

  size_border = 7;
  size_velocity = 3;
#ifdef BOND_GRAN_USE_MOLECULE_ID
  size_data_atom = 6; // with molecule ID
#else
  size_data_atom = 5; // without molecule ID
#endif
  size_data_vel = 4;
  xcol_data = 4;

  atom->molecule_flag = 1;

  fbpg = NULL;
}

/* ----------------------------------------------------------------------
   no additional args by default
------------------------------------------------------------------------- */

void AtomVecBondGran::settings(int narg, char **arg)
{
  if (narg == 0) return; //in case of restart no arguments are given, instead nbondtypes and bond_per_atom are defined by read_restart_settings

  if (narg != 4) error->all(FLERR,"Invalid atom_style bond/gran command,expecting exactly 4 arguments");

  if(strcmp(arg[0],"n_bondtypes"))
    error->all(FLERR,"Illegal atom_style bond/gran command, expecting 'n_bondtypes'");

  atom->nbondtypes = atoi(arg[1]);

  if(strcmp(arg[2],"bonds_per_atom"))
    error->all(FLERR,"Illegal atom_style bond/gran command, expecting 'bonds_per_atom'");

  atom->bond_per_atom = atoi(arg[3]);

  /*NL*/ if(screen) fprintf(screen,"atom->nbondtypes %d atom->bond_per_atom %d\n",atom->nbondtypes,atom->bond_per_atom);
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::init()
{

  if(fbpg == NULL)
  {
      char **fixarg = new char*[3];
      fixarg[0] = (char *) "BOND_PROPAGATE";
      fixarg[1] = (char *) "all";
      fixarg[2] = (char *) "bond/propagate/gran";
      modify->add_fix(3,fixarg);
      delete [] fixarg;
  }
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecBondGran::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;

  tag = memory->grow(atom->tag,nmax,"atom:tag");
  type = memory->grow(atom->type,nmax,"atom:type");
  mask = memory->grow(atom->mask,nmax,"atom:mask");
  image = memory->grow(atom->image,nmax,"atom:image");
  x = memory->grow(atom->x,nmax,3,"atom:x");
  v = memory->grow(atom->v,nmax,3,"atom:v");
  f = memory->grow(atom->f,nmax*comm->nthreads,3,"atom:f");

  x_mol = memory->grow(atom->x_mol,nmax,3,"atom:x_mol");    // modified A.N.
  v_mol = memory->grow(atom->v_mol,nmax,3,"atom:v_mol");    // modified A.N.

  molecule = memory->grow(atom->molecule,nmax,"atom:molecule");

  nspecial = memory->grow(atom->nspecial,nmax,3,"atom:nspecial");
  special = memory->grow(atom->special,nmax,atom->maxspecial,"atom:special");

  num_bond = memory->grow(atom->num_bond,nmax,"atom:num_bond");
  bond_type = memory->grow(atom->bond_type,nmax,atom->bond_per_atom,"atom:bond_type");
  bond_atom = memory->grow(atom->bond_atom,nmax,atom->bond_per_atom,"atom:bond_atom");

  /*NL*/ //if(screen) fprintf(screen,"grow nmax %d, atom->bond_per_atom %d, atom->n_bondhist %d\n",nmax,atom->bond_per_atom,atom->n_bondhist);

  if(0 == atom->bond_per_atom)
    error->all(FLERR,"Bonded particles need bond_per_atom > 0");

  if(atom->n_bondhist < 0)
    error->all(FLERR,"atom->n_bondhist < 0 suggests that 'bond_style gran' has not been called before 'read_restart' command! Please check that.");

  if(atom->n_bondhist)
  {
     bond_hist =
        memory->grow(atom->bond_hist,nmax,atom->bond_per_atom,atom->n_bondhist,"atom:bond_hist");
  }

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecBondGran::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  molecule = atom->molecule;
  nspecial = atom->nspecial; special = atom->special;
  num_bond = atom->num_bond; bond_type = atom->bond_type;
  bond_atom = atom->bond_atom;
  bond_hist = atom->bond_hist;
  x_mol = atom->x_mol; v_mol = atom->v_mol;     // modified A.N.
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::copy(int i, int j, int delflag)
{
  int k,l;

  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  // Modified A.N.
  x_mol[j][0] = x_mol[i][0];
  x_mol[j][1] = x_mol[i][1];
  x_mol[j][2] = x_mol[i][2];
  v_mol[j][0] = v_mol[i][0];
  v_mol[j][1] = v_mol[i][1];
  v_mol[j][2] = v_mol[i][2];

  molecule[j] = molecule[i];

  num_bond[j] = num_bond[i];
  for (k = 0; k < num_bond[j]; k++) {
    bond_type[j][k] = bond_type[i][k];
    bond_atom[j][k] = bond_atom[i][k];
  }

  if(atom->n_bondhist)
  {
      for (k = 0; k < num_bond[j]; k++)
         for (l =0; l < atom->n_bondhist; l++)
            bond_hist[j][k][l] = bond_hist[i][k][l];
  }

  nspecial[j][0] = nspecial[i][0];
  nspecial[j][1] = nspecial[i][1];
  nspecial[j][2] = nspecial[i][2];
  for (k = 0; k < nspecial[j][2]; k++) special[j][k] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j,delflag);
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_comm(int n, int *list, double *buf,
         int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_comm_vel(int n, int *list, double *buf,
             int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = x_mol[j][0];
        buf[m++] = x_mol[j][1];
        buf[m++] = x_mol[j][2];
        buf[m++] = v_mol[j][0];
        buf[m++] = v_mol[j][1];
        buf[m++] = v_mol[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = x_mol[j][0];
        buf[m++] = x_mol[j][1];
        buf[m++] = x_mol[j][2];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
          buf[m++] = v_mol[j][0];
          buf[m++] = v_mol[j][1];
          buf[m++] = v_mol[j][2];
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = v_mol[j][0];
          buf[m++] = v_mol[j][1];
          buf[m++] = v_mol[j][2];
        }
      }
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    x_mol[i][0] = buf[m++];
    x_mol[i][1] = buf[m++];
    x_mol[i][2] = buf[m++];
    v_mol[i][0] = buf[m++];
    v_mol[i][1] = buf[m++];
    v_mol[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    x_mol[i][0] = buf[m++];
    x_mol[i][1] = buf[m++];
    x_mol[i][2] = buf[m++];
    v_mol[i][0] = buf[m++];
    v_mol[i][1] = buf[m++];
    v_mol[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_border(int n, int *list, double *buf,
           int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_border_vel(int n, int *list, double *buf,
         int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz,dvx,dvy,dvz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = molecule[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
      buf[m++] = x_mol[j][0];
      buf[m++] = x_mol[j][1];
      buf[m++] = x_mol[j][2];
      buf[m++] = v_mol[j][0];
      buf[m++] = v_mol[j][1];
      buf[m++] = v_mol[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    if (!deform_vremap) {
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = tag[j];
        buf[m++] = type[j];
        buf[m++] = mask[j];
        buf[m++] = molecule[j];
        buf[m++] = v[j][0];
        buf[m++] = v[j][1];
        buf[m++] = v[j][2];
        buf[m++] = x_mol[j][0];
        buf[m++] = x_mol[j][1];
        buf[m++] = x_mol[j][2];
        buf[m++] = v_mol[j][0];
        buf[m++] = v_mol[j][1];
        buf[m++] = v_mol[j][2];
      }
    } else {
      dvx = pbc[0]*h_rate[0] + pbc[5]*h_rate[5] + pbc[4]*h_rate[4];
      dvy = pbc[1]*h_rate[1] + pbc[3]*h_rate[3];
      dvz = pbc[2]*h_rate[2];
      for (i = 0; i < n; i++) {
        j = list[i];
        buf[m++] = x[j][0] + dx;
        buf[m++] = x[j][1] + dy;
        buf[m++] = x[j][2] + dz;
        buf[m++] = tag[j];
        buf[m++] = type[j];
        buf[m++] = mask[j];
        buf[m++] = molecule[j];
        buf[m++] = x_mol[j][0];
        buf[m++] = x_mol[j][1];
        buf[m++] = x_mol[j][2];
        if (mask[i] & deform_groupbit) {
          buf[m++] = v[j][0] + dvx;
          buf[m++] = v[j][1] + dvy;
          buf[m++] = v[j][2] + dvz;
          buf[m++] = v_mol[j][0];
          buf[m++] = v_mol[j][1];
          buf[m++] = v_mol[j][2];
        } else {
          buf[m++] = v[j][0];
          buf[m++] = v[j][1];
          buf[m++] = v[j][2];
          buf[m++] = v_mol[j][0];
          buf[m++] = v_mol[j][1];
          buf[m++] = v_mol[j][2];
        }
      }
    }
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->pack_border(n,list,&buf[m]);

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_border_hybrid(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = molecule[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::pack_border_one(int i, double *buf)
{
  buf[0] = molecule[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    molecule[i] = static_cast<int> (buf[m++]);
    x_mol[i][0] = buf[m++];
    x_mol[i][1] = buf[m++];
    x_mol[i][2] = buf[m++];
    v_mol[i][0] = buf[m++];
    v_mol[i][1] = buf[m++];
    v_mol[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

void AtomVecBondGran::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    molecule[i] = static_cast<int> (buf[m++]);
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
    x_mol[i][0] = buf[m++];
    x_mol[i][1] = buf[m++];
    x_mol[i][2] = buf[m++];
    v_mol[i][0] = buf[m++];
    v_mol[i][1] = buf[m++];
    v_mol[i][2] = buf[m++];
  }

  if (atom->nextra_border)
    for (int iextra = 0; iextra < atom->nextra_border; iextra++)
      m += modify->fix[atom->extra_border[iextra]]->
        unpack_border(n,first,&buf[m]);
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::unpack_border_hybrid(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++)
    molecule[i] = static_cast<int> (buf[m++]);
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::unpack_border_one(int i, double *buf)
{
  molecule[i] = static_cast<int> (buf[0]);
  return 1;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecBondGran::pack_exchange(int i, double *buf)
{
  int k,l;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m] = 0.0;      // for valgrind
  *((tagint *) &buf[m++]) = image[i];

  buf[m++] = molecule[i];

  buf[m++] = x_mol[i][0];
  buf[m++] = x_mol[i][1];
  buf[m++] = x_mol[i][2];
  buf[m++] = v_mol[i][0];
  buf[m++] = v_mol[i][1];
  buf[m++] = v_mol[i][2];

  buf[m++] = num_bond[i];
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = bond_type[i][k];
    buf[m++] = bond_atom[i][k];
  }

  if(atom->n_bondhist)
  {
      for (k = 0; k < num_bond[i]; k++)
        for (l = 0; l < atom->n_bondhist; l++)
          buf[m++] = bond_hist[i][k][l];
  }

  buf[m++] = nspecial[i][0];
  buf[m++] = nspecial[i][1];
  buf[m++] = nspecial[i][2];
  for (k = 0; k < nspecial[i][2]; k++) buf[m++] = special[i][k];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecBondGran::unpack_exchange(double *buf)
{
  int k,l;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);

  molecule[nlocal] = static_cast<int> (buf[m++]);

  x_mol[nlocal][0] = buf[m++];
  x_mol[nlocal][1] = buf[m++];
  x_mol[nlocal][2] = buf[m++];
  v_mol[nlocal][0] = buf[m++];
  v_mol[nlocal][1] = buf[m++];
  v_mol[nlocal][2] = buf[m++];

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }

  if(atom->n_bondhist)
  {
      for (k = 0; k < num_bond[nlocal]; k++)
        for (l = 0; l < atom->n_bondhist; l++)
          bond_hist[nlocal][k][l] = buf[m++];
  }

  nspecial[nlocal][0] = static_cast<int> (buf[m++]);
  nspecial[nlocal][1] = static_cast<int> (buf[m++]);
  nspecial[nlocal][2] = static_cast<int> (buf[m++]);
  for (k = 0; k < nspecial[nlocal][2]; k++)
    special[nlocal][k] = static_cast<int> (buf[m++]);

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
  unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecBondGran::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 0;
  for (i = 0; i < nlocal; i++)
  {
    n += 13 + 2*num_bond[i];

    if(atom->n_bondhist) n += 1/*num_bondhist*/ + num_bond[i] * atom->n_bondhist/*bond_hist*/;
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
  n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecBondGran::pack_restart(int i, double *buf)
{
  int k,l;

  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m] = 0.0;      // for valgrind
  *((tagint *) &buf[m++]) = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  buf[m++] = molecule[i];

  buf[m++] = num_bond[i];
  for (k = 0; k < num_bond[i]; k++) {
    buf[m++] = MAX(bond_type[i][k],-bond_type[i][k]);
    buf[m++] = bond_atom[i][k];
  }

  if(atom->n_bondhist)
  {
      buf[m++] = atom->n_bondhist;
      for (k = 0; k < num_bond[i]; k++)
         for (l = 0; l < atom->n_bondhist; l++)
            buf[m++] = bond_hist[i][k][l];
  }

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecBondGran::unpack_restart(double *buf)
{
  int k,l;

  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      atom->extra = memory->grow(atom->extra,nmax,atom->nextra_store,"atom:extra");
}

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = *((tagint *) &buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  molecule[nlocal] = static_cast<int> (buf[m++]);

  num_bond[nlocal] = static_cast<int> (buf[m++]);
  for (k = 0; k < num_bond[nlocal]; k++) {
    bond_type[nlocal][k] = static_cast<int> (buf[m++]);
    bond_atom[nlocal][k] = static_cast<int> (buf[m++]);
  }

  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  if(atom->n_bondhist)
  {
      if(atom->n_bondhist != static_cast<int>(buf[m++]))
          error->all(FLERR,"Íncompatibel restart file: file was created using a bond model with a different number of history values");
      for (k = 0; k < num_bond[nlocal]; k++)
         for (l = 0; l < atom->n_bondhist; l++)
            atom->bond_hist[nlocal][k][l] = buf[m++];
  }

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecBondGran::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = ((tagint) IMGMAX << IMG2BITS) |
    ((tagint) IMGMAX << IMGBITS) | IMGMAX;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  molecule[nlocal] = 0;
  num_bond[nlocal] = 0;
  nspecial[nlocal][0] = nspecial[nlocal][1] = nspecial[nlocal][2] = 0;

  x_mol[nlocal][0] = 0;
  x_mol[nlocal][1] = 0;
  x_mol[nlocal][2] = 0;
  v_mol[nlocal][0] = 0;
  v_mol[nlocal][1] = 0;
  v_mol[nlocal][2] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecBondGran::data_atom(double *coord, tagint imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = atoi(values[0]);
  if (tag[nlocal] <= 0)
    error->one(FLERR,"Invalid atom ID in Atoms section of data file");

#ifdef BOND_GRAN_USE_MOLECULE_ID
  molecule[nlocal] = atoi(values[1]);

  type[nlocal] = atoi(values[2]);
#else
  molecule[nlocal] = 0;

  type[nlocal] = atoi(values[1]);
#endif
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one(FLERR,"Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;
  num_bond[nlocal] = 0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecBondGran::data_atom_hybrid(int nlocal, char **values)
{
#ifdef BOND_GRAN_USE_MOLECULE_ID
  molecule[nlocal] = atoi(values[0]);
#else
  molecule[nlocal] = 0;
#endif
  num_bond[nlocal] = 0;

  return 1;
}

/* ----------------------------------------------------------------------
   pack atom info for data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBondGran::pack_data(double **buf)
{
  error->all(FLERR,"needs revision");
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    buf[i][0] = tag[i];
    buf[i][1] = molecule[i];
    buf[i][2] = type[i];
    buf[i][3] = x[i][0];
    buf[i][4] = x[i][1];
    buf[i][5] = x[i][2];
    buf[i][6] = (image[i] & IMGMASK) - IMGMAX;
    buf[i][7] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    buf[i][8] = (image[i] >> IMG2BITS) - IMGMAX;
  }
}

/* ----------------------------------------------------------------------
   pack hybrid atom info for data file
------------------------------------------------------------------------- */

int AtomVecBondGran::pack_data_hybrid(int i, double *buf)
{
  error->all(FLERR,"needs revision");
  buf[0] = molecule[i];
  return 1;
}

/* ----------------------------------------------------------------------
   write atom info to data file including 3 image flags
------------------------------------------------------------------------- */

void AtomVecBondGran::write_data(FILE *fp, int n, double **buf)
{
  error->all(FLERR,"needs revision");
  for (int i = 0; i < n; i++)
    fprintf(fp,"%d %d %d %-1.16e %-1.16e %-1.16e %d %d %d \n",
            (int) buf[i][0],(int) buf[i][1],(int) buf[i][2],
            buf[i][3],buf[i][4],buf[i][5],
            (int) buf[i][6],(int) buf[i][7],(int) buf[i][8]);
}

/* ----------------------------------------------------------------------
   write hybrid atom info to data file
------------------------------------------------------------------------- */

int AtomVecBondGran::write_data_hybrid(FILE *fp, double *buf)
{
  error->all(FLERR,"needs revision");
  fprintf(fp," %d",(int) buf[0]);
  return 1;
}

void AtomVecBondGran::write_restart_settings(FILE *fp)
{
  fwrite(&atom->nbondtypes,sizeof(int),1,fp);
  fwrite(&atom->bond_per_atom,sizeof(int),1,fp);
}

void AtomVecBondGran::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&atom->nbondtypes,sizeof(int),1,fp);
    fread(&atom->bond_per_atom,sizeof(int),1,fp);
  }
  MPI_Bcast(&atom->nbondtypes,1,MPI_INT,0,world);
  MPI_Bcast(&atom->bond_per_atom,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

bigint AtomVecBondGran::memory_usage()
{
  bigint bytes = 0.0;

  if (atom->memcheck("tag")) bytes += memory->usage(tag,nmax);
  if (atom->memcheck("type")) bytes += memory->usage(type,nmax);
  if (atom->memcheck("mask")) bytes += memory->usage(mask,nmax);
  if (atom->memcheck("image")) bytes += memory->usage(image,nmax);
  if (atom->memcheck("x")) bytes += memory->usage(x,nmax,3);
  if (atom->memcheck("v")) bytes += memory->usage(v,nmax,3);
  if (atom->memcheck("f")) bytes += memory->usage(f,nmax*comm->nthreads,3);

  if (atom->memcheck("molecule")) bytes += memory->usage(molecule,nmax);
  if (atom->memcheck("nspecial")) bytes += memory->usage(nspecial,nmax,3);
  if (atom->memcheck("special"))  bytes += memory->usage(special,nmax,atom->maxspecial);

  if (atom->memcheck("num_bond")) bytes += memory->usage(num_bond,nmax);
  if (atom->memcheck("bond_type")) bytes += memory->usage(bond_type,nmax,atom->bond_per_atom);
  if (atom->memcheck("bond_atom")) bytes += memory->usage(bond_atom,nmax,atom->bond_per_atom);
  if (atom->memcheck("x_mol")) bytes += memory->usage(x_mol,nmax,3); // mod A.N
  if (atom->memcheck("v_mol")) bytes += memory->usage(v_mol,nmax,3); // mod A.N
// P.F. not too sure about atom->n_bondhist
  if (atom->n_bondhist) bytes += nmax*sizeof(int);

  return bytes;
}
