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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_nve_sphere_limit.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVESphereLimit::FixNVESphereLimit(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix nve/sphere command");

  time_integrate = 1;

  vlimit = omegalimit = 0.;

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/sphere/limit command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else error->all(FLERR,"Illegal fix nve/sphere/limit command");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"vlimit") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix nve/sphere/limit command");
      vlimit = atof(arg[iarg+1]);
      if(vlimit <= 0.)
         error->all(FLERR,"Illegal fix nve/sphere/limit command, vlimit > 0 required");
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"omegalimit") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix nve/sphere/limit command");
      omegalimit = atof(arg[iarg+1]);
      if(omegalimit <= 0.)
         error->all(FLERR,"Illegal fix nve/sphere/limit command, omegalimit > 0 required");
      iarg += 2;
    }
    error->all(FLERR,"Illegal fix nve/sphere/limit command");
  }

  vlimitsq = vlimit * vlimit;
  omegalimitsq = omegalimit * omegalimit;

  if(vlimitsq <= 0. || omegalimitsq <= 0.)
    error->all(FLERR,"Illegal fix nve/sphere/limit command, vlimit > 0 and omegalimit > 0 required");

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nve/sphere requires atom style sphere");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::init()
{
  FixNVE::init();

  // check that all particles are finite-size spheres
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nve/sphere requires extended particles");
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::initial_integrate(int vflag)
{
  int itype;
  double dtfm,dtirotate,msq,scale,vsq,osq;
  double g[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate; //NP modified GM
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
          vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > vlimitsq) {
          ncount++;
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
      }
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
          osq = omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + omega[i][2]*omega[i][2];
      if (osq > omegalimitsq) {
          ncount++;
          scale = sqrt(omegalimitsq/osq);
          omega[i][0] *= scale;
          omega[i][1] *= scale;
          omega[i][2] *= scale;
      }
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (mu[i][3] > 0.0) {
          g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
          g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
          g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
          msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
          scale = mu[i][3]/sqrt(msq);
          mu[i][0] = g[0]*scale;
          mu[i][1] = g[1]*scale;
          mu[i][2] = g[2]*scale;
        }
      }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::final_integrate()
{
  int itype;
  double dtfm,dtirotate,vsq,osq,scale;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate;
  if (domain->dimension == 2) dtfrotate = dtf / 0.5; // for discs the formula is I=0.5*Mass*Radius^2
  else dtfrotate  = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
          vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
      if (vsq > vlimitsq) {
          ncount++;
          scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
      }

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
          osq = omega[i][0]*omega[i][0] + omega[i][1]*omega[i][1] + omega[i][2]*omega[i][2];
      if (osq > omegalimitsq) {
          ncount++;
          scale = sqrt(omegalimitsq/osq);
          omega[i][0] *= scale;
          omega[i][1] *= scale;
          omega[i][2] *= scale;
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereLimit::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixNVESphereLimit::compute_scalar()
{
  double one = ncount;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}
