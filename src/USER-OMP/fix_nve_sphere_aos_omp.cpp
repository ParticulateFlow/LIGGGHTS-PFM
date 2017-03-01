/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2013-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_sphere_aos_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const double INERTIA=0.4;          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

struct vector3_t {
  double x;
  double y;
  double z;
};

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixNVESphereAOSOMP::initial_integrate(int vflag)
{
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  if(nlocal == 0) return;

  vector3_t * x = (vector3_t*)&atom->x[0][0];
  vector3_t * v = (vector3_t*)&atom->v[0][0];
  vector3_t * f = (vector3_t*)&atom->f[0][0];
  vector3_t * omega  = (vector3_t*)&atom->omega[0][0];
  vector3_t * torque = (vector3_t*)&atom->torque[0][0];
  const double * const radius = atom->radius;
  const double * const rmass = atom->rmass;
  const int * const mask = atom->mask;

  // set timestep here since dt may have changed or come via rRESPA
  const double dtfrotate = dtf / INERTIA;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia
  #pragma omp parallel for
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i].x += dtfm * f[i].x;
      v[i].y += dtfm * f[i].y;
      v[i].z += dtfm * f[i].z;
      x[i].x += dtv * v[i].x;
      x[i].y += dtv * v[i].y;
      x[i].z += dtv * v[i].z;

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i].x += dtirotate * torque[i].x;
      omega[i].y += dtirotate * torque[i].y;
      omega[i].z += dtirotate * torque[i].z;
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double * const * const mu = atom->mu;
    #pragma omp parallel for
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        if (mu[i][3] > 0.0) {
          const double g0 = mu[i][0] + dtv * (omega[i].y*mu[i][2]-omega[i].z*mu[i][1]);
          const double g1 = mu[i][1] + dtv * (omega[i].z*mu[i][0]-omega[i].x*mu[i][2]);
          const double g2 = mu[i][2] + dtv * (omega[i].x*mu[i][1]-omega[i].y*mu[i][0]);
          const double msq = g0*g0 + g1*g1 + g2*g2;
          const double scale = mu[i][3]/sqrt(msq);
          mu[i][0] = g0*scale;
          mu[i][1] = g1*scale;
          mu[i][2] = g2*scale;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereAOSOMP::final_integrate()
{
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  if(nlocal == 0) return;

  vector3_t * v = (vector3_t*)&atom->v[0][0];
  vector3_t * f = (vector3_t*)&atom->f[0][0];
  vector3_t * omega  = (vector3_t*)&atom->omega[0][0];
  vector3_t * torque = (vector3_t*)&atom->torque[0][0];
  const double * const rmass = atom->rmass;
  const double * const radius = atom->radius;
  const int * const mask = atom->mask;

  // set timestep here since dt may have changed or come via rRESPA

  const double dtfrotate = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

  #pragma omp parallel for
  for (int i = 0; i < nlocal; ++i) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i].x += dtfm * f[i].x;
      v[i].y += dtfm * f[i].y;
      v[i].z += dtfm * f[i].z;

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i].x += dtirotate * torque[i].x;
      omega[i].y += dtirotate * torque[i].y;
      omega[i].z += dtirotate * torque[i].z;
    }
  }
}
