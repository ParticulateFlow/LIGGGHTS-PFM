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
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_nve_sphere_partitioned_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"

#include <omp.h>

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixNVESpherePartitionedOMP::initial_integrate(int vflag)
{
  #pragma omp parallel default(none)
  {
    double * const * const x = atom->x;
    double * const * const v = atom->v;
    const double * const * const f = atom->f;
    double * const * const omega = atom->omega;
    const double * const * const torque = atom->torque;
    const double * const radius = atom->radius;
    const double * const rmass = atom->rmass;
    const int * const mask = atom->mask;

    const int tid = omp_get_thread_num();
    const int ifrom = atom->thread_offsets[tid];
    const int ito   = atom->thread_offsets[tid+1];

    // set timestep here since dt may have changed or come via rRESPA
    const double dtfrotate = dtf / INERTIA;

    // update v,x,omega for all particles
    // d_omega/dt = torque / inertia
    for (int i = ifrom; i < ito; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];

        const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
        omega[i][0] += dtirotate * torque[i][0];
        omega[i][1] += dtirotate * torque[i][1];
        omega[i][2] += dtirotate * torque[i][2];
      }
    }

    // update mu for dipoles
    // d_mu/dt = omega cross mu
    // renormalize mu to dipole length

    if (extra == DIPOLE) {
      double * const * const mu = atom->mu;
      for (int i = ifrom; i < ito; ++i) {
        if (mask[i] & groupbit) {
          if (mu[i][3] > 0.0) {
            const double g0 = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
            const double g1 = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
            const double g2 = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
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
}

/* ---------------------------------------------------------------------- */

void FixNVESpherePartitionedOMP::final_integrate()
{
  #pragma omp parallel default(none)
  {
    double * const * const v = atom->v;
    const double * const * const f = atom->f;
    double * const * const omega = atom->omega;
    const double * const * const torque = atom->torque;
    const double * const rmass = atom->rmass;
    const double * const radius = atom->radius;
    const int * const mask = atom->mask;
    const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal; // TODO handle atom->nfirst case

    const int tid = omp_get_thread_num();
    const int ifrom = atom->thread_offsets[tid];
    const int ito   = atom->thread_offsets[tid+1];

    // set timestep here since dt may have changed or come via rRESPA

    const double dtfrotate = dtf / INERTIA;

    // update v,omega for all particles
    // d_omega/dt = torque / inertia

    for (int i = ifrom; i < ito; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
        omega[i][0] += dtirotate * torque[i][0];
        omega[i][1] += dtirotate * torque[i][1];
        omega[i][2] += dtirotate * torque[i][2];
      }
    }
  }
}
