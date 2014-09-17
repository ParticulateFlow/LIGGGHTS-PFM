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
#include "fix_nve_limit_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

void FixNVELimitOMP::initial_integrate(int)
{
  double * const * const x = atom->x;
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  const double * const mass = atom->mass;
  const double * const rmass = atom->rmass;
  const double * const radius = atom->radius;
  const int * const type = atom->type;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int ncount_local = 0;

  if (rmass) {
    if (relflag == 1) //NP modified C.K.
    {
      #pragma omp parallel for reduction(+:ncount_local)
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          const double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          const double rsq = radius[i]*radius[i];
          if (vsq > vlimitsq*rsq) {
            ++ncount_local;
            const double scale = sqrt(vlimitsq*rsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }

          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
    else
    {
      #pragma omp parallel for reduction(+:ncount_local)
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          const double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          if (vsq > vlimitsq) {
            ++ncount_local;
            const double scale = sqrt(vlimitsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }

          x[i][0] += dtv * v[i][0];
          x[i][1] += dtv * v[i][1];
          x[i][2] += dtv * v[i][2];
        }
      }
    }
  } else {
    #pragma omp parallel for reduction(+:ncount_local)
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ++ncount_local;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
    }
  }
  ncount += ncount_local;
}

/* ---------------------------------------------------------------------- */

void FixNVELimitOMP::final_integrate()
{
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  int ncount_local = 0;

  if (rmass) {
    if (relflag == 1) //NP modified C.K.
    {
      #pragma omp parallel for reduction(+:ncount_local)
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          const double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          const double rsq = radius[i]*radius[i];
          if (vsq > vlimitsq*rsq) {
            ++ncount_local;
            const double scale = sqrt(vlimitsq*rsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }
        }
      }
    }
    else
    {
      #pragma omp parallel for reduction(+:ncount_local)
      for (int i = 0; i < nlocal; ++i) {
        if (mask[i] & groupbit) {
          const double dtfm = dtf / rmass[i];
          v[i][0] += dtfm * f[i][0];
          v[i][1] += dtfm * f[i][1];
          v[i][2] += dtfm * f[i][2];

          const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
          if (vsq > vlimitsq) {
            ++ncount_local;
            const double scale = sqrt(vlimitsq/vsq);
            v[i][0] *= scale;
            v[i][1] *= scale;
            v[i][2] *= scale;
          }
        }
      }
    }

  } else {
    #pragma omp parallel for reduction(+:ncount_local)
    for (int i = 0; i < nlocal; ++i) {
      if (mask[i] & groupbit) {
        const double dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ++ncount_local;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }
      }
    }
  }
  ncount += ncount_local;
}
