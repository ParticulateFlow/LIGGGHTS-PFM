/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_nve_sph_limit.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVESphLimit::FixNVESphLimit(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  dtv(0.),
  dtf(0.),
  step_respa(NULL),
  mass_require(false),
  ncount(0),
  xlimit(0.),
  vlimitsq(0.)
{
  if (narg != 5) error->all(FLERR,"Illegal fix nve/sph/limit command");

  time_integrate = 1;
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix sph command requires atom_style with both energy and density");

  if (strcmp(arg[3],"absolute") == 0)
    xlimit = force->numeric(FLERR,arg[4]);
  else
    error->fix_error(FLERR,this,"expecting keyword 'absolute'");
}

/* ---------------------------------------------------------------------- */

int FixNVESphLimit::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESphLimit::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv);
  ncount = 0;

  if (strstr(update->integrate_style,"respa"))
    step_respa = ((Respa *) update->integrate)->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVESphLimit::initial_integrate(int vflag)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];

        // same like src/USER-SPH/fix_meso !!
        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
        rho[i] += 2.0 * dtf * drho[i]; // ... and density by full step

        // euler integration for velocity (estimation for force calculatoin)
        vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
        vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
        vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];
        // half step
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        // check half-step and estimated velocity (only if half step exceeds the limit)
        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ncount++;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;

          const double vsqEst = vest[i][0]*vest[i][0] + vest[i][1]*vest[i][1] + vest[i][2]*vest[i][2];
          if (vsqEst > vlimitsq) {
            const double scaleEst = sqrt(vlimitsq/vsqEst);
            vest[i][0] *= scaleEst;
            vest[i][1] *= scaleEst;
            vest[i][2] *= scaleEst;
          }
        }

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
        rho[i] += 2.0 * dtf * drho[i]; // ... and density

        // euler integration for velocity
        vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
        vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
        vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        // check half-step and estimated velocity (only if half step exceeds the limit)
        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ncount++;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;

          const double vsqEst = vest[i][0]*vest[i][0] + vest[i][1]*vest[i][1] + vest[i][2]*vest[i][2];
          if (vsqEst > vlimitsq) {
            const double scaleEst = sqrt(vlimitsq/vsqEst);
            vest[i][0] *= scaleEst;
            vest[i][1] *= scaleEst;
            vest[i][2] *= scaleEst;
          }
        }

        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphLimit::final_integrate()
{
  double dtfm;

  // update v of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  //double *rho = atom->rho; //NP modified R.B.
  //const double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ncount++;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];

        const double vsq = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
        if (vsq > vlimitsq) {
          ncount++;
          const double scale = sqrt(vlimitsq/vsq);
          v[i][0] *= scale;
          v[i][1] *= scale;
          v[i][2] *= scale;
        }

        e[i] += dtf * de[i]; // half-step update of particle internal energy
        //rho[i] += dtf * drho[i]; // ... and density
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphLimit::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESphLimit::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVESphLimit::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  vlimitsq = (xlimit/dtv) * (xlimit/dtv);
}


/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixNVESphLimit::compute_scalar()
{
  double one = ncount;
  double all;
  MPI_Allreduce(&one,&all,1,MPI_DOUBLE,MPI_SUM,world);
  return all;
}

