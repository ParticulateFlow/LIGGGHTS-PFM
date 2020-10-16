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
#include "fix_nve_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"
#include "domain.h" //NP modified GM
#include "fix_property_atom.h"
#include "fix_cfd_coupling_force_implicit.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVESphere::FixNVESphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg),
  implicitIntegration_(false),
  fix_Ksl_(0),
  fix_cfd_coupling_force_implicit_(0)
{
  if (narg < 3) error->all(FLERR,"Illegal fix nve/sphere command");

  time_integrate = 1;

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix nve/sphere command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else error->all(FLERR,"Illegal fix nve/sphere command");
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"implicit_integration") == 0)
    {
        if(narg < iarg+2)
            error->fix_error(FLERR,this,"not enough arguments for 'transfer_reactant'");
        iarg++;
        if(strcmp(arg[iarg],"yes") == 0)
            implicitIntegration_ = true;
        else if(strcmp(arg[iarg],"no") == 0)
            implicitIntegration_ = false;
        else
            error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'implicit_integration'");
        iarg++;
    }
    else error->all(FLERR,"Illegal fix nve/sphere command");
  }

  // error checks

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nve/sphere requires atom style sphere");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all(FLERR,"Fix nve/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::post_create()
{
    if (implicitIntegration_)
    {
        fix_Ksl_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("Ksl","property/atom","scalar",0,0,style,false));
        if (!fix_Ksl_)
        {
            error->fix_error(FLERR,this,"Fix NVE/sphere could not find fix 'Ksl' for the drag coefficient");
        }

        fix_cfd_coupling_force_implicit_ = static_cast<FixCfdCouplingForceImplicit*>(modify->find_fix_style_strict("couple/cfd/force/implicit",0));
        if (!fix_cfd_coupling_force_implicit_)
        {
            error->fix_error(FLERR,this,"Could not find fix ID 'couple/cfd/force/implicit'");
        }

        if (!fix_cfd_coupling_force_implicit_->implicitIntegration())
        {
            error->fix_error(FLERR,this,"Fix 'couple/cfd/force/implicit' is not in mode for implicit integration.");
        }
    }
}

void FixNVESphere::init()
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

void FixNVESphere::initial_integrate(int vflag)
{
  double dtfm,dtirotate,msq,scale;
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
  double invImpDenom;
  double *Ksl;
  if (implicitIntegration_) Ksl = fix_Ksl_->vector_atom;

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
      if (implicitIntegration_)
      {
          invImpDenom = 1.0/(1.0 + dtfm * Ksl[i]);
          v[i][0] *= invImpDenom;
          v[i][1] *= invImpDenom;
          v[i][2] *= invImpDenom;
      }
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
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

void FixNVESphere::final_integrate()
{
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double invImpDenom;
  double *Ksl;
  if (implicitIntegration_) Ksl = fix_Ksl_->vector_atom;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate; //NP modified GM
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
      if (implicitIntegration_)
      {
          invImpDenom = 1.0/(1.0 + dtfm * Ksl[i]);
          v[i][0] *= invImpDenom;
          v[i][1] *= invImpDenom;
          v[i][2] *= invImpDenom;
      }

      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}
