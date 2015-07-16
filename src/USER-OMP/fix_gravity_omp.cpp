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

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple U), Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "atom.h"
#include "update.h"
#include "modify.h"
#include "variable.h"
#include "fix_multisphere.h"
#include "fix_gravity_omp.h"
#include "thr_omp.h"
#include "domain.h"
#include "input.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{CHUTE,SPHERICAL,GRADIENT,VECTOR};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixGravityOMP::FixGravityOMP(LAMMPS *lmp, int narg, char **arg) :
  FixGravity(lmp, narg, arg) { }

/* ---------------------------------------------------------------------- */

int FixGravityOMP::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGravityOMP::post_force(int vflag)
{
  // update gravity due to variables

  if (varflag != CONSTANT) {
    modify->clearstep_compute();
    if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
    if (vstyle == EQUAL) magnitude = input->variable->compute_equal(vvar);
    if (pstyle == EQUAL) magnitude = input->variable->compute_equal(pvar);
    if (tstyle == EQUAL) magnitude = input->variable->compute_equal(tvar);
    if (xstyle == EQUAL) magnitude = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) magnitude = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) magnitude = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);

    set_acceleration();
  }

  eflag = 0;
  egrav = 0.0;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  double * const rmass = atom->rmass;
  double * const mass = atom->mass;
  int * const mask = atom->mask;
  int * const type = atom->type;
  const double xacc_thr = xacc;
  const double yacc_thr = yacc;
  const double zacc_thr = zacc;

  #pragma omp parallel default(none)
  {
    int ifrom;
    int ito;
    int tid;
    loop_setup_thr(ifrom, ito, tid, nall, nthreads);

    double grav = 0.0;

    if (rmass) {
      for (int i = ifrom; i < ito; i++) {
        if (mask[i] & groupbit) {
          const double massone = rmass[i];
          f[i][0] += massone*xacc_thr;
          f[i][1] += massone*yacc_thr;
          f[i][2] += massone*zacc_thr;
          grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
        }
      }
    } else {
      for (int i = ifrom; i < ito; i++) {
        if (mask[i] & groupbit) {
          const double massone = mass[type[i]];
          f[i][0] += massone*xacc_thr;
          f[i][1] += massone*yacc_thr;
          f[i][2] += massone*zacc_thr;
          grav -= massone * (xacc_thr*x[i][0] + yacc_thr*x[i][1] + zacc_thr*x[i][2]);
        }
      }
    }

    #pragma omp atomic
    egrav += grav;
  }
}
