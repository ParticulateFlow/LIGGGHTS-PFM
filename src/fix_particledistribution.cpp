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

#include "fix_particledistribution.h"
#include "modify.h"
#include "random_park.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixParticledistribution::FixParticledistribution(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  pti(NULL),
  pti_list(NULL),
  n_pti(0),
  n_pti_max(0),
  ninserted(0),
  ninsert(0),
  random(NULL),
  seed(1),
  volexpect(0.),
  massexpect(0.),
  maxtype(0),
  mintype(0),
  maxnspheres(0),
  minrad(0),
  maxrad(0),
  maxrbound(0)
{
  restart_global = 1;
}

/* ---------------------------------------------------------------------- */

FixParticledistribution::~FixParticledistribution()
{
  delete []pti_list;
  delete random;
}

/* ----------------------------------------------------------------------*/

int FixParticledistribution::setmask()
{
  int mask = 0;
  return mask;
}


/* ----------------------------------------------------------------------
   preparations before insertion
------------------------------------------------------------------------- */

void FixParticledistribution::pre_insert(int n, FixPropertyAtom *fp, double val, int idx, int ival, int iidx)
{
  // allow fixes to e.g. update some pointers before set_arrays is called
  // set_arrays called in ParticleToInsert::insert()

  int nfix = modify->nfix;
  Fix **fix = modify->fix;

  for (int j = 0; j < nfix; j++)
    if (fix[j]->create_attribute)
      fix[j]->pre_set_arrays();
}

/* ----------------------------------------------------------------------
   set particle properties - only pti needs to know which properties to set
   loop to n, not n_pti, since not all particles may have been inserted
------------------------------------------------------------------------- */

int FixParticledistribution::insert(int /*n*/)
{
  return 0;
}

/* ----------------------------------------------------------------------
   wrap up insertion
------------------------------------------------------------------------- */

void FixParticledistribution::finalize_insertion()
{
}


/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixParticledistribution::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixParticledistribution::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]) + comm->me;

  random->reset(seed);
}
