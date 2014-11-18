/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include "particleToInsert_fragments.h"
#include "math.h"
#include "math_const.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix_property_atom.h"
#include "vector_liggghts.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ParticleToInsertFragments::ParticleToInsertFragments(LAMMPS* lmp,int ns) : ParticleToInsert(lmp,ns)
{
  collision_factor = 1.0;
}

/* ---------------------------------------------------------------------- */

ParticleToInsertFragments::~ParticleToInsertFragments()
{
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertFragments::insert()
{
  // perform the actual insertion
  // add particles, set coordinate and radius
  // set group mask to "all" plus fix groups

  int inserted = 0;
  int nfix = modify->nfix;
  Fix **fix = modify->fix;
  FixPropertyAtom *fix_pa = NULL;
  if (!fix_property_atom_id.empty()) {
    fix_pa = static_cast<FixPropertyAtom*>(modify->find_fix_id(fix_property_atom_id.c_str()));
  }

  for(int i = 0; i < nspheres; ++i) {
    ++inserted;
    atom->avec->create_atom(atom_type,x_ins[i]);
    int m = atom->nlocal - 1;
    atom->mask[m] = 1 | groupbit;
    vectorCopy3D(v_ins,atom->v[m]);
    vectorCopy3D(omega_ins,atom->omega[m]);
    atom->radius[m] = radius_ins[i];
    atom->density[m] = density_ins;
    atom->rmass[m] = (MY_4PI3) * radius_ins[i]*radius_ins[i]*radius_ins[i] * density_ins;

    for (int j = 0; j < nfix; ++j) {
      if (fix[j]->create_attribute) fix[j]->set_arrays(m);
    }

    if (fix_pa) {
      fix_pa->vector_atom[m] = collision_factor;
    }
  }

  return inserted;
}
