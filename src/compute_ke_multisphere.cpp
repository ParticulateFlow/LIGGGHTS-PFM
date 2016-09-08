/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include "compute_ke_multisphere.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_multisphere.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeKEMultisphere::ComputeKEMultisphere(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  fix_ms_(0)
{
  if (narg != 3) error->compute_error(FLERR,this,"");

  scalar_flag = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

ComputeKEMultisphere::~ComputeKEMultisphere()
{
}

/* ---------------------------------------------------------------------- */

void ComputeKEMultisphere::init()
{
  fix_ms_ =  static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));
  if (!fix_ms_)
    error->compute_error(FLERR,this,"fix multisphere does not exist");
}

/* ---------------------------------------------------------------------- */

double ComputeKEMultisphere::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  scalar = fix_ms_->extract_ke();
  /*NL*/ fprintf(screen,"ke %f\n",scalar);
  return scalar;
}
