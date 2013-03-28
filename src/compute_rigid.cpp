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
#include "string.h"
#include "compute_rigid.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "comm.h"
#include "container.h"
#include "fix_multisphere.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeRigid::ComputeRigid(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  multisphere_(0)
{
  if (narg != 5) error->compute_error(FLERR,this,"Illegal compute rigid command, expecting 5 arguments");

  local_flag = 1;

  //NP get reference to fix rigid

  if(modify->n_fixes_style("multisphere") != 1)
    error->compute_error(FLERR,this,"defining exactly one fix multisphere is required");
  multisphere_ = & (static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0))->data());

  // parse args

  int iarg = 3;
  if(strcmp(arg[iarg++],"property"))
    error->compute_error(FLERR,this,"expecting keyword 'property'");
  property_ = multisphere_->prop().getElementPropertyBase(arg[iarg++]);

  // check arg validity

  if(!property_)
    error->compute_error(FLERR,this,"illegal property name used");
  if(!property_->isDoubleData())
    error->compute_error(FLERR,this,"can only operate on data of type double");
}

/* ---------------------------------------------------------------------- */

ComputeRigid::~ComputeRigid()
{}

/* ---------------------------------------------------------------------- */

void ComputeRigid::update_pointers()
{
    size_local_rows = multisphere_->n_body();

    if(property_->lenVec() > 1)
    {
        size_local_cols = property_->lenVec();
        array_local = (double**) property_->begin_slow_dirty();
    }
    else
    {
        size_local_cols  = 0;
        vector_local = (double*) property_->begin_slow_dirty();
    }
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::init()
{
  //NP make a first call to get size_local_rows and size_local_cols
  update_pointers();
}

/* ---------------------------------------------------------------------- */

void ComputeRigid::compute_local()
{
  double **dmpptr;
  int nprops;

  invoked_local = update->ntimestep;

  update_pointers();
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeRigid::memory_usage()
{
  return 0;
}
