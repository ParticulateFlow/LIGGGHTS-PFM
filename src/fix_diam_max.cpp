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

#include "fix_diam_max.h"
#include "modify.h"
#include "math_extra_liggghts.h"
#include "fix_particledistribution_discrete.h"


using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixDiamMax::FixDiamMax(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  maxrbound_(0.)
{
  scalar_flag = 1;

  FixParticledistributionDiscrete *fpdd;
  int nfix = modify->n_fixes_style("particledistribution/discrete");
  for(int i = 0; i < nfix; i++)
  {
      fpdd = static_cast<FixParticledistributionDiscrete*>(modify->find_fix_style("particledistribution/discrete",0));
      maxrbound_ = max(maxrbound_,fpdd->max_r_bound());
  }
}

/* ---------------------------------------------------------------------- */

FixDiamMax::~FixDiamMax()
{

}

/* ----------------------------------------------------------------------*/

int FixDiamMax::setmask()
{
    int mask = 0;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixDiamMax::init()
{
  maxrbound_ = 0.;

  FixParticledistributionDiscrete *fpdd;
  int nfix = modify->n_fixes_style("particledistribution/discrete");
  for(int i = 0; i < nfix; i++)
  {
      fpdd = static_cast<FixParticledistributionDiscrete*>(modify->find_fix_style("particledistribution/discrete",0));
      maxrbound_ = max(maxrbound_,fpdd->max_r_bound());
  }
}

/* ---------------------------------------------------------------------- */

double FixDiamMax::compute_scalar()
{
    return 2.*maxrbound_;
}
