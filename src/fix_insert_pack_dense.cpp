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
   contributing authors
   Philippe Seil (JKU)
   ---------------------------------------------------------------------- */


#include "fix_insert_pack_dense.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixInsertPackDense::FixInsertPackDense(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  printf("hello i am new fix\n");
}

/* ---------------------------------------------------------------------- */

FixInsertPackDense::~FixInsertPackDense()
{

}

/* ---------------------------------------------------------------------- */

