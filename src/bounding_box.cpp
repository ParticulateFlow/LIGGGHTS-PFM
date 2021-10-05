/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "bounding_box.h"
#include <algorithm>

namespace LAMMPS_NS
{

  BoundingBox::BoundingBox()
  : xLo(0.), xHi(0.), yLo(0.), yHi(0.), zLo(0.), zHi(0.), initGiven(false), dirty(true)
  {}
  BoundingBox::BoundingBox(double bounds[6])
  : xLo(bounds[0]), xHi(bounds[1]), yLo(bounds[2]), yHi(bounds[3]), zLo(bounds[4]), zHi(bounds[5]), initGiven(true), dirty(true)
  {}
  BoundingBox::BoundingBox(double xLo_, double xHi_, double yLo_, double yHi_, double zLo_, double zHi_)
  : xLo(xLo_), xHi(xHi_), yLo(yLo_), yHi(yHi_), zLo(zLo_), zHi(zHi_), initGiven(true), dirty(true)
  {}

  BoundingBox::~BoundingBox()
  {}

  void BoundingBox::reset()
  {
    xLo = 0.; xHi = 0.;
    yLo = 0.; yHi = 0.;
    zLo = 0.; zHi = 0.;
    initGiven = false;
    dirty = true;
  }

  void BoundingBox::extendByDelta(double delta)
  {
    xLo = xLo-delta;
    yLo = yLo-delta;
    zLo = zLo-delta;
    xHi = xHi+delta;
    yHi = yHi+delta;
    zHi = zHi+delta;
  }

  void BoundingBox::extrude(double length, const double * vec)
  {
    xLo = std::min(xLo, (xLo + length * vec[0]));
    yLo = std::min(yLo, (yLo + length * vec[1]));
    zLo = std::min(zLo, (zLo + length * vec[2]));
    xHi = std::max(xHi, (xHi + length * vec[0]));
    yHi = std::max(yHi, (yHi + length * vec[1]));
    zHi = std::max(zHi, (zHi + length * vec[2]));
  }

} /* namespace LAMMPS_NS */
