/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

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
   Thomas Lichtenegger (JKU Linz)
   M.Efe Kinaci (JKU Linz)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(chem/shrink/Arrhenius,FixChemShrinkArrhenius)

#else

#ifndef LMP_FIX_CHEM_SHRINK_ARRHENIUS_H
#define LMP_FIX_CHEM_SHRINK_ARRHENIUS_H

#include "fix_chem_shrink.h"

namespace LAMMPS_NS {

class FixChemShrinkArrhenius : public FixChemShrink  {

public:
  FixChemShrinkArrhenius(class LAMMPS *, int, char **);
  ~FixChemShrinkArrhenius();
  
  void updatePtrs();

 

 protected:

  // values from user                                       
  double T0;			// activation temperature
  
  double reactionRatConst(int);


};

}

#endif
#endif


