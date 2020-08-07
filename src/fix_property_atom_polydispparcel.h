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

#ifdef FIX_CLASS

FixStyle(property/atom/polydispparcel,FixPropertyAtomPolydispParcel)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_POLYDISPPARCEL_H
#define LMP_FIX_PROPERTY_ATOM_POLYDISPPARCEL_H

#include "fix_property_atom.h"
#include "input.h"

namespace LAMMPS_NS {

class FixPropertyAtomPolydispParcel : public FixPropertyAtom {
 friend class Set;
 friend class FixPropertyAtomUpdateFix;
 friend class FixPropertyAtomRandom;
 public:
  FixPropertyAtomPolydispParcel(class LAMMPS *, int, char **,bool parse = true);
  ~FixPropertyAtomPolydispParcel();
  void init();

  Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,const char *caller,bool errflag);

  void pre_set_arrays();
  void set_arrays(int);
  void set_all(double value);
  
  void set_vector(int, double);

 protected:
  void parse_args(int narg, char **arg);
  
  int me;
  int ndefaultvalues;


}; //end class

}
#endif
#endif
