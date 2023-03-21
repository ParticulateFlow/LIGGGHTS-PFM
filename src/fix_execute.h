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

#ifdef FIX_CLASS

FixStyle(execute,FixExecute)

#else

#ifndef LMP_FIX_EXECUTE_H
#define LMP_FIX_EXECUTE_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixExecute : public Fix {
 public:
  FixExecute(class LAMMPS *, int, char **);
  ~FixExecute();
  int setmask();
  void initial_integrate(int);
  void end_of_step();

 private:
  int me;
  char *string;

  int execution_point; // 0 initial_integrate, 1 end_of_step

  // conditional execution
  bool conditional;
  // execute whole input file
  bool file;
  // variable used to identify if insertion should take place
  char *var;
  double var_valid;
  double var_threshold;

  // single execution
  bool once;
  int execution_step;

  void execution_command();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line

*/
