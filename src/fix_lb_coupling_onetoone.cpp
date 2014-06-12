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
   Contributing authors:
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#include "fix_lb_coupling_onetoone.h"

#include "lammps.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "comm.h"

#include "fix_property_atom.h"

#include <iomanip>
#include <iostream>

namespace LAMMPS_NS {
  FixLbCouplingOnetoone::FixLbCouplingOnetoone(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp,narg,arg), fix_dragforce_(0), fix_hdtorque_(0)
  {
    std::cout << "new fix here" << std::endl;
  }

  FixLbCouplingOnetoone::~FixLbCouplingOnetoone()
  {

  }

  int FixLbCouplingOnetoone::setmask()
  {
    int mask = 0;
    mask |= FixConst::POST_FORCE;
    return mask;
  }

  void FixLbCouplingOnetoone::init()
  {
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");


  }

  void FixLbCouplingOnetoone::post_create()
  {
    // register dragforce
    if(!fix_dragforce_)
      {
        const char* fixarg[11];
        fixarg[0]="dragforce";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dragforce";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="yes";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_dragforce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }

    // register hydrodynamic torque
    if(!fix_hdtorque_)
      {
        const char* fixarg[11];
        fixarg[0]="hdtorque";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="hdtorque";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="yes";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_hdtorque_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }
  }


  void FixLbCouplingOnetoone::post_force(int)
  {
    double **f_ext = fix_dragforce_->array_atom;
    //double **t_ext = fix_hdtorque_->array_atom;
    double **f = atom->f;
    //double **t = atom->torque;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts "
    //             << std::setprecision(12) << f_ext[i][2] << std::endl;

    // std::cout << "before" << std::endl;
    // fix_dragforce_->do_reverse_comm();
    // fix_hdtorque_->do_reverse_comm();
    // std::cout << "after" << std::endl;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts_after "
    //             << std::setprecision(12) << f_ext[i][2] << " x " << atom->x[0][2] << std::endl;

    for(int i=0;i<atom->nlocal;i++){
      for(int j=0;j<3;j++){
        f[i][j] += f_ext[i][j];
        // t[i][j] += t_ext[i][j];
      }
    }

  }

  double **FixLbCouplingOnetoone::get_force_ptr()
  {
    return fix_dragforce_->array_atom;
  }
  double **FixLbCouplingOnetoone::get_torque_ptr()
  {
    return fix_hdtorque_->array_atom;
  }
  void FixLbCouplingOnetoone::comm_force_torque()
  {
    fix_dragforce_->do_reverse_comm();
    fix_hdtorque_->do_reverse_comm();
  }


}; /* LAMMPS_NS */
