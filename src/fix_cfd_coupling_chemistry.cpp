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
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_chemistry.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingChemistry::FixCfdCouplingChemistry(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    fix_coupling_(0),
    use_Re_(false)
{
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"transfer_Re) == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_Re'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_Re_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_Re_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_Re'");
            iarg++;
            hasargs = true;
        } else if (strcmp(this->style,"couple/cfd/chemistry") == 0) {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    // flags for vector output
    vector_flag = 1;
    size_vector = 3;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingChemistry::~FixCfdCouplingChemistry()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_create()
{
   
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::pre_delete(bool unfixflag)
{
    
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingChemistry::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/chemistry needs a fix of type couple/cfd");

    //  values to be transfered to OF
    fix_coupling_->add_push_property("reactionheat","scalar-atom");




    // values to come from OF
    if(use_Re_) fix_coupling_->add_pull_property("Re","scalar-atom");


}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingChemistry::post_force(int)
{
  
// for all species names i and reaction heat
// if dc_->pushednow(i)
// clear masschange(i),reactionheat
// pushednow(i)=false
  
  
  
}

