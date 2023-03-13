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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::FixCfdCouplingConvection(LAMMPS *lmp, int narg, char **arg) :  Fix(lmp, narg, arg)
{
    fix_coupling = NULL;
    fix_conductiveFlux = fix_convectiveFlux  = fix_heatFlux = NULL;
    T0 = -1.0;
    gran_field_conduction = false;
    limit_change = false;
    max_change = -1.0;

    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"T0") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'T0'");
            iarg++;
            T0 = atof(arg[iarg]);
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_conduction") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_conduction'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                gran_field_conduction = true;
            else if(strcmp(arg[iarg],"no") == 0)
                gran_field_conduction = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_conduction'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"max_change") == 0)
        {
            iarg++;
            max_change = atof(arg[iarg]);
            limit_change = true;
            iarg++;
            hasargs = true;
        }
    }

    if(T0 < 0.) error->all(FLERR,"Fix couple/cfd/convection: T0 must be >= 0. Specify a meaningful value.");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::~FixCfdCouplingConvection()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::pre_delete(bool unfixflag)
{
    if(fix_conductiveFlux) modify->delete_fix("conductiveHeatFlux");
    if(fix_convectiveFlux) modify->delete_fix("convectiveHeatFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvection::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_create()
{
  //  register convective flux
  if(!fix_convectiveFlux)
  {
        const char* fixarg[11];
        fixarg[0]="convectiveHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="convectiveHeatFlux";
        fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
        fixarg[5]="no";    //NP restart yes
        fixarg[6]="yes";    //NP communicate ghost no
        fixarg[7]="no";    //NP communicate rev yes
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  register conductive flux
  if(!fix_conductiveFlux && gran_field_conduction)
  {
        const char* fixarg[11];
        fixarg[0]="conductiveHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="conductiveHeatFlux";
        fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
        fixarg[5]="no";    //NP restart yes
        fixarg[6]="yes";    //NP communicate ghost no
        fixarg[7]="no";    //NP communicate rev yes
        fixarg[8]="0.";
        fix_conductiveFlux = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  add heat transfer model if not yet active
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  int nargs = 15;
  if (limit_change) nargs += 2;
  if(!fix_ste)
  {
        const char *newarg[nargs];
        newarg[0] = "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = "transportequation/scalar";
        newarg[3] = "equation_id";
        newarg[4] = "heattransfer";
        newarg[5] = "quantity";
        newarg[6] = "Temp";
        newarg[7] = "default_value";
        char arg8[30];
        sprintf(arg8,"%f",T0);
        newarg[8] = arg8;
        newarg[9] = "flux_quantity";
        newarg[10] = "heatFlux";
        newarg[11] = "source_quantity";
        newarg[12] = "heatSource";
        newarg[13] = "capacity_quantity";
        newarg[14] = "thermalCapacity";
        if (limit_change)
        {
            newarg[15] = "max_change";
            char arg16[30];
            sprintf(arg16,"%f",max_change);
            newarg[16] = arg16;
        }
        modify->add_fix(nargs,const_cast<char**>(newarg));
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling)
      error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    //values to send to OF
    fix_coupling->add_push_property("Temp","scalar-atom");

    //values to come from OF
    fix_coupling->add_pull_property("convectiveHeatFlux","scalar-atom");

    // heat transfer added heatFlux, get reference to it
    fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));

    fix_convectiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0,style));

    // granular field conduction?
    if (gran_field_conduction)
    {
      fix_coupling->add_pull_property("conductiveHeatFlux","scalar-atom");
      fix_conductiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("conductiveHeatFlux","property/atom","scalar",0,0,style));
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_force(int)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // communicate convective flux to ghosts, there might be new data
  
  if(0 == neighbor->ago)
  {
      fix_convectiveFlux->do_forward_comm();
      if (gran_field_conduction)
      {
          fix_conductiveFlux->do_forward_comm();
      }
  }

  double *heatFlux = fix_heatFlux->vector_atom;
  double *conductiveFlux = NULL;
  if (gran_field_conduction) conductiveFlux = fix_conductiveFlux->vector_atom;
  double *convectiveFlux = fix_convectiveFlux->vector_atom;


  //NP add convective flux to per-particle heat flux
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
    {
        heatFlux[i] += convectiveFlux[i];
        if (gran_field_conduction)
        {
            heatFlux[i] += conductiveFlux[i];
        }
    }
}
