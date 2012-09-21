/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_balance.h"
#include "balance.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "comm.h"
#include "irregular.h"
#include "force.h"
#include "kspace.h"
#include "error.h"
#include "modify.h" //NP modified C.K.

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix balance command");

  box_change = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  // parse arguments

  int dimension = domain->dimension;

  nevery = atoi(arg[3]);
  if (strlen(arg[4]) > 3) error->all(FLERR,"Illegal fix balance command");
  strcpy(bstr,arg[4]);
  nitermax = atoi(arg[5]);
  thresh = atof(arg[6]);

  if (nevery < 0 || nitermax <= 0 || thresh < 1.0)
    error->all(FLERR,"Illegal fix balance command");

  for (int i = 0; i < strlen(bstr); i++) {
    if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
      error->all(FLERR,"Fix balance string is invalid");
    if (bstr[i] == 'z' && dimension == 2)
      error->all(FLERR,"Fix balance string is invalid for 2d simulation");
    for (int j = i+1; j < strlen(bstr); j++)
      if (bstr[i] == bstr[j])
        error->all(FLERR,"Fix balance string is invalid");
  }

  // optional args

  int outarg = 0;
  fp = NULL;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix balance command");
      outarg = iarg+1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix balance command");
  }

  // create instance of Balance class and initialize it with params
  // create instance of Irregular class

  balance = new Balance(lmp);
  balance->dynamic_setup(bstr,nitermax,thresh);
  irregular = new Irregular(lmp);

  //NP modified C.K.
  //NP disallow balancing for more than one dimension for meshes
  //NP b/c would need irregular() in case particle is at edge
  //NP P_b = before , P_a = after balance
  //NP     P_b|
  //NP  -----------
  //NP        |P_a
  //NP if(strlen(bstr) > 1 && balance->disallow_irregular())
  //NP  error->fix_error(FLERR,this,"Cannot balance in more than one dimention when meshes are present");

  // output file

  if (outarg && comm->me == 0) {
    fp = fopen(arg[outarg],"w");
    if (fp == NULL) error->one(FLERR,"Cannot open fix balance output file");
  }

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;

  // compute initial outputs

  imbfinal = imbprev = balance->imbalance_nlocal(maxperproc);
  itercount = 0;
  pending = 0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  if (fp) fclose(fp);
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  // don't allow PPPM for now

  if (force->kspace && strstr(force->kspace_style,"pppm"))
    error->all(FLERR,"Cannot yet use fix balance with PPPM");

  //NP modified C.K.
  //NP all fixes that insert particles must come before this fix
  //NP this is b/c balancing is done pre_exchange()
  //NP and insertion as well
  int my_i = modify->find_fix(id);

  int nfix = modify->nfix;
  Fix **fix = modify->fix;
  int ntypes = atom->ntypes;
  for(int ifix = 0; ifix < nfix; ifix++)
  {
      for(int itype = 0; itype < ntypes+1; itype++)
      {
          if(my_i < ifix && fix[ifix]->max_rad(itype) > 0.)
          {
              char errstr[200];
              sprintf(errstr,"Fix %s has to come before fix %s",fix[ifix]->style,style);
              error->fix_error(FLERR,this,errstr);
          }
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup(int vflag)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshhold exceeded

  imbnow = balance->imbalance_nlocal(maxperproc);
  if (imbnow > thresh) rebalance();

  // next_reneighbor = next time to force reneighboring

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // return if imbalance < threshhold

  imbnow = balance->imbalance_nlocal(maxperproc);
  if (imbnow <= thresh) {
    if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    return;
  }

  rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occured
------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_nlocal(maxperproc);
  pending = 0;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::rebalance()
{
  imbprev = imbnow;
  itercount = balance->dynamic();

  // output of final result

  if (fp) balance->dumpout(update->ntimestep,fp);

  // reset comm->uniform flag

  comm->uniform = 0;

  // reset proc sub-domains

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();

  // if splits moved further than neighboring processor
  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far,
  // else allow caller's comm->exchange() to do it

  //NP modified C.K.

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (irregular->migrate_check())
  {
      if(balance->disallow_irregular() && 0 == comm->me)
        error->warning(FLERR,"Processor boundaries moved too far - mesh elements might be lost");
      irregular->migrate_atoms();
  }
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // NOTE: still to be implemented
  // check that new sub-domains are valid with KSpace constraints
  // if (kspace_flag) force->kspace->check();

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in caller's comm->exchange()

  pending = 1;
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return (double) maxperproc;
  if (i == 1) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = irregular->memory_usage();
  return bytes;
}
