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
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "comm.h"
#include "error.h"
#include "force.h"
#include "pair_hybrid.h"
#include "suffix.h"
#include <string.h>

#include "pair_gran_proxy.h"
#include "granular_pair_style.h"

using namespace LAMMPS_NS;
using namespace LIGGGHTS::PairStyles;

PairGranProxy::PairGranProxy(LAMMPS * lmp) : PairGran(lmp), impl(NULL)
{
}

PairGranProxy::~PairGranProxy()
{
  delete impl;
}

void PairGranProxy::settings(int nargs, char ** args)
{
  delete impl;

  const char * style = force->pair_style;

  /*NL*/ // figure out our style, unlike Fix this is not stored
  if(force->pair_match("hybrid", 0)) {
    PairHybrid * hybrid = static_cast<PairHybrid*>(force->pair);
    for(int i = 0; i < hybrid->nstyles; i++) {
      if(hybrid->styles[i] == this) {
        const char * pair_style = hybrid->keywords[i];
        int64_t variant = Factory::instance().selectVariant(pair_style, nargs, args);
        impl = Factory::instance().create(pair_style, variant, lmp, this);
        style = pair_style;
        break;
      }
    }
  } else {
    int64_t variant = Factory::instance().selectVariant(style, nargs, args);
    impl = Factory::instance().create(style, variant, lmp, this);
  }

  int length = strlen(style);

  if(length > 4 && (strcmp(&style[length-4], "/omp") == 0)) {
    suffix_flag |= Suffix::OMP;
    respa_enable = 0;
  }

  if(impl) {
    impl->settings(nargs, args);
  } else {
    error->one(FLERR, "unknown contact model");
  }
}

void PairGranProxy::init_granular()
{
  impl->init_granular();
}

void PairGranProxy::write_restart_settings(FILE * fp)
{
  impl->write_restart_settings(fp);
}

void PairGranProxy::read_restart_settings(FILE * fp)
{
  int me = comm->me;

  int64_t selected = -1;
  if(me == 0){
    // read model hashcode, but reset file pointer afterwards.
    // this way read_restart_settings can still read the hashcode (sanity check)
    fread(&selected, sizeof(int64_t), 1, fp);
    long offset = sizeof(int64_t);
    fseek(fp, -offset, SEEK_CUR);
  }
  MPI_Bcast(&selected,8,MPI_CHAR,0,world);

  impl = Factory::instance().create("gran", selected, lmp, this);

  if(impl) {
    impl->read_restart_settings(fp);
  } else {
    error->one(FLERR, "unknown contact model");
  }
}

void PairGranProxy::compute_force(int eflag, int vflag, int addflag)
{
  impl->compute_force(this, eflag, vflag, addflag);
}

void PairGranProxy::compute_single_pair_force(LCM::CollisionData &cdata, LCM::ForceData &i_forces, LCM::ForceData &j_forces)
{
  impl->compute_single_pair_force(cdata, i_forces, j_forces);
}

double PairGranProxy::stressStrainExponent()
{
  return impl->stressStrainExponent();
}

int64_t PairGranProxy::hashcode() {
  return impl->hashcode();
}
