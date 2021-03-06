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
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_FIX_WALL_GRAN_BASE_H
#define LMP_FIX_WALL_GRAN_BASE_H

#include "fix_wall_gran.h"
#include "fix_contact_property_atom_wall.h"
#include "contact_interface.h"
#include "compute_pair_gran_local.h"
#include "settings.h"
#include <string.h>
#include "force.h"
#include <stdlib.h>
#include "contact_models.h"
#include "granular_wall.h"

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_const.h"
#endif

/*NL*/ #define DEBUGMODE_WALL_HOOKE false
/*NL*/ #define DEBUG_OUTP_WALL_HOOKE screen

namespace LIGGGHTS {
using namespace ContactModels;
namespace Walls {

template<typename ContactModel>
class Granular : private Pointers, public IGranularWall {
  ContactModel cmodel;
  FixWallGran * parent;

public:
  Granular(LAMMPS * lmp, FixWallGran * parent) :
    Pointers(lmp),
    cmodel(lmp, parent),
    parent(parent)
  {
  }

  virtual ~Granular() {
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      if(screen) {
        fprintf(screen, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
        force->registry.print_all(screen);
        fprintf(screen, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
      }

      if(logfile) {
        fprintf(logfile, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
        force->registry.print_all(logfile);
        fprintf(logfile, "==== WALL %s GLOBAL PROPERTIES ====\n", parent->id);
      }
    }
#endif
  }


  int64_t hashcode()
  { return cmodel.hashcode(); }

  virtual void settings(int nargs, char ** args) {
    Settings settings(lmp);
    cmodel.registerSettings(settings);
    bool success = settings.parseArguments(nargs, args);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      if(screen) {
        fprintf(screen, "==== WALL %s SETTINGS ====\n", parent->id);
        settings.print_all(screen);
        fprintf(screen, "==== WALL %s SETTINGS ====\n", parent->id);
      }

      if(logfile) {
        fprintf(logfile, "==== WALL %s SETTINGS ====\n", parent->id);
        settings.print_all(logfile);
        fprintf(logfile, "==== WALL %s SETTINGS ====\n", parent->id);
      }
    }
#endif

    if(!success) {
      error->fix_error(FLERR, parent, settings.error_message.c_str());
    }
  }

  inline void force_update(double * const f, double * const torque,
      const ForceData & forces) {
    for (int coord = 0; coord < 3; coord++) {
      f[coord] += forces.delta_F[coord];
      torque[coord] += forces.delta_torque[coord];
    }
  }

  virtual void compute_force(FixWallGran * wg, CollisionData & cdata, bool intersectflag, double *vwall, ForceData & i_forces, ForceData & j_forces)
  {
    const int ip = cdata.i;

    double *f = atom->f[ip];
    double *torque = atom->torque[ip];
    double *v = atom->v[ip];
    double *omega = atom->omega[ip];
    double mass = atom->rmass[ip];
    int *type = atom->type;

    if(wg->fix_rigid() && wg->body(ip) >= 0)
      mass = wg->masstotal(wg->body(ip));

    /*NL*///if(DEBUGMODE_WALL_HOOKE) fprintf(lmp->DEBUG_OUTP_WALL_HOOKE, "WALL_GRAN_compute_force 0\n");
    /*NL*///if (screen) printVec3D(screen,"shear history",c_history);

    const double r = cdata.r;
    const double rinv = 1.0/r;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
    double enx, eny, enz;
    if (atom->superquadric_flag) {
      const double delta_inv = 1.0 / vectorMag3D(cdata.delta);
      enx = cdata.delta[0] * delta_inv;
      eny = cdata.delta[1] * delta_inv;
      enz = cdata.delta[2] * delta_inv;
      cdata.radi = cbrt(0.75 * atom->volume[ip] / M_PI);
    } else { // sphere case
      enx = cdata.delta[0] * rinv;
      eny = cdata.delta[1] * rinv;
      enz = cdata.delta[2] * rinv;
    }
#else // sphere case
    const double enx = cdata.delta[0] * rinv;
    const double eny = cdata.delta[1] * rinv;
    const double enz = cdata.delta[2] * rinv;
#endif
    // copy collision data to struct (compiler can figure out a better way to
    // interleave these stores with the double calculations above.
    cdata.v_i = v;
    cdata.v_j = vwall;
    cdata.omega_i = omega;
    cdata.en[0] = enx;
    cdata.en[1] = eny;
    cdata.en[2] = enz;
    cdata.i = ip;
    cdata.touch = NULL;
    cdata.itype = type[ip];

    cdata.r = r;
    cdata.rinv = rinv;
    cdata.radsum = cdata.radi;
    cdata.mi = mass;

    if (intersectflag)
    {
      cmodel.collision(cdata, i_forces, j_forces);
      cdata.has_force_update = true;
    }
    else
    {
      cdata.has_force_update = false;
      cmodel.noCollision(cdata, i_forces, j_forces);
    }


    if (cdata.computeflag) {
      if (cdata.has_force_update)
        force_update(f, torque, i_forces);
    }

    if (wg->store_force_contact()) {
      wg->add_contactforce_wall(ip,i_forces);
    }

    if(wg->compute_pair_gran_local() && wg->addflag()) {
      wg->cwl_add_wall_2(cdata, i_forces);
    }
  }
};

}

}

#endif
