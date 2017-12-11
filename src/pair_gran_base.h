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

#ifndef PAIR_GRAN_BASE_H_
#define PAIR_GRAN_BASE_H_

#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_nonspherical.h"
#include "math_const.h"
#endif

#include "contact_interface.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "fix_contact_property_atom.h"
#include "os_specific.h"

#include "granular_pair_style.h"

namespace LIGGGHTS {
using namespace ContactModels;

namespace PairStyles {

using namespace LAMMPS_NS;

template<typename ContactModel>
class Granular : private Pointers, public IGranularPairStyle {
  CollisionData * aligned_cdata;
  ForceData * aligned_i_forces;
  ForceData * aligned_j_forces;
  ContactModel cmodel;

  inline void force_update(double * const f, double * const torque,
      const ForceData & forces) {
    for (int coord = 0; coord < 3; coord++) {
      f[coord] += forces.delta_F[coord];
      torque[coord] += forces.delta_torque[coord];
    }
  }

public:
  Granular(class LAMMPS * lmp, PairGran* parent) : Pointers(lmp),
    aligned_cdata(aligned_malloc<CollisionData>(32)),
    aligned_i_forces(aligned_malloc<ForceData>(32)),
    aligned_j_forces(aligned_malloc<ForceData>(32)),
    cmodel(lmp, parent) {
  }

  virtual ~Granular() {
    aligned_free(aligned_cdata);
    aligned_free(aligned_i_forces);
    aligned_free(aligned_j_forces);
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
        fprintf(screen, "==== PAIR SETTINGS ====\n");
        settings.print_all(screen);
        fprintf(screen, "==== PAIR SETTINGS ====\n");
      }

      if(logfile) {
        fprintf(logfile, "==== PAIR SETTINGS ====\n");
        settings.print_all(logfile);
        fprintf(logfile, "==== PAIR SETTINGS ====\n");
      }
    }
#endif

    if(!success) {
      error->all(FLERR,settings.error_message.c_str());
    }
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);

#ifdef LIGGGHTS_DEBUG
    if(comm->me == 0) {
      if(screen) {
        fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");
        force->registry.print_all(screen);
        fprintf(screen, "==== PAIR GLOBAL PROPERTIES ====\n");
      }

      if(logfile) {
        fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
        force->registry.print_all(logfile);
        fprintf(logfile, "==== PAIR GLOBAL PROPERTIES ====\n");
      }
    }
#endif
  }

  virtual void write_restart_settings(FILE * fp)
  {
    int64_t hashcode = ContactModel::STYLE_HASHCODE;
    fwrite(&hashcode, sizeof(int64_t), 1, fp);
  }

  virtual void read_restart_settings(FILE * fp)
  {
    int me = comm->me;
    int64_t hashcode = -1;
    if(me == 0){
      size_t dummy = fread(&hashcode, sizeof(int64_t), 1, fp);
      UNUSED(dummy);
      // sanity check
      if(hashcode != ContactModel::STYLE_HASHCODE)
        error->all(FLERR,"wrong pair style loaded!");
    }
  }

  double stressStrainExponent()
  {
    return cmodel.stressStrainExponent();
  }

  virtual void compute_force(PairGran * pg, int eflag, int vflag, int addflag)
  {
    if (eflag || vflag)
      pg->ev_setup(eflag, vflag);
    else
      pg->evflag = pg->vflag_fdotr = 0;

    //NP update for fix rigid done in PairGran

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    int superquadric_flag = atom->superquadric_flag;
#endif
    const int newton_pair = force->newton_pair;

    int inum = pg->list->inum;
    int * ilist = pg->list->ilist;
    int * numneigh = pg->list->numneigh;

    int ** firstneigh = pg->list->firstneigh;
    int ** firsttouch = pg->listgranhistory ? pg->listgranhistory->firstneigh : NULL;
    double ** firstshear = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

    const int dnum = pg->dnum();
    const bool store_contact_forces = pg->storeContactForces();
    const int freeze_group_bit = pg->freeze_group_bit();

    const double contactDistanceMultiplier = neighbor->contactDistanceFactor*neighbor->contactDistanceFactor;

    // clear data, just to be safe
    memset(aligned_cdata, 0, sizeof(CollisionData));
    memset(aligned_i_forces, 0, sizeof(ForceData));
    memset(aligned_j_forces, 0, sizeof(ForceData));
    aligned_cdata->area_ratio = 1.0;

    CollisionData & cdata = *aligned_cdata;
    ForceData & i_forces = *aligned_i_forces;
    ForceData & j_forces = *aligned_j_forces;
    cdata.is_wall = false;
    cdata.computeflag = pg->computeflag();
    cdata.shearupdate = pg->shearupdate();

    cmodel.beginPass(cdata, i_forces, j_forces);

    // loop over neighbors of my atoms

    for (int ii = 0; ii < inum; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const double radi = radius[i];
      int * const touch = firsttouch ? firsttouch[i] : NULL;
      double * const allshear = firstshear ? firstshear[i] : NULL;
      int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      cdata.i = i;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if (superquadric_flag) {
        cdata.radi = cbrt(0.75 * atom->volume[i] / M_PI);
      } else {
        cdata.radi = radi;
      }
#else
      cdata.radi = radi;
#endif

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double radj = radius[j];
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if (superquadric_flag) {
          cdata.radj = cbrt(0.75 * atom->volume[j] / M_PI);
        } else
          cdata.radj = radj;
#else
        cdata.radj = radj;
#endif
        const double radsum = radi + radj;

        cdata.j = j;
        cdata.delta[0] = delx;
        cdata.delta[1] = dely;
        cdata.delta[2] = delz;
        cdata.rsq = rsq;
        cdata.radsum = radsum;
        cdata.touch = touch ? &touch[jj] : NULL;
        cdata.contact_history = allshear ? &allshear[dnum*jj] : NULL;

        i_forces.reset();
        j_forces.reset();

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if (rmass) {
          cdata.mi = rmass[i];
          cdata.mj = rmass[j];
        } else {
          cdata.mi = mass[type[i]];
          cdata.mj = mass[type[j]];
        }
#endif

        cdata.v_i     = v[i];
        cdata.v_j     = v[j];
        const int itype = type[i];
        const int jtype = type[j];
        cdata.itype = itype;
        cdata.jtype = jtype;
        cdata.omega_i = omega[i];
        cdata.omega_j = omega[j];

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if (rsq < radsum * radsum && cmodel.checkSurfaceIntersect(cdata)) {
#else
        if (rsq < radsum * radsum) {
#endif
          const double r = sqrt(rsq);
          const double rinv = 1.0 / r;

          // unit normal vector
          const double enx = delx * rinv;
          const double eny = dely * rinv;
          const double enz = delz * rinv;

          // meff = effective mass of pair of particles
          // if I or J part of rigid body, use body mass
          // if I or J is frozen, meff is other particle
          double mi, mj;

          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[itype];
            mj = mass[jtype];
          }
          if (pg->fr_pair()) {
            const double * mass_rigid = pg->mr_pair();
            if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
            if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
          }

          double meff = mi * mj / (mi + mj);
          if (mask[i] & freeze_group_bit)
            meff = mj;
          if (mask[j] & freeze_group_bit)
            meff = mi;

          // copy collision data to struct (compiler can figure out a better way to
          // interleave these stores with the double calculations above.
          cdata.r = r;
          cdata.rinv = rinv;
          cdata.meff = meff;
          cdata.mi = mi;
          cdata.mj = mj;
          if (atom->sphere_flag) {
            cdata.en[0]   = enx;
            cdata.en[1]   = eny;
            cdata.en[2]   = enz;
          }

          cmodel.collision(cdata, i_forces, j_forces);

          // if there is a collision, there will always be a force
          cdata.has_force_update = true;

        } else if(rsq < contactDistanceMultiplier * radsum * radsum) {
          // apply force update only if selected contact models have requested it
          cdata.has_force_update = false;
          cmodel.noCollision(cdata, i_forces, j_forces);
        } else {
          cdata.has_force_update = false;
        }

        if(cdata.has_force_update) {
          if (cdata.computeflag) {
            force_update(f[i], torque[i], i_forces);

            if(newton_pair || j < nlocal) {
              force_update(f[j], torque[j], j_forces);
            }
          }

          //NP call to compute_pair_gran_local
          if (pg->cpl() && addflag)
            pg->cpl_add_pair(cdata, i_forces);

          if (pg->evflag)
            pg->ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0,i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2],cdata.delta[0],cdata.delta[1],cdata.delta[2]);

          if (store_contact_forces)
          {
            double forces_torques_i[6],forces_torques_j[6];

            if(pg->fix_contact_forces()->has_partner(i,atom->tag[j]) == -1)
            {
                vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
                vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
                pg->fix_contact_forces()->add_partner(i,atom->tag[j],forces_torques_i);
            }
            if(pg->fix_contact_forces()->has_partner(j,atom->tag[i]) == -1)
            {
                vectorCopy3D(j_forces.delta_F,&(forces_torques_j[0]));
                vectorCopy3D(j_forces.delta_torque,&(forces_torques_j[3]));
                pg->fix_contact_forces()->add_partner(j,atom->tag[i],forces_torques_j);
            }
          }
        }
      }
    }

    cmodel.endPass(cdata, i_forces, j_forces);

    if (pg->vflag_fdotr) {
      pg->virial_fdotr_compute();
    }

    if(store_contact_forces)
        pg->fix_contact_forces()->do_forward_comm();
  }

  virtual void compute_single_pair_force(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
  {
    cmodel.beginPass(cdata, i_forces, j_forces);

    i_forces.reset();
    j_forces.reset();

    if (cdata.rsq < cdata.radsum * cdata.radsum) {
      cmodel.collision(cdata, i_forces, j_forces);
      // if there is a collision, there will always be a force
      cdata.has_force_update = true;
    } else {
      // apply force update only if selected contact models have requested it
      cdata.has_force_update = false;
      cmodel.noCollision(cdata, i_forces, j_forces);
    }

    if(cdata.has_force_update) {
      if (cdata.computeflag) {
        force_update(atom->f[cdata.i], atom->torque[cdata.i], i_forces);

        if(force->newton_pair || cdata.j < atom->nlocal) {
          force_update(atom->f[cdata.j], atom->torque[cdata.j], j_forces);
        }
      }
    }

    cmodel.endPass(cdata, i_forces, j_forces);
  }

};

}

}
#endif /* PAIR_GRAN_BASE_H_ */
