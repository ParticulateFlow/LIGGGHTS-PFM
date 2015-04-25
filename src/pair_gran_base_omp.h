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
#ifndef PAIR_GRAN_BASE_OMP_H_
#define PAIR_GRAN_BASE_OMP_H_

#include "contact_interface.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "thr_omp.h"
#include "suffix.h"
#include "atom.h"
#include "fix_contact_property_atom.h"

#include <omp.h>
#include <pthread.h>
#include <vector>
#include <algorithm>

// defining NDEBUG disables assertions
#include <assert.h>

namespace LIGGGHTS {
using namespace ContactModels;

namespace PairStyles {

using namespace LAMMPS_NS;

static inline void force_update(double * const f, double * const torque, const ForceData & forces) {
  for (int coord = 0; coord < 3; ++coord) {
    f[coord] += forces.delta_F[coord];
    torque[coord] += forces.delta_torque[coord];
  }
}

template<typename ContactModel>
class GranularOMP : private Pointers, public ThrOMP, public IGranularPairStyle {
  ContactModel cmodel;

  bool use_patchup_list;

  int counter;
  double time_setup;
  double time_eval;
  double time_wait_finish;
  double * time_thread_finish;
  double * time_eval_per_thread;
  double time_patchup;
  double time_total;

public:
  GranularOMP(class LAMMPS * lmp, PairGran * parent) : Pointers(lmp),
    ThrOMP(lmp, THR_PAIR), 
    cmodel(lmp, parent),
    use_patchup_list(false)
{
    counter = 0;

    time_setup = 0.0;
    time_eval = 0.0;
    time_wait_finish = 0.0;
    time_thread_finish = NULL;
    time_eval_per_thread = NULL;
    time_patchup = 0.0;
    time_total = 0.0;
  }

  ~GranularOMP()
  {
    delete [] time_eval_per_thread;
    delete [] time_thread_finish;
  }

  int64_t hashcode()
  { return cmodel.hashcode(); }

  virtual void settings(int nargs, char ** args) {
    Settings settings(Pointers::lmp);
    settings.registerOnOff("patchup", use_patchup_list);
    cmodel.registerSettings(settings);
    bool success = settings.parseArguments(nargs, args);

    if(!success) {
      error->all(FLERR,settings.error_message.c_str());
    }
  }

  virtual void init_granular() {
    cmodel.connectToProperties(force->registry);
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
      fread(&hashcode, sizeof(int64_t), 1, fp);
      // sanity check
      if(hashcode != ContactModel::STYLE_HASHCODE)
        error->all(FLERR,"wrong pair style loaded!");
    }
  }

  double stressStrainExponent()
  {
    return cmodel.stressStrainExponent();
  }

  void apply_patchup(Pair * pair, ThrData* thr) {
    std::vector<ForceUpdate> & forceUpdates = thr->patchupForceUpdates;
    std::vector<VatomUpdate> & vatomUpdates = thr->patchupVatomUpdates;

    for(std::vector<ForceUpdate>::iterator it = forceUpdates.begin(); it != forceUpdates.end(); ++it) {
      ForceUpdate & update = *it;
      force_update(update.f, update.torque, update.forces);
    }

    if (pair->evflag && pair->vflag_fdotr) {
      for(std::vector<VatomUpdate>::iterator it = vatomUpdates.begin(); it != vatomUpdates.end(); ++it) {
        VatomUpdate & update = *it;
        v_tally_patchup(pair, update.index, update.v);
      }
    }
  }

  virtual void compute_force(PairGran * pg, int eflag, int vflag, int addflag)
  {
#ifdef PAIR_OMP_TIMING
    double totalStart = MPI_Wtime();
    double startTime = totalStart;
    double endTime;
#endif

    if (eflag || vflag)
      pg->ev_setup(eflag, vflag);
    else
      pg->evflag = pg->vflag_fdotr = 0;

    const int nall = atom->nlocal + atom->nghost;
    const int nthreads = comm->nthreads;

    //NP update for fix rigid done in PairGran

#ifdef PAIR_OMP_TIMING
    if(!time_eval_per_thread) {
      time_eval_per_thread = new double[nthreads];
      for(int tid = 0; tid < nthreads; ++tid) time_eval_per_thread[tid] = 0.0;
    }
    
    if(!time_thread_finish) {
      time_thread_finish = new double[nthreads];
      for(int tid = 0; tid < nthreads; ++tid) time_thread_finish[tid] = 0.0;
    }

    endTime = MPI_Wtime();
    time_setup += endTime - startTime;
    startTime = endTime;
#endif

#ifdef PAIR_OMP_TIMING
    #pragma omp parallel default(none) shared(pg,eflag,vflag,addflag,startTime,pg)
    {
      double threadStart = startTime;
#else
    #pragma omp parallel default(none) shared(eflag,vflag,addflag,pg)
    {
#endif
      const int tid = omp_get_thread_num();
      const int iifrom = atom->thread_offsets[tid];
      const int iito   = atom->thread_offsets[tid+1];

      ThrData *thr = fix->get_thr(tid);

      if(use_patchup_list) {
        thr->reset_patchup();
      }

      if(pg->evflag)
      {
        if(fix->use_reduction()) {
          // use reductions for compute
          ev_setup_thr(eflag, vflag, nall, pg->eatom, pg->vatom, thr);

          if (force->newton_pair)
            eval<1,1,1>(pg, iifrom, iito, thr, addflag);
          else
            eval<1,0,1>(pg, iifrom, iito, thr, addflag);

          reduce_thr(this, eflag, vflag, thr);
        } else {
          // no ev_setup_thr needed because we do not use atom arrays per thread
          // only the thread local accumulators are used

          if (force->newton_pair)
            eval<1,1,0>(pg, iifrom, iito, thr, addflag);
          else
            eval<1,0,0>(pg, iifrom, iito, thr, addflag);

          // missing contributions if patchup list is used, must do this after force updates
          if (pg->vflag_fdotr && !use_patchup_list) {
            #pragma omp barrier
            if (neighbor->includegroup == 0) {
              thr->virial_fdotr_compute_omp(atom->x, atom->f, nthreads, atom->nlocal, atom->nghost, -1);
            } else {
              thr->virial_fdotr_compute_omp(atom->x, atom->f, nthreads, atom->nlocal, atom->nghost, atom->nfirst);
            }
          }
        }
      }
      else
      {
        if (force->newton_pair)
          eval<0,1,0>(pg, iifrom, iito, thr, addflag);
        else
          eval<0,0,0>(pg, iifrom, iito, thr, addflag);
      }
#ifdef PAIR_OMP_TIMING
      time_thread_finish[tid]   = MPI_Wtime();
      time_eval_per_thread[tid] += time_thread_finish[tid] - threadStart;
#endif
    } // end of omp parallel region
#ifdef PAIR_OMP_TIMING
    endTime = MPI_Wtime();
    time_wait_finish += endTime - *std::min_element(&time_thread_finish[0], &time_thread_finish[nthreads-1]);
    time_eval += endTime - startTime;
    startTime = endTime;
#endif

    // reduction of accumulators only
    // force->reduction stands for full reduction code
    // if a full reduction is used, this code is unnecessary
    if(pg->evflag && !fix->use_reduction())
    {
      for(int tid = 0; tid < nthreads; ++tid){
        ThrData * const thr = fix->get_thr(tid);
        reduce_accumulators_only(eflag, vflag, thr);
      }
    }

    if(use_patchup_list) {
      // apply remaining force updates in patchupList
      for(int tid = 0; tid < nthreads; ++tid) {
        ThrData * const thr = fix->get_thr(tid);
        apply_patchup(pg, thr);
      }

      if (pg->evflag && pg->vflag_fdotr) {
        pg->virial_fdotr_compute();
      }
    } else {
      if (pg->evflag && pg->vflag_fdotr) {
        // prevent multiple calls to update the virial
        // virial_fdotr_compute was already done in parallel region
        pg->vflag_fdotr = 0;
      }
    }
#ifdef PAIR_OMP_TIMING
    time_patchup += MPI_Wtime() - startTime;
#endif

    if(pg->storeContactForces())
      pg->fix_contact_forces()->do_forward_comm();

#ifdef PAIR_OMP_TIMING
    time_total += MPI_Wtime() - totalStart;

    if(counter < 1000) {
      ++counter;
    } else {
      counter = 0;
      printf("pair time stats:\n");
      printf("----------------\n");
      printf("time_setup: %g\n", time_setup);
      printf("time_eval: %g\n", time_eval);
      printf("time_wait_finish: %g\n", time_wait_finish);
      for(int tid = 0; tid < nthreads; ++tid) {
        printf("time_eval_per_thread[%d]: %g\n", tid, time_eval_per_thread[tid]);
      }
      if(use_patchup_list) printf("time_patchup: %g\n", time_patchup);
      printf("time_total: %g\n", time_total);
    }
#endif
  }

  template <int EVFLAG, int NEWTON_PAIR, int REDUCTION>
  void eval(PairGran * pg, const int iifrom, const int iito, ThrData * const thr, int addflag)
  {
    double * const * const x = atom->x;
    double * const * const v = atom->v;
    double * const * const omega = atom->omega;

    double * const * const f = (EVFLAG && REDUCTION) ? thr->get_f() : atom->f;
    double * const * const torque = (EVFLAG && REDUCTION) ? thr->get_torque() : atom->torque;

    const double * const radius = atom->radius;
    const double * const rmass = atom->rmass;
    const double * const mass = atom->mass;
    const int * const type = atom->type;
    const int * const mask = atom->mask;
    const int nlocal = atom->nlocal;
    const int dnum = pg->dnum();
    const bool store_contact_forces = pg->storeContactForces();
    const int freeze_group_bit = pg->freeze_group_bit();

    const int tid = thr->get_tid();
    std::vector<ForceUpdate> * updateList = use_patchup_list ?  &thr->patchupForceUpdates : NULL;

    const int * ilist = pg->list->ilist;
    const int * numneigh = pg->list->numneigh;
    const int * const * const firstneigh = pg->list->firstneigh;
    int ** firsttouch = pg->listgranhistory ? pg->listgranhistory->firstneigh : NULL;
    double ** firstshear = pg->listgranhistory ? pg->listgranhistory->firstdouble : NULL;

    FixRigid * const rigid = pg->fr_pair();

    CollisionData cdata;
    ForceData i_forces;
    ForceData j_forces;
    cdata.is_wall = false;
    cdata.computeflag = pg->computeflag();
    cdata.shearupdate = pg->shearupdate();

    cmodel.beginPass(cdata, i_forces, j_forces);

    // loop over neighbors of my atoms

    for (int ii = iifrom; ii < iito; ii++) {
      const int i = ilist[ii];
      const double xtmp = x[i][0];
      const double ytmp = x[i][1];
      const double ztmp = x[i][2];
      const double radi = radius[i];
      int * const touch = firsttouch ? firsttouch[i] : NULL;
      double * const allshear = firstshear ? firstshear[i] : NULL;
      const int * const jlist = firstneigh[i];
      const int jnum = numneigh[i];

      cdata.i = i;
      cdata.radi = radi;

      assert(ii == i); //TODO break this assumption
      assert(atom->thread[i] == tid);
      //assert(atom->in_thread_region(tid, i));

      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        
        // use patchup list instead of duplicated work in conflict case
        if(j <= i && use_patchup_list) continue; 

        const double delx = xtmp - x[j][0];
        const double dely = ytmp - x[j][1];
        const double delz = ztmp - x[j][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double radj = radius[j];
        const double radsum = radi + radj;

        cdata.j = j;
        cdata.delta[0] = delx;
        cdata.delta[1] = dely;
        cdata.delta[2] = delz;
        cdata.rsq = rsq;
        cdata.radj = radj;
        cdata.radsum = radsum;
        cdata.touch = touch ? &touch[jj] : NULL;
        cdata.contact_history = allshear ? &allshear[dnum*jj] : NULL;

        i_forces.reset();
        j_forces.reset();

        if (rsq < radsum * radsum) {
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
          const int itype = type[i];
          const int jtype = type[j];

          if (rmass) {
            mi = rmass[i];
            mj = rmass[j];
          } else {
            mi = mass[itype];
            mj = mass[jtype];
          }
          if (rigid) {
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
          cdata.itype = itype;
          cdata.jtype = jtype;
          cdata.r = r;
          cdata.rinv = rinv;
          cdata.meff = meff;
          cdata.mi = mi;
          cdata.mj = mj;
          cdata.en[0]   = enx;
          cdata.en[1]   = eny;
          cdata.en[2]   = enz;
          cdata.v_i     = v[i];
          cdata.v_j     = v[j];
          cdata.omega_i = omega[i];
          cdata.omega_j = omega[j];

          cmodel.collision(cdata, i_forces, j_forces);

          // if there is a collision, there will always be a force
          cdata.has_force_update = true;
        } else {
          // apply force update only if selected contact models have requested it
          cdata.has_force_update = false;
          cmodel.noCollision(cdata, i_forces, j_forces);
        }

        if(cdata.has_force_update) {

          if (EVFLAG && REDUCTION)
          {
            if (cdata.computeflag) {
              force_update(f[i], torque[i], i_forces);

              if(NEWTON_PAIR || j < nlocal) {
                force_update(f[j], torque[j], j_forces);
              }
            }

            ev_tally_xyz_thr(pg,i,j,nlocal,NEWTON_PAIR, 0.0,0.0,i_forces.delta_F[0],i_forces.delta_F[1],i_forces.delta_F[2],delx,dely,delz,thr);
          }
          else
          {
            // Partial Newton OFF is equivalent to same_thread = false 
            bool same_thread = atom->in_thread_region(tid, j);

            if(cdata.computeflag)
            {
              force_update(f[i], torque[i], i_forces);

              if(NEWTON_PAIR || j < nlocal)
              {
                if(same_thread) {
                  // same thread, use partial newton
                  force_update(f[j], torque[j], j_forces);
                } else if(use_patchup_list) {
                  // schedule for later execution
                  updateList->push_back(ForceUpdate(f[j], torque[j], j_forces));
                }
              }
            }

            if(EVFLAG) {
              ev_tally_xyz_omp_new(pg, i, j, nlocal, NEWTON_PAIR, i_forces.delta_F[0], i_forces.delta_F[1], i_forces.delta_F[2], delx, dely, delz, same_thread, use_patchup_list, thr);
            }

          }

          //NP call to compute_pair_gran_local
          if (pg->cpl() && addflag) {
            #pragma omp critical
            pg->cpl_add_pair(cdata, i_forces);
          }
          if (store_contact_forces)
          {
            #pragma omp critical
            {
              double forces_torques_i[6],forces_torques_j[6];

              if(!pg->fix_contact_forces()->has_partner(i,atom->tag[j]))
              {
                vectorCopy3D(i_forces.delta_F,&(forces_torques_i[0]));
                vectorCopy3D(i_forces.delta_torque,&(forces_torques_i[3]));
                pg->fix_contact_forces()->add_partner(i,atom->tag[j],forces_torques_i);
              }
              if(!pg->fix_contact_forces()->has_partner(j,atom->tag[i]))
              {
                vectorCopy3D(j_forces.delta_F,&(forces_torques_j[0]));
                vectorCopy3D(j_forces.delta_torque,&(forces_torques_j[3]));
                pg->fix_contact_forces()->add_partner(j,atom->tag[i],forces_torques_j);
              }
            }
          }
        }
      }
    }

    cmodel.endPass(cdata, i_forces, j_forces);
  }
};

}

}
#endif /* PAIR_GRAN_BASE_OMP_H_ */
