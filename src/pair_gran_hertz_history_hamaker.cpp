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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pair_gran_hertz_history_hamaker.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"
#include "neigh_list.h"
#include "compute_pair_gran_local.h"
#include "fix_rigid.h"
#include "mech_param_gran.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;

// #define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875
// #define MAX(a,b) ((a) > (b) ? (a) : (b))


/* ---------------------------------------------------------------------- */

PairGranHertzHistoryHamaker::PairGranHertzHistoryHamaker(LAMMPS *lmp) :
  PairGranHertzHistory(lmp)
{
  aHamakerEff = NULL;
  hCutEff = NULL;
  hMaxEff = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHertzHistoryHamaker::~PairGranHertzHistoryHamaker()
{
    memory->destroy(aHamakerEff);
    memory->destroy(hCutEff);
    memory->destroy(hMaxEff);
}

/* ---------------------------------------------------------------------- */

void PairGranHertzHistoryHamaker::compute_force(int eflag, int vflag,int addflag)
{
  //calculated from the material properties
  double kn,kt,gamman,gammat,xmu,rmu;
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = type[j];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

        // unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[dnum()*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

        r = sqrt(rsq);
        rinv = 1.0/r;

        Fn_coh = addCohesionForce(i,j,r);
        ccel = - Fn_coh*rinv;

        // forces & torques

        fx = delx*ccel;
        fy = dely*ccel;
        fz = delz*ccel;

        if(computeflag)
        {
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
        }

        if (j < nlocal && computeflag) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
        }
        // TODO: also for cohesion?!; NULL?!
        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,0.0,0.0,0.0,NULL);

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);


      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        double deltan=radsum-r;
        cri = radi-0.5*deltan;
        crj = radj-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        double mi,mj;
        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          itype = type[i];
          jtype = type[j];
          mi = mass[itype];
          mj = mass[jtype];
        }
        if (fix_rigid)
        {
          if(body[i] >= 0) mi = masstotal[body[i]];
          if(body[j] >= 0) mj = masstotal[body[j]];
        }

        meff = mi*mj/(mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr);   //modified C.K

        damp = gamman*vnnr*rsqinv;
        ccel = kn*(radsum-r)*rinv - damp;

        Fn_coh = addCohesionForce(i,j,0.0); // r = 0.0, because particles are in contact
        ccel -= Fn_coh*rinv;

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

        shear = &allshear[dnum()*jj];

        if (shearupdate && computeflag)
        {
          shear[0] += vtr1*dt;
          shear[1] += vtr2*dt;
          shear[2] += vtr3*dt;

          // rotate shear displacements

          rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
          rsht *= rsqinv;
          shear[0] -= rsht*delx;
          shear[1] -= rsht*dely;
          shear[2] -= rsht*delz;
        }

        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +  shear[2]*shear[2]);

        // tangential forces = shear + tangential velocity damping

        fs1 = - (kt*shear[0]);
        fs2 = - (kt*shear[1]);
        fs3 = - (kt*shear[2]);

        // rescale frictional displacements and forces if needed
        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        // energy loss from sliding or damping
        if (fs > fn) {
          if (shrmag != 0.0) {
            fs1 *= fn/fs;
            fs2 *= fn/fs;
            fs3 *= fn/fs;
            shear[0] = -fs1/kt;
            shear[1] = -fs2/kt;
            shear[2] = -fs3/kt;
          }
          else fs1 = fs2 = fs3 = 0.0;
        }
        else
        {
          fs1 -= (gammat*vtr1);
          fs2 -= (gammat*vtr2);
          fs3 -= (gammat*vtr3);
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);

        // add rolling friction torque
        vectorZeroize3D(r_torque);
        if(rollingflag)
        {
          vectorSubtract3D(omega[i],omega[j],wr_roll);
          wr_rollmag = vectorMag3D(wr_roll);

          if(wr_rollmag > 0.)
          {
            // calculate torque
            reff=radi*radj/(radi+radj);
            vectorScalarMult3D(wr_roll,rmu*kn*deltan*reff/wr_rollmag,r_torque);

            // remove normal (torsion) part of torque
            double rtorque_dot_delta = r_torque[0]*delx + r_torque[1]*dely + r_torque[2]*delz;
            r_torque_n[0] = delx * rtorque_dot_delta * rsqinv;
            r_torque_n[1] = dely * rtorque_dot_delta * rsqinv;
            r_torque_n[2] = delz * rtorque_dot_delta * rsqinv;
            vectorSubtract3D(r_torque,r_torque_n,r_torque);
          }
        }

        if(computeflag)
        {
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
          torque[i][0] -= cri*tor1 + r_torque[0];
          torque[i][1] -= cri*tor2 + r_torque[1];
          torque[i][2] -= cri*tor3 + r_torque[2];
        }

        if (j < nlocal && computeflag) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= crj*tor1 - r_torque[0];
          torque[j][1] -= crj*tor2 - r_torque[1];
          torque[j][2] -= crj*tor3 - r_torque[2];
        }

        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

}

/* ----------------------------------------------------------------------
   insert Hamaker model for van der Waals forces
------------------------------------------------------------------------- */

inline double PairGranHertzHistoryHamaker::addCohesionForce(int ip, int jp,double r)
{
  double Fn_coh;
  //r is the distance between the sphere's centers
  double ri = atom->radius[ip];
  double rj = atom->radius[jp];
  int itype = atom->type[ip];
  int jtype = atom->type[jp];
  double hIJ = (r-ri-rj);

  if (hIJ > hCutEff[itype][jtype])
    Fn_coh = aHamakerEff[itype][jtype]/(6*hIJ*hIJ)*(ri*rj/(ri+rj));
  else
    Fn_coh = aHamakerEff[itype][jtype]/(6*hCutEff[itype][jtype]*hCutEff[itype][jtype])*(ri*rj/(ri+rj));

  return Fn_coh;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHertzHistoryHamaker::settings(int narg, char **arg)
{

  PairGranHertzHistory::settings(narg, arg);

  iarg_ = 0;

  // set defaults
  cohesionflag = 1;

  // parse args

  bool hasargs = true;
  while(iarg_ < narg && hasargs)
  {
    hasargs = false;
    if (force->pair_match("gran/hooke/history",1) || force->pair_match("gran/hertz/history",1))
      error->all(FLERR,"Illegal pair_style gran command, illegal keyword");
  }

  if(cohesionflag && domain->dimension!=3)
    error->all(FLERR,"Cohesion model valid for 3d simulations only");

}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHertzHistoryHamaker::init_granular()
{
  double skin = neighbor->skin;
  double maxhMaxEff = -1.0;
  int max_type = mpg->max_type();

  PairGranHertzHistory::init_granular();

  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties

  aHamaker_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("hamakerConstant","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  hCut_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("minParticleDist","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
    for(int j=1;j<max_type+1;j++)
    {

      aHamakerEff[i][j] = aHamaker_->compute_array(i-1,j-1);
      hCutEff[i][j] = hCut_->compute_array(i-1,j-1);
      hMaxEff[i][j] = hCutEff[i][j]*100; // force f(hMax)/f(hCut) = 0.0001

      maxhMaxEff = MAX(maxhMaxEff,hMaxEff[i][j]); // get maximum hMaxEff for skin check


    }
  }

  // check skin
  if(skin < maxhMaxEff) {
    if(screen) fprintf(screen,"Maximum cutoff distance (~ minParticleDist) = %f. Skin = %f\n",maxhMaxEff,skin);
    error->all(FLERR,"Skin is too small for Hamaker model.\n");
  }
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHertzHistoryHamaker::allocate_properties(int size)
{
  memory->destroy(aHamakerEff);
  memory->destroy(hCutEff);
  memory->destroy(hMaxEff);
  memory->create(aHamakerEff,size+1,size+1,"aHamakerEff");
  memory->create(hCutEff,size+1,size+1,"hCutEff");
  memory->create(hMaxEff,size+1,size+1,"hMaxEff");
}

/* ---------------------------------------------------------------------- */
