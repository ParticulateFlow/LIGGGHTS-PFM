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
   Patrick Fodor (JKU Linz)
   Christian Richter (Otto-von-Guericke-University Magdeburg)
   Matthew Schramm (Iowa State University)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "bond_gran.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "error.h"
#include "update.h"
#include "vector_liggghts.h"
#include "math_const.h"

using namespace LAMMPS_NS;

#define FLEXIBLE_BONDS
/*NP
large TODO list for granular bonds:  (could be a diploma thesis?)

+ need a better dissipative formulation than the hardcoded
  'dissipate' value which produces plastic deformation
  need some vel-dependant damping
+ need to carefully debug and validate this bond style
  valiation against fix rigid
+ check whether run this bond style w/ or w/o gran pair style active,
  (neigh_modify command)
+ need to store bond radii per particle, not per type
+ parallel debugging and testing not yet done
+ need evtally implemetation
*/

enum{
     BREAKSTYLE_SIMPLE,
     BREAKSTYLE_STRESS,
     BREAKSTYLE_STRESS_TEMP
    };

/* ---------------------------------------------------------------------- */

BondGran::BondGran(LAMMPS *lmp) : Bond(lmp)
{
    // we need 12 history values - the 6 forces and 6 torques from the last time-step
#ifdef FLEXIBLE_BONDS
    n_granhistory(13);
#else
    n_granhistory(12);
#endif
    // number of entries in bondhistlist. bondhistlist[number of bond][number of value (from 0 to number given here)]
    // so with this number you can modify how many pieces of information you savae with every bond
    // following dependencies and processes for saving,copying,growing the bondhistlist:

    // neighbor.cpp:           memory->create(bondhistlist,maxbond,atom->n_bondhist,"neigh:bondhistlist");
    // neigh_bond.cpp:         memory->grow(bondhistlist,maxbond,atom->n_bondhist,"neighbor:bondhistlist");
    // bond.cpp:               void Bond::n_granhistory(int nhist) {ngranhistory = nhist; atom->n_bondhist = ngranhistory; if(){FLERR}}
    // atom_vec_bond_gran.cpp: memory->grow(atom->bond_hist,nmax,atom->bond_per_atom,atom->n_bondhist,"atom:bond_hist");

    if(!atom->style_match("bond/gran"))
        error->all(FLERR,"A granular bond style can only be used together with atom style bond/gran");
    if(comm->me == 0)
        error->warning(FLERR,"Bond granular: This is a beta version - be careful!");
    fix_Temp = NULL;
}

/* ---------------------------------------------------------------------- */

BondGran::~BondGran()
{
  if (allocated) {
    memory->destroy(setflag);
#ifdef FLEXIBLE_BONDS
    memory->destroy(ro);
    memory->destroy(ri);
    memory->destroy(lb);
    memory->destroy(damp);
    memory->destroy(bnl);
#else
    memory->destroy(rb);
#endif
    memory->destroy(Sn);
    memory->destroy(St);
    memory->destroy(r_break);
    memory->destroy(sigman_break);
    memory->destroy(tau_break);
    memory->destroy(T_break);
  }
}

/* ---------------------------------------------------------------------- */

void  BondGran::init_style()
{
    if(breakmode == BREAKSTYLE_STRESS_TEMP)
       fix_Temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,"bond gran"));
}

/* ---------------------------------------------------------------------- */

void BondGran::compute(int eflag, int vflag)
{

  double rsq,r,rinv,rsqinv;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3;
  double wr1,wr2,wr3,vtr1,vtr2,vtr3,tor1,tor2,tor3;
  double wnnr,wn1,wn2,wn3,wt1,wt2,wt3;

  int i1,i2,n,type;
  double delx,dely,delz;
#ifndef FLEXIBLE_BONDS
  double dnforce[3];
#endif
  double dtforce[3];
  double dntorque[3],dttorque[3];
  double rot;
  double A,J;

#ifdef FLEXIBLE_BONDS
  double force_damp_n[3],force_damp_t[3];
  double torque_damp_n[3],torque_damp_t[3];
#else
  double rbmin; // parallel-bond radius
#endif

#ifdef FLEXIBLE_BONDS
  double Ip, Me, I, Js; // MS
  double *density = atom->density; //MS
  double Kn,Kt,K_ben,K_tor; //MS
  double d_fn_sqrt_2_Me_Sn, d_ft_sqrt_2_Me_St, d_mn_sqrt_2_Js_Ktor, d_mt_sqrt_2_Js_Kben; //MS
  double rin,rout,m1,m2; //MS
  double J1, J2, bondLength; //YG
  double fn_bond[3], vel_temp[3]; //YG
  double vel_norm, f_norm; //YG
#endif
  double sndt, stdt, K_tor_dt, K_ben_dt; //MS

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *radius = atom->radius;
  double **torque = atom->torque;
  double **omega = atom->omega;
  int **bondlist = neighbor->bondlist;
  double **bondhistlist = neighbor->bondhistlist;

  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;
  double dt = update->dt;
  double cutoff = neighbor->skin;

  //modified A.N.
  int xperiodic = domain->xperiodic;
  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;

  if(breakmode == BREAKSTYLE_STRESS_TEMP)
  {
      if(!fix_Temp) error->all(FLERR,"Internal error in BondGran");
      Temp = fix_Temp->vector_atom;
  }

  for (n = 0; n < nbondlist; n++) {
    /*NL*/ //if (screen) fprintf(screen,"bondlist[n][3] %d, exec at ts %d\n",bondlist[n][3],update->ntimestep);
    // continue if bond is broken
    if(bondlist[n][3])
    {
        //NP if (screen) fprintf(screen,"i1 %d i2 %d bondlist[n][3] %d\n",i1,i2,bondlist[n][3]);
        //NP if (screen) fprintf(screen,"bondlist[n][3] %d, broken ts %d\n",bondlist[n][3],update->ntimestep);
        //NP error->all(FLERR,"broken");
        continue;
    }

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    //printf("nbondlist = %d, i1 = %d, i2 = %d \n", nbondlist, i1, i2);

    // check if bond overlap the box-borders
    // consider the periodicity of the boundaries also - mod by A.N.
    if (       x[i1][0] < (domain->boxlo[0]+cutoff) && !xperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i1][0] > (domain->boxhi[0]-cutoff) && !xperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i1][1] < (domain->boxlo[1]+cutoff) && !yperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i1][1] > (domain->boxhi[1]-cutoff) && !yperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i1][2] < (domain->boxlo[2]+cutoff) && !zperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i1][2] > (domain->boxhi[2]-cutoff) && !zperiodic) {
      bondlist[n][3] = 1;
      continue;
    }

    if (       x[i2][0] < (domain->boxlo[0]+cutoff) && !xperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i2][0] > (domain->boxhi[0]-cutoff) && !xperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i2][1] < (domain->boxlo[1]+cutoff) && !yperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i2][1] > (domain->boxhi[1]-cutoff) && !yperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i2][2] < (domain->boxlo[2]+cutoff) && !zperiodic) {
      bondlist[n][3] = 1;
      continue;
    } else if (x[i2][2] > (domain->boxhi[2]-cutoff) && !zperiodic) {
      bondlist[n][3] = 1;
      continue;
    }

    /*NL*/ //if (screen) fprintf(screen,"ts %d: handling id %d and %d\n",update->ntimestep,tag[i1],tag[i2]);

    type = bondlist[n][2]; // Get current bond type properties

#ifdef FLEXIBLE_BONDS
    rin = ri[type]*MIN(radius[i1],radius[i2]);
    rout= ro[type]*MIN(radius[i1],radius[i2]);

    A = M_PI * (rout*rout - rin*rin); // area of parallel bond cross-section
    J = A * 0.5 * (rout*rout - rin*rin); // polar moment of inertia of parallel bond cross-section

    m1 = MathConst::MY_4PI3*density[i1]*radius[i1]*radius[i1]*radius[i1];
    m2 = MathConst::MY_4PI3*density[i2]*radius[i2]*radius[i2]*radius[i2];
    Me = m1*m2/(m1+m2);

    Ip = 0.5*M_PI*(rout*rout*rout*rout - rin*rin*rin*rin); // MS
    I  = 0.5*Ip;

    J1 = 0.4 * m1 * radius[i1]*radius[i1];
    J2 = 0.4 * m2 * radius[i2]*radius[i2];
    Js = J1*J2/(J1+J2);
#else
    rbmin = rb[type]*MIN(radius[i1],radius[i2]); //lamda * min(rA,rB) see Potyondy and Cundall, "A bonded-particle model for rock" (2004)

    A = M_PI * rbmin* rbmin; // area of parallel bond cross-section
    J = A * 0.5 * rbmin * rbmin; // polar moment of inertia of parallel bond cross-section
#endif

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    rsqinv = 1./rsq;
    r = sqrt(rsq);
    rinv = 1./r;

#ifdef FLEXIBLE_BONDS
    // set bond length
    bondLength = lb[type]*(radius[i1]+radius[i2]);

    // set stiffness values
    Kn = Sn[type]*A/bondLength;
    Kt = St[type]*A/bondLength;
    K_tor = St[type]*Ip/bondLength;
    K_ben = Sn[type]*I/bondLength;
#endif

    // relative translational velocity

    vr1 = v[i1][0] - v[i2][0];
    vr2 = v[i1][1] - v[i2][1];
    vr3 = v[i1][2] - v[i2][2];

    // normal component

    vnnr = vr1*delx + vr2*dely + vr3*delz;
    vn1 = delx*vnnr * rsqinv;
    vn2 = dely*vnnr * rsqinv;
    vn3 = delz*vnnr * rsqinv;

    // tangential component

    vt1 = vr1 - vn1;
    vt2 = vr2 - vn2;
    vt3 = vr3 - vn3;

    // relative rotational velocity for shear

    wr1 = (radius[i1]*omega[i1][0] + radius[i2]*omega[i2][0]) * rinv;
    wr2 = (radius[i1]*omega[i1][1] + radius[i2]*omega[i2][1]) * rinv;
    wr3 = (radius[i1]*omega[i1][2] + radius[i2]*omega[i2][2]) * rinv;

    // relative velocities for shear

    vtr1 = vt1- (delz*wr2-dely*wr3);
    vtr2 = vt2- (delx*wr3-delz*wr1);
    vtr3 = vt3 - (dely*wr1-delx*wr2);

    // relative rotational velocity for torsion and bending
#ifdef FLEXIBLE_BONDS
    wr1 = omega[i1][0] - omega[i2][0];
    wr2 = omega[i1][1] - omega[i2][1];
    wr3 = omega[i1][2] - omega[i2][2];
#else
    wr1 = (radius[i1]*omega[i1][0] - radius[i2]*omega[i2][0]) * rinv;
    wr2 = (radius[i1]*omega[i1][1] - radius[i2]*omega[i2][1]) * rinv;
    wr3 = (radius[i1]*omega[i1][2] - radius[i2]*omega[i2][2]) * rinv;
#endif

    // normal component

    wnnr =wr1*delx + wr2*dely + wr3*delz;
    wn1 = delx*wnnr * rsqinv;
    wn2 = dely*wnnr * rsqinv;
    wn3 = delz*wnnr * rsqinv;

    //if (screen) fprintf(screen,"omega[i1] %f %f %f, omega[i2] %f %f %f, wn %f %f %f\n",omega[i1][0],omega[i1][1],omega[i1][2],omega[i2][0],omega[i2][1],omega[i2][2],wn1,wn2,wn3);

    // tangential component

    wt1 = wr1 - wn1;
    wt2 = wr2 - wn2;
    wt3 = wr3 - wn3;

    // calc change in normal forces
#ifdef FLEXIBLE_BONDS
    double eps = (r-bondLength)*rinv;
    sndt = Kn * eps * exp(eps*bnl[type]);   // F = k * change in length

    fn_bond[0] = - sndt*delx;           // To get the components F = k * change in lenth * change in the given co-ordinate / new bond length
    fn_bond[1] = - sndt*dely;
    fn_bond[2] = - sndt*delz;

#else
    sndt = Sn[type] * A * dt;
    dnforce[0] = - vn1 * sndt;
    dnforce[1] = - vn2 * sndt;
    dnforce[2] = - vn3 * sndt;
#endif

    // calc change in shear forces
#ifdef FLEXIBLE_BONDS
    stdt = Kt*dt;
#else
    stdt = St[type] * A * dt;
#endif
    dtforce[0] = - vtr1 * stdt;
    dtforce[1] = - vtr2 * stdt;
    dtforce[2] = - vtr3 * stdt;

    // calc change in normal torque
#ifdef FLEXIBLE_BONDS
    K_tor_dt = K_tor*dt;
#else
    K_tor_dt = St[type] * J * dt;
#endif
    dntorque[0] = - wn1 * K_tor_dt;
    dntorque[1] = - wn2 * K_tor_dt;
    dntorque[2] = - wn3 * K_tor_dt;


    // calc change in tang torque
#ifdef FLEXIBLE_BONDS
    K_ben_dt = K_ben*dt; // K_ben will become an input parameter
#else
    K_ben_dt = Sn[type] * J*0.5 * dt;
#endif
    dttorque[0] = - wt1 * K_ben_dt;
    dttorque[1] = - wt2 * K_ben_dt;
    dttorque[2] = - wt3 * K_ben_dt;

#ifdef FLEXIBLE_BONDS
    // damping forces
    // normal force dampening
    d_fn_sqrt_2_Me_Sn = 2.0*damp[type] * sqrt(Me*Kn);
    force_damp_n[0] = d_fn_sqrt_2_Me_Sn*(-vn1);
    force_damp_n[1] = d_fn_sqrt_2_Me_Sn*(-vn2);
    force_damp_n[2] = d_fn_sqrt_2_Me_Sn*(-vn3);

    // tangential force dampening
    d_ft_sqrt_2_Me_St = 2.0*damp[type] * sqrt(Me*Kt);
    force_damp_t[0] = d_ft_sqrt_2_Me_St*(-vtr1);
    force_damp_t[1] = d_ft_sqrt_2_Me_St*(-vtr2);
    force_damp_t[2] = d_ft_sqrt_2_Me_St*(-vtr3);

    // normal moment dampening
    d_mn_sqrt_2_Js_Ktor = 2.0*damp[type] * sqrt(Js*K_tor);
    torque_damp_n[0] = d_mn_sqrt_2_Js_Ktor*(-wn1);
    torque_damp_n[1] = d_mn_sqrt_2_Js_Ktor*(-wn2);
    torque_damp_n[2] = d_mn_sqrt_2_Js_Ktor*(-wn3);

    // tangential moment dampening
    d_mt_sqrt_2_Js_Kben = 2.0*damp[type] * sqrt(Js*K_ben);
    torque_damp_t[0] = d_mt_sqrt_2_Js_Kben*(-wt1);
    torque_damp_t[1] = d_mt_sqrt_2_Js_Kben*(-wt2);
    torque_damp_t[2] = d_mt_sqrt_2_Js_Kben*(-wt3);
#endif

    // rotate forces

#ifndef FLEXIBLE_BONDS
    //rotate normal force
    rot = bondhistlist[n][0]*delx + bondhistlist[n][1]*dely + bondhistlist[n][2]*delz;
    rot *= rsqinv;
    bondhistlist[n][0] = rot*delx;
    bondhistlist[n][1] = rot*dely;
    bondhistlist[n][2] = rot*delz;
#endif

    //rotate tangential force
    rot = bondhistlist[n][3]*delx + bondhistlist[n][4]*dely + bondhistlist[n][5]*delz;
    rot *= rsqinv;
#ifdef FLEXIBLE_BONDS
    vel_temp[0] = bondhistlist[n][3] - rot*delx;
    vel_temp[1] = bondhistlist[n][4] - rot*dely;
    vel_temp[2] = bondhistlist[n][5] - rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][3]*bondhistlist[n][3] + bondhistlist[n][4]*bondhistlist[n][4] + bondhistlist[n][5]*bondhistlist[n][5];
    if (vel_norm == 0) f_norm = 0.;
    else f_norm = sqrt (f_norm) /vel_norm;

    bondhistlist[n][3] = f_norm*vel_temp[0];
    bondhistlist[n][4] = f_norm*vel_temp[1];
    bondhistlist[n][5] = f_norm*vel_temp[2];
#else
    bondhistlist[n][3] -= rot*delx;
    bondhistlist[n][4] -= rot*dely;
    bondhistlist[n][5] -= rot*delz;
#endif

    //rotate normal torque
    rot = bondhistlist[n][6]*delx + bondhistlist[n][7]*dely + bondhistlist[n][8]*delz;
    rot *= rsqinv;
#ifdef FLEXIBLE_BONDS
    vel_temp[0] = rot*delx;
    vel_temp[1] = rot*dely;
    vel_temp[2] = rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][6]*bondhistlist[n][6] + bondhistlist[n][7]*bondhistlist[n][7] + bondhistlist[n][8]*bondhistlist[n][8];
    if (vel_norm == 0) f_norm =0;
    else f_norm = sqrt(f_norm) / vel_norm;

    bondhistlist[n][6] = f_norm*vel_temp[0];
    bondhistlist[n][7] = f_norm*vel_temp[1];
    bondhistlist[n][8] = f_norm*vel_temp[2];
#else
    bondhistlist[n][6] = rot*delx;
    bondhistlist[n][7] = rot*dely;
    bondhistlist[n][8] = rot*delz;
#endif

    //rotate tangential torque
    rot = bondhistlist[n][9]*delx + bondhistlist[n][10]*dely + bondhistlist[n][11]*delz;
    rot *= rsqinv;
#ifdef FLEXIBLE_BONDS
    vel_temp[0] = bondhistlist[n][9] - rot*delx;
    vel_temp[1] = bondhistlist[n][10] - rot*dely;
    vel_temp[2] = bondhistlist[n][11] - rot*delz;
    vel_norm = sqrt (vel_temp[0]*vel_temp[0]+vel_temp[1]*vel_temp[1]+vel_temp[2]*vel_temp[2]);
    f_norm = bondhistlist[n][9]*bondhistlist[n][9] + bondhistlist[n][10]*bondhistlist[n][10] + bondhistlist[n][11]*bondhistlist[n][11];
    if (vel_norm == 0) f_norm =0;
    else f_norm = sqrt (f_norm) /vel_norm;

    bondhistlist[n][ 9] = f_norm*vel_temp[0];
    bondhistlist[n][10] = f_norm*vel_temp[1];
    bondhistlist[n][11] = f_norm*vel_temp[2];
#else
    bondhistlist[n][ 9] -= rot*delx;
    bondhistlist[n][10] -= rot*dely;
    bondhistlist[n][11] -= rot*delz;
#endif

    //increment normal and tangential force and torque
    const double dissipate = 1.0;
#ifdef FLEXIBLE_BONDS
    bondhistlist[n][0] = fn_bond[0];
    bondhistlist[n][1] = fn_bond[1];
    bondhistlist[n][2] = fn_bond[2];
#else
    bondhistlist[n][0] = dissipate * bondhistlist[n][0] + dnforce[0];
    bondhistlist[n][1] = dissipate * bondhistlist[n][1] + dnforce[1];
    bondhistlist[n][2] = dissipate * bondhistlist[n][2] + dnforce[2];
#endif
    bondhistlist[n][ 3] = bondhistlist[n][ 3] + dtforce[0];
    bondhistlist[n][ 4] = bondhistlist[n][ 4] + dtforce[1];
    bondhistlist[n][ 5] = bondhistlist[n][ 5] + dtforce[2];
    bondhistlist[n][ 6] = bondhistlist[n][ 6] + dntorque[0];
    bondhistlist[n][ 7] = bondhistlist[n][ 7] + dntorque[1];
    bondhistlist[n][ 8] = bondhistlist[n][ 8] + dntorque[2];
    bondhistlist[n][ 9] = bondhistlist[n][ 9] + dttorque[0];
    bondhistlist[n][10] = bondhistlist[n][10] + dttorque[1];
    bondhistlist[n][11] = bondhistlist[n][11] + dttorque[2];

    //torque due to tangential bond force
    tor1 = - rinv * (dely*bondhistlist[n][5] - delz*bondhistlist[n][4]);
    tor2 = - rinv * (delz*bondhistlist[n][3] - delx*bondhistlist[n][5]);
    tor3 = - rinv * (delx*bondhistlist[n][4] - dely*bondhistlist[n][3]);

    //flag breaking of bond if criterion met
    if(breakmode == BREAKSTYLE_SIMPLE)
    {
        if(r > 2. * r_break[type])
        {
            /*NL*///if (screen) fprintf(screen,"step " BIGINT_FORMAT " broke bond between atom tags %d %d r %f, 2. * r_break[type] %f \n",
            /*NL*///         update->ntimestep,atom->tag[i1],atom->tag[i2],r,2. * r_break[type]);
            bondlist[n][3] = 1;
            //NP error->all(FLERR,"broken");
        }
    }
    else //NP stress or stress_temp
    {
#ifdef FLEXIBLE_BONDS
        double nforce_mag  = sqrt(fn_bond[0]*fn_bond[0] + fn_bond[1]*fn_bond[1] + fn_bond[2]*fn_bond[2]);
#else
        double nforce_mag  = vectorMag3D(&bondhistlist[n][0]);
#endif
        double tforce_mag  = vectorMag3D(&bondhistlist[n][3]);
        double ntorque_mag = vectorMag3D(&bondhistlist[n][6]);
        double ttorque_mag = vectorMag3D(&bondhistlist[n][9]);

#ifdef FLEXIBLE_BONDS
        bool nstress = sigman_break[type] < (nforce_mag/A + 2.*ttorque_mag/J*(rout-rin));
        bool tstress = tau_break[type]    < (tforce_mag/A +    ntorque_mag/J*(rout-rin));

        //printf("Sigma = %f, Tau = %f \n",(nforce_mag/A + 2.*ttorque_mag/J*(rout-rin)),(tforce_mag/A +    ntorque_mag/J*(rout-rin)));
#else
        bool nstress = sigman_break[type] < (nforce_mag/A + 2.*ttorque_mag/J*rbmin);
        bool tstress = tau_break[type]    < (tforce_mag/A +    ntorque_mag/J*rbmin);
#endif
        bool toohot = false;

        if(breakmode == BREAKSTYLE_STRESS_TEMP)
        {
            toohot = 0.5 * (Temp[i1] + Temp[i2]) > T_break[type];
            /*NL*/ //if (screen) fprintf(screen,"Temp[i1] %f Temp[i2] %f, T_break[type] %f\n",Temp[i1],Temp[i2],T_break[type]);
        }

        if(nstress || tstress || toohot)
        {
            bondlist[n][3] = 1;
            if (screen) fprintf(screen,"broken bond between atoms %d and %d at time %ld \n",atom->tag[i1],atom->tag[i2],update->ntimestep);
            /*NL*/ //if(toohot && screen)fprintf(screen,"   it was too hot\n");
            if(nstress && screen)fprintf(screen,"   it was nstress\n");
            if(tstress && screen)fprintf(screen,"   it was tstress\n");
        }
    }

    //NP if (screen) fprintf(screen,"ts %d, particles %d %d - shear %f %f %f - tor %f %f %f\n",update->ntimestep,tag[i1],tag[i2],bondhistlist[n][3],bondhistlist[n][4],bondhistlist[n][5],tor1,tor2,tor3);

    // energy
    //if (eflag) error->all(FLERR,"Granular bonds currently do not support energy calculation");

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
#ifdef FLEXIBLE_BONDS
      f[i1][0] += (fn_bond[0] + bondhistlist[n][3]) + (force_damp_n[0] + force_damp_t[0]);
      f[i1][1] += (fn_bond[1] + bondhistlist[n][4]) + (force_damp_n[1] + force_damp_t[1]);
      f[i1][2] += (fn_bond[2] + bondhistlist[n][5]) + (force_damp_n[2] + force_damp_t[2]);

      torque[i1][0] += radius[i1]*tor1 + (bondhistlist[n][6] + bondhistlist[n][ 9])+(torque_damp_n[0]+torque_damp_t[0]);
      torque[i1][1] += radius[i1]*tor2 + (bondhistlist[n][7] + bondhistlist[n][10])+(torque_damp_n[1]+torque_damp_t[1]);
      torque[i1][2] += radius[i1]*tor3 + (bondhistlist[n][8] + bondhistlist[n][11])+(torque_damp_n[2]+torque_damp_t[2]);
#else
      f[i1][0] += (bondhistlist[n][0] + bondhistlist[n][3]);
      f[i1][1] += (bondhistlist[n][1] + bondhistlist[n][4]);
      f[i1][2] += (bondhistlist[n][2] + bondhistlist[n][5]);

      torque[i1][0] += radius[i1]*tor1 + (bondhistlist[n][6] + bondhistlist[n][ 9]);
      torque[i1][1] += radius[i1]*tor2 + (bondhistlist[n][7] + bondhistlist[n][10]);
      torque[i1][2] += radius[i1]*tor3 + (bondhistlist[n][8] + bondhistlist[n][11]);
#endif
    }

    if (newton_bond || i2 < nlocal) {
#ifdef FLEXIBLE_BONDS
      f[i2][0] -= (fn_bond[0] + bondhistlist[n][3]) + (force_damp_n[0] + force_damp_t[0]);
      f[i2][1] -= (fn_bond[1] + bondhistlist[n][4]) + (force_damp_n[1] + force_damp_t[1]);
      f[i2][2] -= (fn_bond[2] + bondhistlist[n][5]) + (force_damp_n[2] + force_damp_t[2]);

      torque[i2][0] += radius[i2]*tor1 - (bondhistlist[n][6]+bondhistlist[n][ 9]) - (torque_damp_n[0]+torque_damp_t[0]);
      torque[i2][1] += radius[i2]*tor2 - (bondhistlist[n][7]+bondhistlist[n][10]) - (torque_damp_n[1]+torque_damp_t[1]);
      torque[i2][2] += radius[i2]*tor3 - (bondhistlist[n][8]+bondhistlist[n][11]) - (torque_damp_n[2]+torque_damp_t[2]);
#else
      f[i2][0] -= (bondhistlist[n][0] + bondhistlist[n][3]);
      f[i2][1] -= (bondhistlist[n][1] + bondhistlist[n][4]);
      f[i2][2] -= (bondhistlist[n][2] + bondhistlist[n][5]);

      torque[i2][0] += radius[i2]*tor1 - (bondhistlist[n][6] + bondhistlist[n][ 9]);
      torque[i2][1] += radius[i2]*tor2 - (bondhistlist[n][7] + bondhistlist[n][10]);
      torque[i2][2] += radius[i2]*tor3 - (bondhistlist[n][8] + bondhistlist[n][11]);
#endif
    }

    //if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,0./*fbond*/,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondGran::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

#ifdef FLEXIBLE_BONDS
  memory->create(ro,n+1,"bond:ro");
  memory->create(ri,n+1,"bond:ri");
  memory->create(lb,n+1,"bond:lb");
  memory->create(damp,n+1,"bond:damp");
  memory->create(bnl,n+1,"bond:bnl");
#else
  memory->create(rb,n+1,"bond:rb");
#endif
  memory->create(Sn,n+1,"bond:Sn");
  memory->create(St,n+1,"bond:St");

  memory->create(r_break,(n+1),"bond:r_break");
  memory->create(sigman_break,(n+1),"bond:sigman_break");
  memory->create(tau_break,(n+1),"bond:tau_break");
  memory->create(T_break,(n+1),"bond:T_break");

  memory->create(setflag,(n+1),"bond:setflag");
  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondGran::coeff(int narg, char **arg)
{
#ifdef FLEXIBLE_BONDS
  if(narg < 8) error->all(FLERR,"Incorrect args for bond coefficients (ro, ri, lb, sn, st, damp, bnl, b_density)"); // Matt Schramm
#else
  if(narg < 4) error->all(FLERR,"Incorrect args for bond coefficients");
#endif

#ifdef FLEXIBLE_BONDS
  double ro_one = force->numeric(FLERR,arg[1]);
  double ri_one = force->numeric(FLERR,arg[2]);
  double lb_one = force->numeric(FLERR,arg[3]);
  double Sn_one = force->numeric(FLERR,arg[4]);
  double St_one = force->numeric(FLERR,arg[5]);
  double damp_one = force->numeric(FLERR,arg[6]);
  double bnl_one = force->numeric(FLERR,arg[7]);
#else
  double rb_one = force->numeric(FLERR,arg[1]);
  double Sn_one = force->numeric(FLERR,arg[2]);
  double St_one = force->numeric(FLERR,arg[3]);
#endif

#ifdef FLEXIBLE_BONDS
  if (ro_one <= ri_one)
    error->all(FLERR,"ro must be greater than ri");
#endif

  if(Sn_one < 0. || St_one < 0.)
    error->all(FLERR,"Sn, St must be > 0 (if values > 0 were provided, they are probably too large)");

  /*NL*///if (screen) fprintf(screen,"Sn %f, St%f\n",Sn_one,St_one);

#ifdef FLEXIBLE_BONDS
  int iarg = 8;
#else
  int iarg = 4;
#endif

  if(force->numeric(FLERR,arg[iarg]) == 0. )
  {
      breakmode = BREAKSTYLE_SIMPLE;
      if (narg != iarg+2) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else if(force->numeric(FLERR,arg[iarg]) == 1. )
  {
      breakmode = BREAKSTYLE_STRESS;
      if (narg != iarg+3) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else if(force->numeric(FLERR,arg[iarg]) == 2. )
  {
      breakmode = BREAKSTYLE_STRESS_TEMP;
      if (narg != iarg+4) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else  error->all(FLERR,"Incorrect args for bond coefficients");

  if (!allocated) allocate();

  double r_break_one,sigman_break_one,tau_break_one,T_break_one;
  r_break_one = sigman_break_one = tau_break_one = T_break_one = 0.0;

  if(breakmode == BREAKSTYLE_SIMPLE) r_break_one = force->numeric(FLERR,arg[iarg+1]);
  else
  {
      sigman_break_one = force->numeric(FLERR,arg[iarg+1]);
      tau_break_one = force->numeric(FLERR,arg[iarg+2]);
      if(breakmode == BREAKSTYLE_STRESS_TEMP) T_break_one = force->numeric(FLERR,arg[iarg+3]);
  }

  int ilo,ihi;
  /*NL*/ //if (screen) fprintf(screen,"atom->nbondtypes %d\n",atom->nbondtypes);
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
#ifdef FLEXIBLE_BONDS
    ro[i] = ro_one;
    ri[i] = ri_one;
    lb[i] = lb_one;
    damp[i] = damp_one;
    bnl[i] = bnl_one;
#else
    rb[i] = rb_one;
#endif
    Sn[i] = Sn_one;
    St[i] = St_one;

    if(breakmode == BREAKSTYLE_SIMPLE) r_break[i] = r_break_one;
    else
    {
        sigman_break[i] = sigman_break_one;
        tau_break[i] = tau_break_one;
        if(breakmode == BREAKSTYLE_STRESS_TEMP) T_break[i] = T_break_one;
    }
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients - or the bonds are not initialized in create_atoms");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondGran::equilibrium_distance(int i)
{
  //NP ATTENTION: this is _not_ correct - and has implications on fix shake, pair_lj_cut_coul_long and pppm
  //NP it is not possible to define a general equilibrium distance for this bond model
  //NP as rotational degree of freedom is present
  return 0.;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondGran::write_restart(FILE *fp)
{
#ifdef FLEXIBLE_BONDS
  fwrite(&ro[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&ri[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&lb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&bnl[1],sizeof(double),atom->nbondtypes,fp);
#else
  fwrite(&rb[1],sizeof(double),atom->nbondtypes,fp);
#endif
  fwrite(&Sn[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&St[1],sizeof(double),atom->nbondtypes,fp);
#ifdef FLEXIBLE_BONDS
  fwrite(&damp[1],sizeof(double),atom->nbondtypes,fp);
#endif
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGran::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
#ifdef FLEXIBLE_BONDS
    fread(&ro[1],sizeof(double),atom->nbondtypes,fp);
    fread(&ri[1],sizeof(double),atom->nbondtypes,fp);
    fread(&lb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&bnl[1],sizeof(double),atom->nbondtypes,fp);
#else
    fread(&rb[1],sizeof(double),atom->nbondtypes,fp);
#endif
    fread(&Sn[1],sizeof(double),atom->nbondtypes,fp);
    fread(&St[1],sizeof(double),atom->nbondtypes,fp);
#ifdef FLEXIBLE_BONDS
    fread(&damp[1],sizeof(double),atom->nbondtypes,fp); //MS
#endif
  }
#ifdef FLEXIBLE_BONDS
  MPI_Bcast(&ro[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ri[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&lb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bnl[1],atom->nbondtypes,MPI_DOUBLE,0,world);
#else
  MPI_Bcast(&rb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
#endif
  MPI_Bcast(&Sn[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&St[1],atom->nbondtypes,MPI_DOUBLE,0,world);
#ifdef FLEXIBLE_BONDS
  MPI_Bcast(&damp[1],atom->nbondtypes,MPI_DOUBLE,0,world); //MS
#endif

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondGran::write_data(FILE *fp)
{
  //NP needs revision
  error->all(FLERR,"Bond granular does not support this feature");
  //for (int i = 1; i <= atom->nbondtypes; i++)
  //  fprintf(fp,"%d %g %g %g %g %g\n",i,k[i],b1[i],b2[i],rc[i],u0[i]);
}

/* ---------------------------------------------------------------------- */

double BondGran::single(int type, double rsq, int i, int j,double &fforce)
{
  error->all(FLERR,"Bond granular does not support this feature");
  /*double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = ????
  return rk*dr;*/
  return 0.;
}
