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

using namespace LAMMPS_NS;

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
    n_granhistory(12);
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
    memory->destroy(rb);
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
  double dnforce[3],dtforce[3];
  double dntorque[3],dttorque[3];
  double rot;
  double A,J;

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

  if(breakmode == BREAKSTYLE_STRESS_TEMP)
  {
      if(!fix_Temp) error->all(FLERR,"Internal error in BondGran");
      Temp = fix_Temp->vector_atom;
  }

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    /*NL*/ //if (screen) fprintf(screen,"ts %d: handling id %d and %d\n",update->ntimestep,tag[i1],tag[i2]);
    type = bondlist[n][2];

    A = M_PI * rb[type] * rb[type];
    J = A * 0.5 * rb[type] * rb[type];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    rsqinv = 1./rsq;
    r = sqrt(rsq);
    rinv = 1./r;

    /*NL*/ //if (screen) fprintf(screen,"bondlist[n][3] %d, exec at ts %d\n",bondlist[n][3],update->ntimestep);
    //NP continue if bond is broken
    if(bondlist[n][3])
    {
        //NP if (screen) fprintf(screen,"i1 %d i2 %d bondlist[n][3] %d\n",i1,i2,bondlist[n][3]);
        //NP if (screen) fprintf(screen,"bondlist[n][3] %d, broken ts %d\n",bondlist[n][3],update->ntimestep);
        //NP error->all(FLERR,"broken");
        continue;
    }

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

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);

        // relative rotational velocity for torsion and bending

        wr1 = (radius[i1]*omega[i1][0] - radius[i2]*omega[i2][0]) * rinv;
        wr2 = (radius[i1]*omega[i1][1] - radius[i2]*omega[i2][1]) * rinv;
        wr3 = (radius[i1]*omega[i1][2] - radius[i2]*omega[i2][2]) * rinv;

        // normal component

        wnnr = wr1*delx + wr2*dely + wr3*delz;
        wn1 = delx*wnnr * rsqinv;
        wn2 = dely*wnnr * rsqinv;
        wn3 = delz*wnnr * rsqinv;

    //if (screen) fprintf(screen,"omega[i1] %f %f %f, omega[i2] %f %f %f, wn %f %f %f\n",omega[i1][0],omega[i1][1],omega[i1][2],omega[i2][0],omega[i2][1],omega[i2][2],wn1,wn2,wn3);

        // tangential component

        wt1 = wr1 - wn1;
        wt2 = wr2 - wn2;
        wt3 = wr3 - wn3;

    // calc change in normal forces
    dnforce[0] = - vn1 * Sn[type] * A * dt;
    dnforce[1] = - vn2 * Sn[type] * A * dt;
    dnforce[2] = - vn3 * Sn[type] * A * dt;

        // calc change in shear forces
        dtforce[0] = - vtr1 * St[type] * A * dt;
        dtforce[1] = - vtr2 * St[type] * A * dt;
        dtforce[2] = - vtr3 * St[type] * A * dt;

    // calc change in normal torque
    dntorque[0] = - wn1 * St[type] * J * dt;
    dntorque[1] = - wn2 * St[type] * J * dt;
    dntorque[2] = - wn3 * St[type] * J * dt;

    // calc change in tang torque
    dttorque[0] = - wt1 * Sn[type] * J*0.5 * dt;
    dttorque[1] = - wt2 * Sn[type] * J*0.5 * dt;
    dttorque[2] = - wt3 * Sn[type] * J*0.5 * dt;

    // rotate forces

    //rotate normal force
        rot = bondhistlist[n][0]*delx + bondhistlist[n][1]*dely + bondhistlist[n][2]*delz;
        rot *= rsqinv;
    bondhistlist[n][0] = rot*delx;
    bondhistlist[n][1] = rot*dely;
    bondhistlist[n][2] = rot*delz;

    //rotate tangential force
    rot = bondhistlist[n][3]*delx + bondhistlist[n][4]*dely + bondhistlist[n][5]*delz;
    rot *= rsqinv;
    bondhistlist[n][3] -= rot*delx;
    bondhistlist[n][4] -= rot*dely;
    bondhistlist[n][5] -= rot*delz;

    //rotate normal torque
        rot = bondhistlist[n][6]*delx + bondhistlist[n][7]*dely + bondhistlist[n][8]*delz;
        rot *= rsqinv;
    bondhistlist[n][6] = rot*delx;
    bondhistlist[n][7] = rot*dely;
    bondhistlist[n][8] = rot*delz;

    //rotate tangential torque
    rot = bondhistlist[n][9]*delx + bondhistlist[n][10]*dely + bondhistlist[n][11]*delz;
    rot *= rsqinv;
    bondhistlist[n][ 9] -= rot*delx;
    bondhistlist[n][10] -= rot*dely;
    bondhistlist[n][11] -= rot*delz;

    //increment normal and tangential force and torque
    double dissipate = 0.995;
    bondhistlist[n][0] = dissipate * bondhistlist[n][0] + dnforce[0];
    bondhistlist[n][1] = dissipate * bondhistlist[n][1] + dnforce[1];
    bondhistlist[n][2] = dissipate * bondhistlist[n][2] + dnforce[2];
    bondhistlist[n][3] = dissipate * bondhistlist[n][3] + dtforce[0];
    bondhistlist[n][4] = dissipate * bondhistlist[n][4] + dtforce[1];
    bondhistlist[n][5] = dissipate * bondhistlist[n][5] + dtforce[2];
    bondhistlist[n][6] = dissipate * bondhistlist[n][6] + dntorque[0];
    bondhistlist[n][7] = dissipate * bondhistlist[n][7] + dntorque[1];
    bondhistlist[n][8] = dissipate * bondhistlist[n][8] + dntorque[2];
    bondhistlist[n][ 9] = dissipate * bondhistlist[n][ 9] + dttorque[0];
    bondhistlist[n][10] = dissipate * bondhistlist[n][10] + dttorque[1];
    bondhistlist[n][11] = dissipate * bondhistlist[n][11] + dttorque[2];

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
        double nforce_mag = vectorMag3D(&bondhistlist[n][0]);
        double tforce_mag = vectorMag3D(&bondhistlist[n][3]);
        double ntorque_mag = vectorMag3D(&bondhistlist[n][6]);
        double ttorque_mag = vectorMag3D(&bondhistlist[n][9]);

        bool nstress = sigman_break[type] < (nforce_mag/A + 2.*ttorque_mag/J*rb[type]);
        bool tstress = tau_break[type]    < (tforce_mag/A +    ntorque_mag/J*rb[type]);
        bool toohot = false;

        if(breakmode == BREAKSTYLE_STRESS_TEMP)
        {
            toohot = 0.5 * (Temp[i1] + Temp[i2]) > T_break[type];
            /*NL*/ //if (screen) fprintf(screen,"Temp[i1] %f Temp[i2] %f, T_break[type] %f\n",Temp[i1],Temp[i2],T_break[type]);
        }

        if(nstress || tstress || toohot)
        {
            bondlist[n][3] = 1;
            /*NL*/ //if (screen) fprintf(screen,"broken bond at step %d\n",update->ntimestep);
            /*NL*/ //if(toohot && screen)fprintf(screen,"   it was too hot\n");
            /*NL*/ //if(nstress && screen)fprintf(screen,"   it was nstress\n");
            /*NL*/ //if(tstress && screen)fprintf(screen,"   it was tstress\n");
        }
    }


        //NP if (screen) fprintf(screen,"ts %d, particles %d %d - shear %f %f %f - tor %f %f %f\n",update->ntimestep,tag[i1],tag[i2],bondhistlist[n][3],bondhistlist[n][4],bondhistlist[n][5],tor1,tor2,tor3);

    // energy
    //if (eflag) error->all(FLERR,"Granular bonds currently do not support energy calculation");

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += (bondhistlist[n][0] + bondhistlist[n][3]);
      f[i1][1] += (bondhistlist[n][1] + bondhistlist[n][4]);
      f[i1][2] += (bondhistlist[n][2] + bondhistlist[n][5]);
      torque[i1][0] += radius[i1]*tor1 + (bondhistlist[n][6] + bondhistlist[n][ 9]);
      torque[i1][1] += radius[i1]*tor2 + (bondhistlist[n][7] + bondhistlist[n][10]);
      torque[i1][2] += radius[i1]*tor3 + (bondhistlist[n][8] + bondhistlist[n][11]);
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= (bondhistlist[n][0] + bondhistlist[n][3]);
      f[i2][1] -= (bondhistlist[n][1] + bondhistlist[n][4]);
      f[i2][2] -= (bondhistlist[n][2] + bondhistlist[n][5]);
      torque[i2][0] += radius[i2]*tor1 - (bondhistlist[n][6] + bondhistlist[n][ 9]);
      torque[i2][1] += radius[i2]*tor2 - (bondhistlist[n][7] + bondhistlist[n][10]);
      torque[i2][2] += radius[i2]*tor3 - (bondhistlist[n][8] + bondhistlist[n][11]);
    }

    //if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,0./*fbond*/,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondGran::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(rb,n+1,"bond:rb");
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
  if(narg < 4)  error->all(FLERR,"Incorrect args for bond coefficients");

  double rb_one = force->numeric(FLERR,arg[1]);
  double Sn_one = atof(arg[2]);
  double St_one = atof(arg[3]);

  if(Sn_one < 0. || St_one < 0.)
    error->all(FLERR,"Sn, St must be > 0 (if values > 0 were provided, they are probably too large)");

  /*NL*///if (screen) fprintf(screen,"Sn %f, St%f\n",Sn_one,St_one);

  if(force->numeric(FLERR,arg[4]) == 0. )
  {
      breakmode = BREAKSTYLE_SIMPLE;
      if (narg != 6) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else if(force->numeric(FLERR,arg[4]) == 1. )
  {
      breakmode = BREAKSTYLE_STRESS;
      if (narg != 7) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else if(force->numeric(FLERR,arg[4]) == 2. )
  {
      breakmode = BREAKSTYLE_STRESS_TEMP;
      if (narg != 8) error->all(FLERR,"Incorrect args for bond coefficients");
  }
  else  error->all(FLERR,"Incorrect args for bond coefficients");

  if (!allocated) allocate();

  double r_break_one,sigman_break_one,tau_break_one,T_break_one;
  r_break_one = sigman_break_one = tau_break_one = T_break_one = 0.0;

  if(breakmode == BREAKSTYLE_SIMPLE) r_break_one = force->numeric(FLERR,arg[5]);
  else
  {
      sigman_break_one = force->numeric(FLERR,arg[5]);
      tau_break_one = force->numeric(FLERR,arg[6]);
      if(breakmode == BREAKSTYLE_STRESS_TEMP) T_break_one = force->numeric(FLERR,arg[7]);
  }

  int ilo,ihi;
  /*NL*/ //if (screen) fprintf(screen,"atom->nbondtypes %d\n",atom->nbondtypes);
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    rb[i] = rb_one;
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
  fwrite(&rb[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&Sn[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&St[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondGran::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&rb[1],sizeof(double),atom->nbondtypes,fp);
    fread(&Sn[1],sizeof(double),atom->nbondtypes,fp);
    fread(&St[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&rb[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&Sn[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&St[1],atom->nbondtypes,MPI_DOUBLE,0,world);

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
