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
   Contributing authors for original version: Andreas Aigner (JKU)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_jkr_history.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"
#include "compute_pair_gran_local.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"

#include "limits.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define sqrtFiveOverSix 0.91287092917527685576161630466800355658790782499663875

/* ---------------------------------------------------------------------- */

PairGranJKRHistory::PairGranJKRHistory(LAMMPS *lmp) : PairGran(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum_pairgran = 3;

    Yeff = NULL;
    Geff = NULL;
    betaeff = NULL;
    veff = NULL;
    coeffRestLog = NULL;
    coeffFrict = NULL;
    coeffRollFrict = NULL;
    coeffMu = NULL;
    coeffRestMax = NULL;
    coeffStc = NULL;

    coefArea = NULL;
    coefForce1 = NULL;
    coefForce2 = NULL;

    charVelflag = 0;

    force_off = false;
    sanity_checks = true;
}

/* ---------------------------------------------------------------------- */

PairGranJKRHistory::~PairGranJKRHistory()
{
    memory->destroy(Yeff);
    memory->destroy(Geff);
    memory->destroy(betaeff);
    memory->destroy(veff);
    memory->destroy(cohEnergyDens);
    memory->destroy(coeffRestLog);
    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);

    memory->destroy(coeffMu);
    memory->destroy(coeffRestMax);
    memory->destroy(coeffStc);

    memory->destroy(coefArea);
    memory->destroy(coefForce1);
    memory->destroy(coefForce2);
}

/* ---------------------------------------------------------------------- */

void PairGranJKRHistory::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value
    args[0] = (char *) "shearx";
    args[1] = (char *) "1";
    args[2] = (char *) "sheary";
    args[3] = (char *) "1";
    args[4] = (char *) "shearz";
    args[5] = (char *) "1";
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model //NP modified C.K.
------------------------------------------------------------------------- */

inline void PairGranJKRHistory::deriveContactModelParams(int &ip, int &jp,double &meff, double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu,double &vnnr, double &area) //NP modified C.K.
{
    /*NL*/ //sprintf(testmsg,"Y=%f, v=%f, cr=%f, cf=%f\n",Y->compute_vector(0),v->compute_vector(0),coeffRest->compute_array(0,0),coeffFrict->compute_array(0,0));
 /*NP   char *testmsg=new char[200];
    sprintf(testmsg,"ip=%d, jp=%d, meff=%f, deltan=%f\n", ip,jp,meff,deltan);
    error->warning(testmsg);
    delete []testmsg;
*/
    int itype = atom->type[ip];
    int jtype = atom->type[jp];

    double Sn=2.*Yeff[itype][jtype]*area;
    double St=8.*Geff[itype][jtype]*area;

    /*NL*/ //kn=4./3.*Yeff[itype][jtype]*area; //original hertz
    kt=St;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
    gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    /*NL*/ //kn /= force->nktv2p; //original hertz
    kt /= force->nktv2p;
    /* NP
    char *testmsg=new char[200];
    sprintf(testmsg,"Yeff=%f,deltan=%f, kt=%f, gamman=%f, gammat=%f, xmu=%f\n",Yeff,deltan,kt,gamman,gammat,xmu);
    error->warning(testmsg);
    delete []testmsg;*/

    /*NL*/ //fprintf(logfile,"pair yeff %f geff %f betaeff %f\n",mpg->Yeff[itype][jtype],mpg->Geff[itype][jtype],mpg->betaeff[itype][jtype]);
    return;
}

/* ---------------------------------------------------------------------- */

void PairGranJKRHistory::compute_force(int eflag, int vflag,int addflag)
{
  //calculated from the material properties //NP modified C.K.
  double kn,kt,gamman,gammat,xmu,rmu; //NP modified C.K.

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wr_roll[3],wr_rollmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double mi,mj,meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
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
        shear = &allshear[dnum_pairgran*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        reff = radi*radj/(radsum);

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

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

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

        // area calculation //NP modified A.A.
        double area,areaOld,areaOld2,coef;
        itype = type[i];
        jtype = type[j];
        area = areaOld = reff;
        for (int k=0; k<10; k++) {
          areaOld = area;
          areaOld2 = areaOld*areaOld;
          coef = coefArea[itype][jtype]*reff*sqrt(areaOld);
          area = (areaOld*areaOld2 + areaOld*coef + deltan*areaOld*reff)/(2*areaOld2 - coef);
          if (abs(area-areaOld) < 10*__DBL_EPSILON__) break;
        }

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu,vnnr,area);         //NP modified A.A.

        // normal forces = JKR contact + normal velocity damping //NP modified A.A.

        damp = gamman*vnnr*rsqinv;  //NP modified C.K.
        //ccel = kn*(radsum-r)*rinv - damp;
        double area3;
        area3 = area*area*area;
        ccel = coefForce1[itype][jtype]*area3/reff - sqrt(coefForce2[itype][jtype]*area3);
        ccel -= damp;
        /*NL*/ //ccel = MAX(0.,ccel);//NP modified C.K. remove artificial force effect (Poeschel, Schwager: Comp. Gran. Dynamics)

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

        /*NL*///fprintf(screen,"PAIR tags %d %d, mol %d %d, \n",atom->tag[i],atom->tag[j],atom->molecule[i],atom->molecule[j]);
        /*NL*///fprintf(screen,"radii %f, %f, rsq %f, radsum*radsum %f\n",radius[i],radius[j],rsq,radsum*radsum);
        /*NL*///printVec3D(screen,"xi",x[i]);
        /*NL*///printVec3D(screen,"xj",x[j]);
        /*NL*///error->one(FLERR,"touch");

        shear = &allshear[dnum_pairgran*jj];

        /*NL*///printVec3D(screen,"shear",shear);

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

        /*NP this is the old way of doing it
        tangential forces = shear + tangential velocity damping

        fs1 = - (kt*shear[0] + gammat*vtr1); //NP modified C.K.
        fs2 = - (kt*shear[1] + gammat*vtr2); //NP modified C.K.
        fs3 = - (kt*shear[2] + gammat*vtr3); //NP modified C.K.

        rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        if (fs > fn) {
          if (shrmag != 0.0) {
            shear[0] = (fn/fs) * (shear[0] + gammat*vtr1/kt) -   //NP modified C.K.
              gammat*vtr1/kt;
            shear[1] = (fn/fs) * (shear[1] + gammat*vtr2/kt) -   //NP modified C.K.
              gammat*vtr2/kt;
            shear[2] = (fn/fs) * (shear[2] + gammat*vtr3/kt) -   //NP modified C.K.
              gammat*vtr3/kt;
            fs1 *= fn/fs;
            fs2 *= fn/fs;
            fs3 *= fn/fs;
          } else fs1 = fs2 = fs3 = 0.0;
        }*/

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

        //NP call to compute_pair_gran_local
        if(cpl && addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);

        /*NL*/ //if(update->ntimestep > 33000 && (atom->tag[i] == 1 || atom->tag[i] == 37) && (atom->tag[j] == 1 || atom->tag[j] == 37)) fprintf(screen,"contact at %d, overlap %f\n",update->ntimestep,deltan);
        /*NL*/ //fprintf(screen,"contact at step %d, force %f %f %f tags %d %d\n",update->ntimestep,fx,fy,fz,atom->tag[i],atom->tag[j]);
        /*NL*/ //error->one(FLERR,"end");

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranJKRHistory::settings(int narg, char **arg) //NP modified C.K.
{
    iarg_ = 0;

    // set defaults
    dampflag = 1;
    rollingflag = 0;
    viscousflag = 0;
    force_off = false;

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;
        if (strcmp(arg[iarg_],"force") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'force'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                force_off = false;
            else if(strcmp(arg[iarg_],"off") == 0)
                force_off = true;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'force'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"sanity_checks") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'sanity_checks'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                sanity_checks = true;
            else if(strcmp(arg[iarg_],"off") == 0)
                sanity_checks = false;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'sanity_checks'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"rolling_friction") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'rolling_friction'");
            iarg_++;
            if(strcmp(arg[iarg_],"cdt") == 0)
                rollingflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                rollingflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'cdt' or 'off' after keyword 'rolling_friction'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"tangential_damping") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'tangential_damping'");
            iarg_++;
            if(strcmp(arg[iarg_],"on") == 0)
                dampflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                dampflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'on' or 'off' after keyword 'tangential_damping'");
            iarg_++;
            hasargs = true;
        } else if (strcmp(arg[iarg_],"viscous") == 0) {
            if (narg < iarg_+2) error->all(FLERR,"Pair gran: not enough arguments for 'viscous'");
            iarg_++;
            if(strcmp(arg[iarg_],"stokes") == 0)
                viscousflag = 1;
            else if(strcmp(arg[iarg_],"off") == 0)
                viscousflag = 0;
            else
                error->all(FLERR,"Illegal pair_style gran command, expecting 'stokes' or 'off' after keyword 'viscous'");
            iarg_++;
            hasargs = true;
        } else if (force->pair_match("gran/jkr/history",1))
            error->all(FLERR,"Illegal pair_style gran command, illegal keyword");
    }

}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranJKRHistory::init_granular()
{
  int max_type = mpg->max_type();

  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties
  //NP get the mechanical properties for all models here, since other model classes are derived from
  //NP the PairGranJKRHistory class

  Y1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,force->pair_style));
  v1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,force->pair_style));

  coeffRest1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  if(viscousflag)
  {
    coeffMu1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("FluidViscosity","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    coeffRestMax1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("MaximumRestitution","property/global","peratomtypepair",max_type,max_type,force->pair_style));
    coeffStc1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("CriticalStokes","property/global","peratomtypepair",max_type,max_type,force->pair_style));
  }

  cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type,force->pair_style));

  if(charVelflag) charVel1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,force->pair_style));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          double Yi=Y1->compute_vector(i-1);
          double Yj=Y1->compute_vector(j-1);
          double vi=v1->compute_vector(i-1);
          double vj=v1->compute_vector(j-1);
          double cor = coeffRest1->compute_array(i-1,j-1);

          // error checks on Y, v, e

          if(sanity_checks)
          {
              if(strcmp(update->unit_style,"si") == 0  && Yi < 5e6)
                 error->all(FLERR,"youngsModulus >= 5e6 required for SI units");
              if(strcmp(update->unit_style,"cgs") == 0 && Yi < 5e5)
                error->all(FLERR,"youngsModulus >= 5e5 required for CGS units");

              if(vi < 0. || vi > 0.5)
                error->all(FLERR,"0 <= poissonsRatio <= 0.5 required");

              if(cor <= 0.05 || cor > 1)
                error->all(FLERR,"0.05 < coefficientRestitution <= 1 required");
          }


          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);

          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

         if(viscousflag)
         {
           coeffMu[i][j] = coeffMu1->compute_array(i-1,j-1);
           coeffRestMax[i][j] = coeffRestMax1->compute_array(i-1,j-1);
           coeffStc[i][j] = coeffStc1->compute_array(i-1,j-1);

           // error check
           if(sanity_checks)
           {
               if(coeffRestMax[i][j] <= 0. || coeffRestMax[i][j] > 1)
                 error->all(FLERR,"0 < MaximumRestitution <= 1 required");
               if(coeffMu[i][j] <= 0.)
                 error->all(FLERR,"coeffMu > 0 required");
               if(coeffStc[i][j] <= 0.)
                 error->all(FLERR,"CriticalStokes > 0 required");
           }
         }

          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          /*NL*/ //fprintf(screen,"");

          cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
          coefArea[i][j] = sqrt(M_PI*cohEnergyDens[i][j]/Yeff[i][j]);
          coefForce1[i][j] = 4/3 * Yeff[i][j];
          coefForce2[i][j] = 16*M_PI*cohEnergyDens[i][j]*Yeff[i][j];
          //omitting veff here

          /*NL*/ //fprintf(screen,"coeff rest for %i %i : %f\n",i,j,coeffRest1->compute_array(i-1,j-1));
      }
  }

  if(charVelflag) charVel = charVel1->compute_scalar();
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranJKRHistory::allocate_properties(int size)
{
    memory->destroy(Yeff);
    memory->destroy(Geff);
    memory->destroy(betaeff);
    memory->destroy(veff);
    memory->destroy(cohEnergyDens);
    memory->destroy(coeffRestLog);
    memory->destroy(coeffFrict);
    memory->destroy(coeffRollFrict);
    memory->destroy(coeffMu);
    memory->destroy(coeffRestMax);
    memory->destroy(coeffStc);

    memory->destroy(coefArea);
    memory->destroy(coefForce1);
    memory->destroy(coefForce2);

    memory->create(Yeff,size+1,size+1,"Yeff");
    memory->create(Geff,size+1,size+1,"Geff");
    memory->create(betaeff,size+1,size+1,"betaeff");
    memory->create(veff,size+1,size+1,"veff");
    memory->create(cohEnergyDens,size+1,size+1,"cohEnergyDens");
    memory->create(coeffRestLog,size+1,size+1,"coeffRestLog");
    memory->create(coeffFrict,size+1,size+1,"coeffFrict");
    memory->create(coeffRollFrict,size+1,size+1,"coeffRollFrict");
    memory->create(coeffMu,size+1,size+1,"coeffMu");
    memory->create(coeffRestMax,size+1,size+1,"coeffRestMax");
    memory->create(coeffStc,size+1,size+1,"coeffStc");

    memory->create(coefArea,size+1,size+1,"coefArea");
    memory->create(coefForce1,size+1,size+1,"coefForce1");
    memory->create(coefForce2,size+1,size+1,"coefForce2");

}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file //NP modified C.K.
------------------------------------------------------------------------- */

void PairGranJKRHistory::write_restart_settings(FILE *fp) //NP modified C.K.
{
  int writeflag = dampflag + rollingflag * 2;
  fwrite(&writeflag,sizeof(int),1,fp);
  fwrite(&viscousflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts //NP modified C.K.
------------------------------------------------------------------------- */

void PairGranJKRHistory::read_restart_settings(FILE *fp) //NP modified C.K.
{
  if (comm->me == 0) {
    int readflag;
    fread(&readflag,sizeof(int),1,fp);
    fread(&viscousflag,sizeof(int),1,fp);
    dampflag = readflag & 1;
    rollingflag = readflag & 2;
  }
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
  MPI_Bcast(&rollingflag,1,MPI_INT,0,world);
  MPI_Bcast(&viscousflag,1,MPI_INT,0,world);
}
