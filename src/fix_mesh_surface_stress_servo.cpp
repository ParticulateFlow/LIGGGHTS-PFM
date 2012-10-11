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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_mesh_surface_stress_servo.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "fix_gravity.h"
#include "fix_rigid.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "mpi.h"
#include "vector_liggghts.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-7

/*NL*/ #define DEBUG_MESH_SURFACE_STRESS_SERVO false

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::FixMeshSurfaceStressServo(LAMMPS *lmp, int narg, char **arg) :
  FixMeshSurfaceStress(lmp, narg, arg),

  xcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes",3)),
  vcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes",1)),
  xcm_orig_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes",3)),
  vel_max_(  0.),
  acc_max_(  0.),
  mass_(     0.),
  f_servo_(  0.),
  fflag_(    *mesh()->prop().addGlobalProperty< VectorContainer<bool,3> > ("fflag","comm_none","frame_invariant","restart_yes",1)),
  v_(        *mesh()->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_none","frame_invariant","restart_no",1)),
  int_flag_( true),
  kp_( 		 0.)
{
    if(!trackStress())
        error->fix_error(FLERR,this,"stress = 'on' required");

    if(manipulated())
        error->warning(FLERR,"Mesh has been scaled, moved, or rotated.\n"
                             "Please note that values for 'com', 'vel' refer to the scaled, moved, or rotated configuration");

    // set defaults

    init_defaults();

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if(strcmp(arg[iarg_],"com") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'com'");
          iarg_++;
          double _com[3];
          _com[0] = force->numeric(arg[iarg_++]);
          _com[1] = force->numeric(arg[iarg_++]);
          _com[2] = force->numeric(arg[iarg_++]);
          xcm_.add(_com);
          set_p_ref(xcm_(0));
          hasargs = true;
      } else if(strcmp(arg[iarg_],"vel_max") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'vel'");
          iarg_++;
          vel_max_ = force->numeric(arg[iarg_++]);
          if(vel_max_ <= 0.)
            error->fix_error(FLERR,this,"vel_max > 0 required");
          hasargs = true;
      } else if(strcmp(arg[iarg_],"f_servo") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'f_servo'");
          iarg_++;
          f_servo_ = force->numeric(arg[iarg_++]);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"acc_max") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          acc_max_ = force->numeric(arg[iarg_++]);
          if(acc_max_ <= 0.)
            error->fix_error(FLERR,this,"acc_max > 0 required");
          if(f_servo_ == 0.)
            error->fix_error(FLERR,this,"please define 'f_servo' before 'acc_max'");
          mass_ = fabs(f_servo_)/acc_max_;
          /*NL*/ //fprintf(screen,"mass %f\n",mass_);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"dim") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'forceflags'");
          iarg_++;
          bool flags[3] = {false,false,false};
          if     (strcmp("x",arg[iarg_]) == 0) flags[0] = true;
          else if(strcmp("y",arg[iarg_]) == 0) flags[1] = true;
          else if(strcmp("z",arg[iarg_]) == 0) flags[2] = true;
          else error->fix_error(FLERR,this,"'x', 'y' or 'z' expected after keyword 'dim'");
          iarg_++;
          fflag_.add(flags);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"kp") == 0) {
    	  if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
    	  kp_ = force->numeric(arg[iarg_+1]);
    	  iarg_ = iarg_+2;
    	  hasargs = true;
      } else if(strcmp(arg[iarg_],"ki") == 0) {
    	  if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
    	  ki_ = force->numeric(arg[iarg_+1]);
    	  iarg_ = iarg_+2;
    	  hasargs = true;
      } else if(strcmp(style,"mesh/surface/stress/servo") == 0) {
          char *errmsg = new char[strlen(arg[iarg_])+20];
          sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
          error->fix_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }

    error_checks();

    // store original position
    xcm_orig_.add(xcm_(0));

    // construct servo force vector

    vectorConstruct3D
    (
        f_servo_vec_,
        fflag_(0)[0] ? f_servo_ : 0.,
        fflag_(0)[1] ? f_servo_ : 0.,
        fflag_(0)[2] ? f_servo_ : 0.
    );

    vectorNegate3D(f_servo_vec_); // my desired force points in opposite direction

    kp_ = kp_/vectorMag3D(f_servo_vec_); // to normalize the error
    ki_ = ki_/vectorMag3D(f_servo_vec_);

/*    vectorConstruct3D
    (
    	sign_servo_vec_,
    	f_servo_vec_[0] >= 0. ? 1.0 : -1.0,
    	f_servo_vec_[1] >= 0. ? 1.0 : -1.0,
    	f_servo_vec_[2] >= 0. ? 1.0 : -1.0
    );*/

    /*NL*/ //if(screen) fprintf(screen,"f_servo_vec_ = %f %f %f\n",f_servo_vec_[0],f_servo_vec_[1],f_servo_vec_[2]);
    /*NL*/ //if(screen) fprintf(screen,"sign_servo_vec_ = %f %f %f\n",sign_servo_vec_[0],sign_servo_vec_[1],sign_servo_vec_[2]);

    vectorZeroize3D(error_vec_);
    vectorZeroize3D(sum_error_vec_);
    vectorZeroize3D(old_error_vec_);

    //NP inform mesh of upcoming movement
    mesh()->registerMove(false,true,false);
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::~FixMeshSurfaceStressServo()
{
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::post_create()
{
    FixMeshSurfaceStress::post_create();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::init_defaults()
{
    double zerovec[3];
    vectorZeroize3D(zerovec);
    vcm_.add(zerovec);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::error_checks()
{
    /*NL*/ //fprintf(screen,"sizes xcm %d mass %d moi %d\n",xcm_.size(),mass_.size(),moi_.size());

    if(!xcm_.size())
        error->fix_error(FLERR,this,"please define 'com' for the mesh");
    if(f_servo_ == 0.)
        error->fix_error(FLERR,this,"please define 'f_servo' for the mesh");
    if(mass_ == 0. )
        error->fix_error(FLERR,this,"please define 'acc_max' for the mesh");
    if(vel_max_ == 0.)
        error->fix_error(FLERR,this,"please define 'vel_max' for the mesh");
    if(!fflag_.size())
        error->fix_error(FLERR,this,"please define 'dim' for the mesh");
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::init()
{
    FixMeshSurfaceStress::init();

    dtv_ = update->dt;
    dtf_ = 0.5 * update->dt * force->ftm2v;
    dtfm_ = dtf_ / mass_;

    if (strcmp(update->integrate_style,"respa") == 0)
        error->fix_error(FLERR,this,"not respa-compatible");
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceStressServo::setmask()
{
    int mask = FixMeshSurfaceStress::setmask();
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::initial_integrate(int vflag)
{
    double dX[3],dx[3];

/*    if (int_flag_) {
    	// update vcm by 1/2 step

    	if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * f_total(0);
    	if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * f_total(1);
    	if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * f_total(2);

    	limit_vel();
    	set_v_node();

    	// update xcm by full step

    	dx[0] = dtv_ * vcm_(0)[0];
    	dx[1] = dtv_ * vcm_(0)[1];
    	dx[2] = dtv_ * vcm_(0)[2];
    	vectorAdd3D(xcm_(0),dx,xcm_(0));
    	vectorSubtract3D(xcm_(0),xcm_orig_(0),dX);

    	mesh()->move(dX,dx);

    	// update reference point to COM
    	//NP would not be necessary b/c p_ref_ is moved, rotated automatically
    	//NP do it anyway to avoid long-term divergence
    	//NP which could happen b/c move, rotate is done incrementally

    	set_p_ref(xcm_(0));
    	NL //fprintf(screen,"p_ref %g %g %g\n",p_ref(0),p_ref(1),p_ref(2));
    	NL //printVec3D(screen,"xcm",xcm_(0));
    }*/

    // simple P-controller
    if (int_flag_) {

    	for (int i=0;i<3;i++) {
    		if(fflag_(0)[i]) {
    			error_vec_[i] = f_servo_vec_[i]-f_total(i);
    			if (sgn(error_vec_[i]) != sgn(old_error_vec_[i])) {
    				vectorZeroize3D(sum_error_vec_);
    				/*NL*/ //fprintf(screen,"TEST: Reset integrator at timestep %d\n",update->ntimestep);
    				/*NL*/ //fprintf(screen,"sgn new = %d sgn old = %d\n",sgn(vcm_(0)[i]),sgn(old_error_vec_[1]));
    			} else sum_error_vec_[i] += error_vec_[i];

    			vcm_(0)[i] = -vel_max_ * (error_vec_[i] * kp_ + sum_error_vec_[i] * ki_); // vel points opposite to force vector
    			old_error_vec_[i] = error_vec_[i];
    		}
    	}

    	limit_vel();
    	set_v_node();

    	/*NL*/ //printVec3D(screen,"f_servo_vec",f_servo_vec_);
    	/*NL*/ //fprintf(screen," vector f_total: %g %g %g\n",f_total(0),f_total(1),f_total(2));
    	/*NL*/ //printVec3D(screen,"vcm",vcm_(0));

    	// update xcm by full step

    	dx[0] = dtv_ * vcm_(0)[0];
    	dx[1] = dtv_ * vcm_(0)[1];
    	dx[2] = dtv_ * vcm_(0)[2];
    	vectorAdd3D(xcm_(0),dx,xcm_(0));
    	vectorSubtract3D(xcm_(0),xcm_orig_(0),dX);

    	mesh()->move(dX,dx);

    	// update reference point to COM
    	//NP would not be necessary b/c p_ref_ is moved, rotated automatically
    	//NP do it anyway to avoid long-term divergence
    	//NP which could happen b/c move, rotate is done incrementally

    	set_p_ref(xcm_(0));
    	/*NL*/ //fprintf(screen,"p_ref %g %g %g\n",p_ref(0),p_ref(1),p_ref(2));
    	/*NL*/ //printVec3D(screen,"xcm",xcm_(0));
    }
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::final_integrate()
{
    //NP update forces
    FixMeshSurfaceStress::final_integrate();

/*    if (int_flag_) {
    	// add servo force
    	//NP add same force in 3 dims, will be integrated in one dim only anyway
    	add_external_contribution(f_servo_vec_);

    	NL //double ft[3]; f_total(ft);
    	NL //printVec3D(screen,"f_total",ft);

    	// update vcm by 1/2 step

    	if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * f_total(0);
    	if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * f_total(1);
    	if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * f_total(2);
    	limit_vel();
    	set_v_node();
    }*/

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::limit_vel()
{
/*    double vmag, factor;

    double dot = f_servo_vec_[0]*f_total(0)+f_servo_vec_[1]*f_total(1)+f_servo_vec_[2]*f_total(2);
    if(dot < 0.) // f_servo_vec_ and f_total point opposite direction
    {
        //NP set velocity to 5% max vel
        factor = 0.001/f_servo_*vel_max_;
        vcm_(0)[0] = -f_servo_vec_[0] * factor;
        vcm_(0)[1] = -f_servo_vec_[1] * factor;
        vcm_(0)[2] = -f_servo_vec_[2] * factor;
    }
    else
    {
        vmag = vectorMag3D(vcm_(0));

        // re-scale velocity if needed

        if(vmag > vel_max_)
        {
            factor = vel_max_ / vmag;*/
            /*NL*/ //fprintf(screen,"factor %f\n",factor);
/*            vcm_(0)[0] *= factor;
            vcm_(0)[1] *= factor;
            vcm_(0)[2] *= factor;
        }
    }*/
	double vmag, factor;
	vmag = vectorMag3D(vcm_(0));

	if(vmag > vel_max_ && vmag != 0) {
		factor = vel_max_ / vmag;
		vcm_(0)[0] *= factor;
		vcm_(0)[1] *= factor;
		vcm_(0)[2] *= factor;
	}
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::set_v_node()
{
    int nall = mesh()->size();
    int nnodes = mesh()->numNodes();

    for(int i = 0; i < nall; i++)
        for(int j = 0; j < nnodes; j++)
            v_.set(i,j,vcm_(0));
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceStressServo::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"integrate") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");

    if (strcmp(arg[1],"start") == 0) {
    	int_flag_ = true;
    } else if (strcmp(arg[1],"stop") == 0) {
    	int_flag_ = false;
    } else
    	error->all(FLERR,"Illegal fix_modify command");

    return 2;

  }

  return 0;
}
