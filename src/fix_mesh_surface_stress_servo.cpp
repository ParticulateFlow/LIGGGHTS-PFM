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

//#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_mesh_surface_stress_servo.h"
#include "force.h"
#include "update.h"
//#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "variable.h"
//#include "memory.h"
#include "error.h"
#include "vector_liggghts.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-7

enum{NONE,CONSTANT,EQUAL,ATOM};

/*NL*/ #define DEBUG_MESH_SURFACE_STRESS_SERVO false

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::FixMeshSurfaceStressServo(LAMMPS *lmp, int narg, char **arg) :
  FixMeshSurfaceStress(lmp, narg, arg),

  xcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes",3)),
  vcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes",1)),
  xcm_orig_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes",3)),
  vel_max_(  0.),
  f_servo_(  0.),
  fflag_(    *mesh()->prop().addGlobalProperty< VectorContainer<bool,3> > ("fflag","comm_none","frame_invariant","restart_yes",1)),
  v_(        *mesh()->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_none","frame_invariant","restart_no",1)),
  int_flag_( true),
  kp_( 		 0.),
  ki_(       0.),
  kd_(       0.)
{
    if(!trackStress())
        error->fix_error(FLERR,this,"stress = 'on' required");

    if(manipulated())
        error->warning(FLERR,"Mesh has been scaled, moved, or rotated.\n"
                             "Please note that values for 'com', 'vel' refer to the scaled, moved, or rotated configuration");

    // override default from base
    size_vector = 9;

    // set defaults

    init_defaults();
    fstr = NULL;

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
          if (strstr(arg[iarg_],"v_") == arg[iarg_]) {
        	  int n = strlen(&arg[iarg_][2]) + 1;
        	  fstr = new char[n];
        	  strcpy(fstr,&arg[iarg_][2]);
          } else {
        	  f_servo_ = force->numeric(arg[iarg_]);
        	  fstyle = CONSTANT;
          }
          iarg_++;
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
      } else if(strcmp(arg[iarg_],"kd") == 0) {
    	  if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
    	  kd_ = force->numeric(arg[iarg_+1]);
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
    vectorNegate3D(f_servo_vec_); // my desired wall force points in opposite direction

    //NP TODO: Right scaling of the controller parameters?
    //kp_ = kp_/vectorMag3D(f_servo_vec_);// to normalize the error
    //ki_ = ki_/vectorMag3D(f_servo_vec_);
    //kd_ = kd_/vectorMag3D(f_servo_vec_);

    /*NL*/ //if(screen) fprintf(screen,"f_servo_vec_ = %f %f %f\n",f_servo_vec_[0],f_servo_vec_[1],f_servo_vec_[2]);

    //NP inform mesh of upcoming movement
    mesh()->registerMove(false,true,false);
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::~FixMeshSurfaceStressServo()
{
	delete [] fstr;
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

    vectorZeroize3D(error_vec_);
    vectorZeroize3D(sum_error_vec_);
    vectorZeroize3D(old_f_total_);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::error_checks()
{
    /*NL*/ //fprintf(screen,"sizes xcm %d moi %d\n",xcm_.size(),moi_.size());

    if(!xcm_.size())
        error->fix_error(FLERR,this,"please define 'com' for the mesh");
    if(fstyle == CONSTANT && f_servo_ == 0.)
        error->fix_error(FLERR,this,"please define 'f_servo' for the mesh");
    if(vel_max_ == 0.)
        error->fix_error(FLERR,this,"please define 'vel_max' for the mesh");
    if(!fflag_.size())
        error->fix_error(FLERR,this,"please define 'dim' for the mesh");
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::init()
{
    FixMeshSurfaceStress::init();

    // get timestep
    reset_dt();

    // check variables
    if (fstr) {
    	fvar = input->variable->find(fstr);
    	if (fvar < 0)
    		error->all(FLERR,"Variable name for fix mesh/surface/stress/servo does not exist");
    	if (input->variable->equalstyle(fvar)) fstyle = EQUAL;
    	else if (input->variable->atomstyle(fvar)) fstyle = ATOM;
    	else error->all(FLERR,"Variable for fix mesh/surface/stress/servo is invalid style");
    }

    // final error checks
    if (fstyle == ATOM)
    	error->fix_error(FLERR,this,"Force variable of style ATOM does not make any sense");

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

void FixMeshSurfaceStressServo::setup_pre_force(int vflag)
{
    FixMeshSurfaceStress::setup_pre_force(vflag);

    // set xcm_orig_
    xcm_orig_.set(0,xcm_(0));
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::initial_integrate(int vflag)
{
    double dX[3],dx[3],dfdt;

    // only if the wall should move
    if (int_flag_) {

    	// variable force, wrap with clear/add
    	if (fstyle == EQUAL) {

    		modify->clearstep_compute();

    		f_servo_ = input->variable->compute_equal(fvar);
    		/*NL*/ //fprintf(screen,"f_servo_ = %e \n",f_servo_);

    	    vectorConstruct3D
    	    (
    	        f_servo_vec_,
    	        fflag_(0)[0] ? f_servo_ : 0.,
    	        fflag_(0)[1] ? f_servo_ : 0.,
    	        fflag_(0)[2] ? f_servo_ : 0.
    	    );
    	    vectorNegate3D(f_servo_vec_); // my desired wall force points in opposite direction

    		modify->addstep_compute(update->ntimestep + 1);

    	}

    	// simple PID-controller
    	for (int i=0;i<3;i++) {
    		if(fflag_(0)[i]) {

    			// calc error and sum of the errors
    			error_vec_[i] = f_servo_vec_[i]-f_total(i);
    			sum_error_vec_[i] += error_vec_[i]*dtv_;
    			// derivative term
    			// force used instead of the error in order to avoid signal spikes in case of change of the set point
    			// de()/dt = - dforce()/dt for constant set point
    			dfdt = -(f_total(i)-old_f_total_[i])/dtv_;

    			// vel points opposite to force vector
    			vcm_(0)[i] = -vel_max_ * (error_vec_[i] * kp_ + sum_error_vec_[i] * ki_ + dfdt * kd_);

    			// save force for next timestep
    			old_f_total_[i] = f_total(i);

    		}
    	}

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
    	/*NL*/ //fprintf(screen,"p_ref %g %g %g\n",p_ref(0),p_ref(1),p_ref(2));
    	/*NL*/ //printVec3D(screen,"xcm",xcm_(0));
    }
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::final_integrate()
{
    //NP update forces
    FixMeshSurfaceStress::final_integrate();

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::limit_vel()
{

	double vmag, factor;
	vmag = vectorMag3D(vcm_(0));

	/*NL*/ //if(screen) fprintf(screen,"Test for velocity: vel_max_ = %g; vmag = %g\n",vel_max_,vmag);

	// saturation of the velocity
	if(vmag > vel_max_ && vmag != 0) {
		factor = vel_max_ / vmag;

    	for (int i=0;i<3;i++) {
    		if(fflag_(0)[i]) {

    			vcm_(0)[i] *= factor;

    			// anti-windup of the integral part (only if ki_>0)
    			if (ki_ > 0) {
    				sum_error_vec_[i] = (-sgn(vcm_(0)[i])-error_vec_[i]*kp_)/ki_; //inverted controller equation
    				/*NL*/ //fprintf(screen,"vcm_(0)[%d] = %g; sum_error_vec_[%d] = %g; error_vec = %g\n",i,vcm_(0)[i],i,sum_error_vec_[i],error_vec_[i]);
    			}

    		}
    	}
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

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::reset_dt()
{
  dtv_ = update->dt;
  dtf_ = 0.5 * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   return total force or torque component on body
------------------------------------------------------------------------- */

double FixMeshSurfaceStressServo::compute_vector(int n)
{
  if(n < 6) return FixMeshSurfaceStress::compute_vector(n);
  else      return xcm_(0)[n-6];
}

