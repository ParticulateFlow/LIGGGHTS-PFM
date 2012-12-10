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

#include <stdlib.h>
#include <string.h>
#include "fix_mesh_surface_stress_servo.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "variable.h"
#include "error.h"
#include "vector_liggghts.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 1.0e-7

enum{NONE,CONSTANT,EQUAL,ATOM};
enum{FORCE,TORQUE,SHEAR};

/*NL*/ #define DEBUG_MESH_SURFACE_STRESS_SERVO false

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::FixMeshSurfaceStressServo(LAMMPS *lmp, int narg, char **arg) :
  FixMeshSurfaceStress(lmp, narg, arg),

  xcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes",3)),
  vcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes",1)),
  xcm_orig_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes",3)),
  vel_max_(  0.),
  int_flag_( true),
  kp_(			 0.),
  ki_(       0.),
  kd_(       0.),
  v_(        *mesh()->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_none","frame_invariant","restart_no",1))
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
    sp_str_ = NULL;

    //TEST AIGNER
    mode_flag_ = false;

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
      } else if(strcmp(arg[iarg_],"ctrlPV") == 0) {
      	  if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'ctrlPV'");
      	  iarg_++;
      	  if 			(strcmp(arg[iarg_++],"force") == 0) pv_flag_ = FORCE;
      	  else if (strcmp(arg[iarg_++],"torque") == 0) pv_flag_ = TORQUE;
      	  else if (strcmp(arg[iarg_++],"shear") == 0 ) pv_flag_ = SHEAR;
      	  else error->fix_error(FLERR,this,"only force, torque and shear are valid arguments for ctrlPV");
      	  hasargs = true;
      } else if(strcmp(arg[iarg_],"vel_max") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'vel'");
          iarg_++;
          vel_max_ = force->numeric(arg[iarg_++]);
          if(vel_max_ <= 0.)
            error->fix_error(FLERR,this,"vel_max > 0 required");
          hasargs = true;
      } else if(strcmp(arg[iarg_],"set_point") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'set_point'");
          iarg_++;
          if (strstr(arg[iarg_],"v_") == arg[iarg_]) {
        	  int n = strlen(&arg[iarg_][2]) + 1;
        	  sp_str_ = new char[n];
        	  strcpy(sp_str_,&arg[iarg_][2]);
          } else {
    				set_point_ = -force->numeric(arg[iarg_]); // the resultant force/torque/shear acts in opposite direction --> negative value
    				if (set_point_ == 0.) error->fix_error(FLERR,this,"Set point (desired force/torque/shear) has to be != 0.0");
    				set_point_inv_ = 1./set_point_;
        	  sp_style_ = CONSTANT;
          }
          iarg_++;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"dim") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments for 'forceflags'");
          iarg_++;
          if     (strcmp("x",arg[iarg_]) == 0) dim_ = 0;
          else if(strcmp("y",arg[iarg_]) == 0) dim_ = 1;
          else if(strcmp("z",arg[iarg_]) == 0) dim_ = 2;
          else error->fix_error(FLERR,this,"'x', 'y' or 'z' expected after keyword 'dim'");
          iarg_++;
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
      } else if(strcmp(arg[iarg_],"mode") == 0) { //NP Test Aigner
      	if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
      	iarg_++;
      	if (strcmp("auto",arg[iarg_]) == 0) {
      		mode_flag_ = true;
      	}	else error->fix_error(FLERR,this,"mode supports only auto");
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

    //NP inform mesh of upcoming movement
    mesh()->registerMove(false,true,false);

    // create modified andrew instance
    mod_andrew_ = new ModifiedAndrew;
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServo::~FixMeshSurfaceStressServo()
{
	delete [] sp_str_;
	delete mod_andrew_;
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
    /*NL*/ //fprintf(screen,"sizes xcm %d moi %d\n",xcm_.size(),moi_.size());

    if(!xcm_.size())
        error->fix_error(FLERR,this,"please define 'com' for the mesh");
    if(sp_style_ == CONSTANT && set_point_ == 0.)
        error->fix_error(FLERR,this,"please define 'set_point' for the mesh");
    if(vel_max_ == 0.)
        error->fix_error(FLERR,this,"please define 'vel_max' for the mesh");

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::init()
{
    FixMeshSurfaceStress::init();

    // get timestep
    reset_dt();

    // check variables
    if (sp_str_) {
    	sp_var_ = input->variable->find(sp_str_);
    	if (sp_var_ < 0)
    		error->all(FLERR,"Variable name for fix mesh/surface/stress/servo does not exist");
    	if (input->variable->equalstyle(sp_var_)) sp_style_ = EQUAL;
    	else if (input->variable->atomstyle(sp_var_)) sp_style_ = ATOM;
    	else error->all(FLERR,"Variable for fix mesh/surface/stress/servo is invalid style");
    }

    // set pointers for controller
    if (pv_flag_ == FORCE) {
    	process_value_ = &f_total_[dim_];
    	control_output_ = &vcm_(0)[dim_];
    } else if (pv_flag_ == TORQUE) {
    	process_value_ = &torque_total_[dim_];
    	control_output_ = NULL; //NP TODO: set right angular velocity
    } else error->all(FLERR,"TEST: WRONG STYLE");

    // final error checks
    if (sp_style_ == ATOM)
    	error->fix_error(FLERR,this,"Control variable of style ATOM does not make any sense for a wall");

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
	double dX[3],dx[3];

	// only if the wall should move
	if (int_flag_) {

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

	// for area calculation
	// set vector of touching particles to zero
	contacts_.clear();

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::final_integrate()
{
	//NP update forces
	FixMeshSurfaceStress::final_integrate();

	// calcualte area
	double area = mod_andrew_->area(contacts_);

	double dfdt;

	// only if the wall should move
	if (int_flag_) {

		// auto mode
		//NP TEST AIGNER
		if (mode_flag_) {

			// variable force, wrap with clear/add
			if (sp_style_ == EQUAL) {

				modify->clearstep_compute();

				set_point_ = -input->variable->compute_equal(sp_var_);
				if (set_point_ == 0.) error->fix_error(FLERR,this,"Set point (desired force/torque/shear) has to be != 0.0");
				set_point_inv_ = 1./set_point_;
				/*NL*/ //fprintf(screen,"set_point_ = %e \n",set_point_);

				modify->addstep_compute(update->ntimestep + 1);

			}

			err_ = (set_point_ - *process_value_)*set_point_inv_;

			if (abs(err_) <= 0.5) {
				*control_output_ = kp_ * err_;
			} else {
				*control_output_ = sgn(err_) * (kp_*0.5 + (2-kp_)*(abs(err_)-0.5));
			}

		} else {

			// variable force, wrap with clear/add
			if (sp_style_ == EQUAL) {

				modify->clearstep_compute();

				set_point_ = -input->variable->compute_equal(sp_var_);
				if (set_point_ == 0.) error->fix_error(FLERR,this,"Set point (desired force/torque/shear) has to be != 0.0");
				set_point_inv_ = 1./set_point_;
				/*NL*/ //fprintf(screen,"set_point_ = %e \n",set_point_);

				modify->addstep_compute(update->ntimestep + 1);

			}

			// simple PID-controller

			// calc error and sum of the errors
			err_ = (set_point_ - *process_value_);
			sum_err_ += err_*dtv_;
			// derivative term
			// force used instead of the error in order to avoid signal spikes in case of change of the set point
			// de()/dt = - dforce()/dt for constant set point
			dfdt = -( *process_value_ - old_process_value_)/dtv_;

			// vel points opposite to force vector
			*control_output_ = -vel_max_ * (err_ * kp_ + sum_err_ * ki_ + dfdt * kd_) * set_point_inv_;

			// save process value for next timestep
			old_process_value_ = *process_value_;

			/*NL*/ //fprintf(screen,"TEST: process_value_ = %g and f_total_ = %g\n",*process_value_,f_total(2));
			/*NL*/ //fprintf(screen,"TEST: control_output_ = %g and vcm_(0)[2] = %g \n",*control_output_,vcm_(0)[2]);

			/*NL*/ //fprintf(screen,"Vector: f_total = %f %f %f\n",f_total(0),f_total(1),f_total(2));
			/*NL*/ //fprintf(screen,"Vector: vcm_(0)[%d] = %g \n",i,vcm_(0)[i]);

		}

		limit_vel();
		set_v_node();

	}

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::limit_vel()
{

	double vmag, factor;
	vmag = *control_output_; //vectorMag3D(vcm_(0)); //NP TODO: Is this also okey?

	/*NL*/ //if(screen) fprintf(screen,"Test for velocity: vel_max_ = %g; vmag = %g\n",vel_max_,vmag);

	// saturation of the velocity
	if(vmag > vel_max_ && vmag != 0) {
		factor = vel_max_ / vmag;

		*control_output_ *= factor;

		// anti-windup of the integral part (only if ki_>0)
		if (ki_ > 0) {
			sum_err_ = (-sgn(*control_output_) * set_point_ -err_*kp_)/ki_; //inverted controller equation
			/*NL*/ //fprintf(screen,"vcm_(0)[%d] = %g; sum_error_vec_[%d] = %g; error_vec = %g\n",i,*control_output_,i,sum_err_,err_);
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

    //NP TODO: Torque?! Per node velocity
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

/* ----------------------------------------------------------------------
  called during wall force calc

  detected contacts are registered to contribute to the area of the servo
------------------------------------------------------------------------- */

void FixMeshSurfaceStressServo::add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall)
{
	FixMeshSurfaceStress::add_particle_contribution(ip,frc,delta,iTri,v_wall);

	double *x = atom->x[ip];
	double r = atom->radius[ip];

	Circle c = {x[(1+dim_)%3], x[(2+dim_)%3], r};
	contacts_.push_back(c);

}

