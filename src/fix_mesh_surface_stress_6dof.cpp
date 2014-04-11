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
#include "fix_mesh_surface_stress_6dof.h"
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

/*NL*/ #define DEBUG_MESH_SURFACE_STRESS_6DOF false

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStress6DOF::FixMeshSurfaceStress6DOF(LAMMPS *lmp, int narg, char **arg) :
  FixMeshSurfaceStress(lmp, narg, arg),

  xcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm","comm_none","frame_invariant","restart_yes",3)),
  quat_(     *mesh()->prop().addGlobalProperty< VectorContainer<double,4> > ("quat","comm_none","frame_invariant","restart_yes",1)),
  vcm_(      *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("vcm","comm_none","frame_invariant","restart_yes",1)),
  omega_(    *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("omega","comm_none","frame_invariant","restart_yes",1)),
  angmom_(   *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("angmom","comm_none","frame_invariant","restart_yes",1)),
  mass_(     *mesh()->prop().addGlobalProperty< ScalarContainer<double> >   ("mass","comm_none","frame_invariant","restart_yes",3)),
  moi_(      *mesh()->prop().addGlobalProperty< MultiVectorContainer<double,3,3> > ("moi","comm_none","frame_invariant","restart_yes",2)),
  inertia_(  *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("inertia","comm_none","frame_invariant","restart_yes",2)),
  ex_space_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("ex_space","comm_none","frame_invariant","restart_yes",2)),
  ey_space_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("ey_space","comm_none","frame_invariant","restart_yes",2)),
  ez_space_( *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("ez_space","comm_none","frame_invariant","restart_yes",2)),

  xcm_orig_(  *mesh()->prop().addGlobalProperty< VectorContainer<double,3> > ("xcm_orig","comm_none","frame_invariant","restart_yes",3)),
  quat_orig_( *mesh()->prop().addGlobalProperty< VectorContainer<double,4> > ("quat_orig","comm_none","frame_invariant","restart_yes",1)),

  fflag_(    *mesh()->prop().addGlobalProperty< VectorContainer<bool,3> > ("fflag","comm_none","frame_invariant","restart_yes",1)),
  tflag_(    *mesh()->prop().addGlobalProperty< VectorContainer<bool,3> > ("tflag","comm_none","frame_invariant","restart_yes",1)),

  displace_( *mesh()->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("displace","comm_none","frame_invariant","restart_yes",3)),

  v_(        *mesh()->prop().addElementProperty< MultiVectorContainer<double,3,3> > ("v","comm_exchange_borders","frame_invariant","restart_no",1)),

  suspension_flag_(false),
  k_t_(0.),
  c_t_(0.),
  k_r_(0.),
  c_r_(0.),
  rot_flip_flag_(false),
  rot_flip_angle_(0.)
{
    if(!trackStress())
        error->fix_error(FLERR,this,"stress = 'on' required");

    if(manipulated())
        error->warning(FLERR,"Mesh has been scaled, moved, or rotated.\n"
                             "Please note that values for 'com', 'vel', 'mass', 'moi' refer to the scaled, moved, or rotated configuration");

    // set defaults

    init_defaults();

    // parse further args
    //NP may NOT change quat since MultiNodeMesh::setRotation() assumes is called with
    //NP unit quat the first time

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if(strcmp(arg[iarg_],"com") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'com'");
          iarg_++;
          double _com[3];
          _com[0] = force->numeric(FLERR,arg[iarg_++]);
          _com[1] = force->numeric(FLERR,arg[iarg_++]);
          _com[2] = force->numeric(FLERR,arg[iarg_++]);
          xcm_.add(_com);
          set_p_ref(xcm_(0));
          hasargs = true;
      } else if(strcmp(arg[iarg_],"vel") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'vel'");
          iarg_++;
          double _vel[3];
          _vel[0] = force->numeric(FLERR,arg[iarg_++]);
          _vel[1] = force->numeric(FLERR,arg[iarg_++]);
          _vel[2] = force->numeric(FLERR,arg[iarg_++]);
          vcm_.set(0,_vel);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"angmom") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'angmom'");
          iarg_++;
          double _angmom[3];
          _angmom[0] = force->numeric(FLERR,arg[iarg_++]);
          _angmom[1] = force->numeric(FLERR,arg[iarg_++]);
          _angmom[2] = force->numeric(FLERR,arg[iarg_++]);
          angmom_.set(0,_angmom);
          hasargs = true;
      } else if (strcmp(arg[iarg_],"mass") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          double _mass = force->numeric(FLERR,arg[iarg_++]);
          if(_mass <= 0.)
            error->fix_error(FLERR,this,"mass > 0 required");
          mass_.add(_mass);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"moi") == 0) {
          if (narg < iarg_+7) error->fix_error(FLERR,this,"not enough arguments for 'moi'");
          iarg_++;
          double **_moi;
          memory->create<double>(_moi,3,3,"6dof:_moi");
          _moi[0][0] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][1] = force->numeric(FLERR,arg[iarg_++]);
          _moi[2][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[0][1] = force->numeric(FLERR,arg[iarg_++]);
          _moi[0][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][2] = force->numeric(FLERR,arg[iarg_++]);
          _moi[1][0] = _moi[0][1];
          _moi[2][0] = _moi[0][2];
          _moi[2][1] = _moi[1][2];
          moi_.add(_moi);
          memory->destroy<double>(_moi);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"forceflags") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'forceflags'");
          iarg_++;
          bool flags[3];
          if(strcmp("0",arg[iarg_++]) == 0) flags[0] = false;
          else                              flags[0] = true;
          if(strcmp("0",arg[iarg_++]) == 0) flags[1] = false;
          else                              flags[1] = true;
          if(strcmp("0",arg[iarg_++]) == 0) flags[2] = false;
          else                              flags[2] = true;
          fflag_.set(0,flags);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"torqueflags") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments for 'torqueflags'");
          iarg_++;
          bool flags[3];
          if(strcmp("0",arg[iarg_++]) == 0) flags[0] = false;
          else                              flags[0] = true;
          if(strcmp("0",arg[iarg_++]) == 0) flags[1] = false;
          else                              flags[1] = true;
          if(strcmp("0",arg[iarg_++]) == 0) flags[2] = false;
          else                              flags[2] = true;
          tflag_.set(0,flags);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"rot_flip") == 0) {
          if (narg < iarg_+3) error->fix_error(FLERR,this,"not enough arguments for 'rot_flip'");
          iarg_++;
          rot_flip_flag_ = true;
          if(strcmp(arg[iarg_++], "angle"))
            error->fix_error(FLERR,this,"expecting keyword 'angle'");
          else
            rot_flip_angle_ = force->numeric(FLERR,arg[iarg_++]);
          if(rot_flip_angle_ < 0. || rot_flip_angle_ > 90.)
            error->fix_error(FLERR,this,"0° < angle < 90° required");
          rot_flip_angle_ = rot_flip_angle_ * M_PI / 180.;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"suspension") == 0) {
          if (narg < iarg_+9) error->fix_error(FLERR,this,"not enough arguments for 'suspension'");
          iarg_++;
          suspension_flag_ = true;
          if(strcmp(arg[iarg_++], "k_t"))
            error->fix_error(FLERR,this,"expecting keyword 'k_t'");
          k_t_ = force->numeric(FLERR,arg[iarg_++]);
          if(k_t_ < 0.)
            error->fix_error(FLERR,this,"k_t >= 0 required");
          if(strcmp(arg[iarg_++], "c_t"))
            error->fix_error(FLERR,this,"expecting keyword 'c_t'");
          c_t_ = force->numeric(FLERR,arg[iarg_++]);
          if(c_t_ < 0.)
            error->fix_error(FLERR,this,"c_t >= 0 required");
          if(strcmp(arg[iarg_++], "k_r"))
            error->fix_error(FLERR,this,"expecting keyword 'k_r'");
          k_r_ = force->numeric(FLERR,arg[iarg_++]);
          if(k_r_ < 0.)
            error->fix_error(FLERR,this,"k_r >= 0 required");
          if(strcmp(arg[iarg_++], "c_r"))
            error->fix_error(FLERR,this,"expecting keyword 'c_r'");
          c_r_ = force->numeric(FLERR,arg[iarg_++]);
          if(c_r_ < 0.)
            error->fix_error(FLERR,this,"c_r >= 0 required");
          hasargs = true;
      } else if(strcmp(style,"mesh/surface/stress/6dof") == 0) {
          char *errmsg = new char[strlen(arg[iarg_])+50];
          sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
          error->fix_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }

    error_checks();

    // store original position and rotation state
    xcm_orig_.add(xcm_(0));
    quat_orig_.add(quat_(0));

    //NP inform mesh of upcoming movement
    mesh()->registerMove(false,true,true);
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStress6DOF::~FixMeshSurfaceStress6DOF()
{
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::post_create()
{
    FixMeshSurfaceStress::post_create();

    init_rotation_props();

    calc_displace();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::init_defaults()
{
    bool truevec[] ={true,true,true};
    fflag_.add(truevec);
    tflag_.add(truevec);
    double unitquat[4];
    quatUnitize4D(unitquat);
    quat_.add(unitquat);
    double zerovec[3];
    vectorZeroize3D(zerovec);
    vcm_.add(zerovec);
    angmom_.add(zerovec);
    omega_.add(zerovec);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::error_checks()
{
    /*NL*/ //fprintf(screen,"sizes xcm %d mass %d moi %d\n",xcm_.size(),mass_.size(),moi_.size());

    if(!xcm_.size())
        error->fix_error(FLERR,this,"please define 'com' for the 6dof mesh");
    if(!mass_.size())
        error->fix_error(FLERR,this,"please define 'mass' for the 6dof mesh");
    if(!moi_.size())
        error->fix_error(FLERR,this,"please define 'moi' for the 6dof mesh");
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::init_rotation_props()
{
  double **evectors;
  memory->create(evectors,3,3,"FixMeshSurfaceStress6DOF:evectors");
  /*NL*/ if(DEBUG_MESH_SURFACE_STRESS_6DOF) fprintf(screen,"Performing jacobi calc\n");

  int ierror = MathExtra::jacobi(moi_(0),inertia_(0),evectors);
  if (ierror) error->fix_error(FLERR,this,"Insufficient Jacobi rotations for rigid body");

  /*NL*/ if(DEBUG_MESH_SURFACE_STRESS_6DOF) fprintf(screen,"Jacobi calc finished\n");

  ex_space_(0)[0] = evectors[0][0];  ex_space_(0)[1] = evectors[1][0];  ex_space_(0)[2] = evectors[2][0];
  ey_space_(0)[0] = evectors[0][1];  ey_space_(0)[1] = evectors[1][1];  ey_space_(0)[2] = evectors[2][1];
  ez_space_(0)[0] = evectors[0][2];  ez_space_(0)[1] = evectors[1][2];  ez_space_(0)[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0
  /*NL*/ if(DEBUG_MESH_SURFACE_STRESS_6DOF) fprintf(screen,"Removing unnecessary intertia terms\n");

  double max;
  max = MathExtraLiggghts::max(inertia_(0)[0],inertia_(0)[1]);
  max = MathExtraLiggghts::max(max,inertia_(0)[2]);

  if (inertia_(0)[0] < EPSILON*max) inertia_(0)[0] = 0.0;
  if (inertia_(0)[1] < EPSILON*max) inertia_(0)[1] = 0.0;
  if (inertia_(0)[2] < EPSILON*max) inertia_(0)[2] = 0.0;

  //NP enforce 3 evectors as a right-handed coordinate system
  //NP flip 3rd evector if needed
  double ez0 = ex_space_(0)[1]*ey_space_(0)[2] -ex_space_(0)[2]*ey_space_(0)[1];
  double ez1 = ex_space_(0)[2]*ey_space_(0)[0] -ex_space_(0)[0]*ey_space_(0)[2];
  double ez2 = ex_space_(0)[0]*ey_space_(0)[1] -ex_space_(0)[1]*ey_space_(0)[0];

  if (ez0*ez_space_(0)[0] + ez1*ez_space_(0)[1] + ez2*ez_space_(0)[2] < 0.0)
  {
      ez_space_(0)[0] = -ez_space_(0)[0];
      ez_space_(0)[1] = -ez_space_(0)[1];
      ez_space_(0)[2] = -ez_space_(0)[2];
  }

  /*NL*/ if(DEBUG_MESH_SURFACE_STRESS_6DOF) fprintf(screen,"Creating initial quat\n");

  // create initial quaternion
  MathExtra::exyz_to_q(ex_space_(0),ey_space_(0),ez_space_(0),quat_(0));

  /*NL*/ if(DEBUG_MESH_SURFACE_STRESS_6DOF) fprintf(screen,"calculating dsplc\n");

  memory->destroy(evectors);

  /*NL*/ //print_results();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::calc_displace()
{
    //NP is called out of post_create() so mesh is parallel already
    //NP have to calculate displace for owned and ghost
    //NP displace values are then transported further on by container classes

    int nall = mesh()->size();
    int nnodes = mesh()->numNodes();
    double nodepos[3],tmp[3];

    double **dsplc;
    memory->create<double>(dsplc,nnodes,3,"6dof:dsplc");

    for(int i = 0; i < nall; i++)
    {
        for(int j = 0; j < nnodes; j++)
        {
            mesh()->node_slow(i,j,nodepos);
            vectorSubtract3D(nodepos,xcm_(0),tmp);
            dsplc[j][0] = vectorDot3D(tmp,ex_space_(0));
            dsplc[j][1] = vectorDot3D(tmp,ey_space_(0));
            dsplc[j][2] = vectorDot3D(tmp,ez_space_(0));
        }
        displace_.set(i,dsplc);
    }

    memory->destroy<double>(dsplc);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::init()
{
    FixMeshSurfaceStress::init();

    dtv_ = update->dt;
    dtf_ = 0.5 * update->dt * force->ftm2v;
    dtfm_ = dtf_ / mass_(0);

    dtq_ = 0.5 * update->dt;

    if (strcmp(update->integrate_style,"respa") == 0)
        error->fix_error(FLERR,this,"not respa-compatible");
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceStress6DOF::setmask()
{
    int mask = FixMeshSurfaceStress::setmask();
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::initial_integrate(int vflag)
{
    double dX[3],dx[3], qOld[4],dQ[4], dq[4];

    // update vcm by 1/2 step

    if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * f_total(0);
    if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * f_total(1);
    if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * f_total(2);

    // update xcm by full step

    dx[0] = dtv_ * vcm_(0)[0];
    dx[1] = dtv_ * vcm_(0)[1];
    dx[2] = dtv_ * vcm_(0)[2];
    vectorAdd3D(xcm_(0),dx,xcm_(0));
    vectorSubtract3D(xcm_(0),xcm_orig_(0),dX);

    // update angular momentum by 1/2 step

    if(tflag_(0)[0]) angmom_(0)[0] += dtf_ * torque_total(0);
    if(tflag_(0)[1]) angmom_(0)[1] += dtf_ * torque_total(1);
    if(tflag_(0)[2]) angmom_(0)[2] += dtf_ * torque_total(2);

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    vectorCopy4D(quat_(0),qOld);
    MathExtra::angmom_to_omega(angmom_(0),ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),omega_(0));
    MathExtra::richardson(quat_(0),angmom_(0),omega_(0),inertia_(0),dtq_);
    MathExtra::q_to_exyz(quat_(0),ex_space_(0),ey_space_(0),ez_space_(0));

    // calculate quaternion difference
    //NP this is numerically ok since original quaternion is stored
    //NP so nodes are updates accurately

    MathExtraLiggghts::quat_diff(quat_(0),qOld,         dq);
    MathExtraLiggghts::quat_diff(quat_(0),quat_orig_(0),dQ);

    // move and rotate mesh, set velocity

    /*NL*/ //printVec3D(screen,"dX",dX);
    /*NL*/ //printVec3D(screen,"dx",dx);
    /*NL*/ //printVec4D(screen,"dQ",dQ);
    /*NL*/ //printVec4D(screen,"dq",dq);

    mesh()->move(dX,dx);
    mesh()->rotate(dQ,dq,xcm_(0));
    set_vel();

    // update reference point to COM
    //NP would not be necessary b/c p_ref_ is moved, rotated automatically
    //NP do it anyway to avoid long-term divergence
    //NP which could happen b/c move, rotate is done increnmentally

    set_p_ref(xcm_(0));
    /*NL*/ //fprintf(screen,"p_ref %g %g %g\n",p_ref(0),p_ref(1),p_ref(2));
    /*NL*/ //printVec3D(screen,"xcm",xcm_(0));
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::final_integrate()
{
    // add extra forces

    add_gravity();

    if(suspension_flag_)
        add_suspension_force();

    if(rot_flip_flag_)
        rot_flip();

    //NP update forces
    FixMeshSurfaceStress::final_integrate();

    /*NL*/ //fprintf(screen,"Integration step "BIGINT_FORMAT", added force %f %f %f\n",update->ntimestep, f_total(0), f_total(1), f_total(2));

    // update vcm by 1/2 step

    if(fflag_(0)[0]) vcm_(0)[0] += dtfm_ * f_total(0);
    if(fflag_(0)[1]) vcm_(0)[1] += dtfm_ * f_total(1);
    if(fflag_(0)[2]) vcm_(0)[2] += dtfm_ * f_total(2);

    // update angular momentum by 1/2 step

    if(tflag_(0)[0]) angmom_(0)[0] += dtf_ * torque_total(0);
    if(tflag_(0)[1]) angmom_(0)[1] += dtf_ * torque_total(1);
    if(tflag_(0)[2]) angmom_(0)[2] += dtf_ * torque_total(2);

    MathExtra::angmom_to_omega(angmom_(0),ex_space_(0),ey_space_(0),ez_space_(0),inertia_(0),omega_(0));

    //NP update velocities of mesh points accordingly
    /*NL*/ //constrain_mesh();
    set_vel();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::add_gravity()
{
    double gravity[3],f_grav[3];

    vectorZeroize3D(gravity);
    if(modify->n_fixes_style_strict("gravity") > 0)
        static_cast<FixGravity*>(modify->find_fix_style_strict("gravity",0))->get_gravity(gravity);

    vectorScalarMult3D(gravity,mass_(0),f_grav);
    /*NL*/ //printVec3D(screen,"f_grav",f_grav);
    add_global_external_contribution(f_grav);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::rot_flip()
{
    double dangle;

    dangle = 2.*acos(quat_(0)[0]) - rot_flip_angle_;

    //NP flip velocity, angular momentum and angular velocity
    if(dangle > 0.)
    {
        vectorFlip3D(vcm_(0));
        vectorFlip3D(angmom_(0));
        vectorFlip3D(omega_(0));
    }
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::add_suspension_force()
{
    double dX[3], dQ[4], axis[3], angle, force[3], torque[3];
    double force_spring[3], force_damper[3], torque_spring[3], torque_damper[3];

    vectorSubtract3D(xcm_(0),xcm_orig_(0),dX);
    MathExtraLiggghts::quat_diff(quat_(0),quat_orig_(0),dQ);

    // force components

    vectorScalarMult3D(dX,-k_t_,force_spring);
    vectorScalarMult3D(vcm_(0),-c_t_,force_damper);

    // torque components
    //NP spring torque = k*angle*axis/abs(axis)
    //NP see http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/index.htm

    angle = 2.*acos(dQ[0]);
    MathExtraLiggghts::vec_from_quat(dQ,axis);
    vectorNormalize3D(axis);
    vectorScalarMult3D(axis,-k_r_*angle,torque_spring);
    vectorScalarMult3D(omega_(0),-c_r_,torque_damper);

    // add forces and contribute them

    vectorAdd3D(force_spring,force_damper,force);
    vectorAdd3D(torque_spring,torque_damper,torque);
    /*NL*/ //printVec3D(screen,"force",force);
    /*NL*/ //printVec3D(screen,"torque",torque);
    /*NL*/ //printVec4D(screen,"dQ",dQ);
    add_global_external_contribution(force,torque);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress6DOF::set_vel()
{
    int nall = mesh()->size();
    int nnodes = mesh()->numNodes();
    double dx,dy,dz;
    double **vNodes;

    memory->create<double>(vNodes,nnodes,3,"6dof:vNodes");

    /*NL*/// printVec3D(screen,"omega",omega_(0));
    /*NL*/// printVec3D(screen,"vcm",vcm_(0));
    /*NL*/// printVec3D(screen,"ex_space_",ex_space_(0));
    /*NL*/// printVec3D(screen,"ey_space_",ey_space_(0));
    /*NL*/// printVec3D(screen,"ez_space_",ez_space_(0));

    for(int i = 0; i < nall; i++)
    {
        for(int j = 0; j < nnodes; j++)
        {
            dx = ex_space_(0)[0]*displace_(i)[j][0] + ey_space_(0)[0]*displace_(i)[j][1] +  ez_space_(0)[0]*displace_(i)[j][2];
            dy = ex_space_(0)[1]*displace_(i)[j][0] + ey_space_(0)[1]*displace_(i)[j][1] +  ez_space_(0)[1]*displace_(i)[j][2];
            dz = ex_space_(0)[2]*displace_(i)[j][0] + ey_space_(0)[2]*displace_(i)[j][1] +  ez_space_(0)[2]*displace_(i)[j][2];
            vNodes[j][0] = omega_(0)[1]*dz - omega_(0)[2]*dy + vcm_(0)[0];
            vNodes[j][1] = omega_(0)[2]*dx - omega_(0)[0]*dz + vcm_(0)[1];
            vNodes[j][2] = omega_(0)[0]*dy - omega_(0)[1]*dx + vcm_(0)[2];
        }
        v_.set(i,vNodes);
    }

    memory->destroy<double>(vNodes);
}
