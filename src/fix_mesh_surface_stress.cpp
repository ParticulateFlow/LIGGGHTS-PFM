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
   Philippe Seil (JKU Linz)
   Evan Smuts (U Cape Town, surface velocity rotation)
------------------------------------------------------------------------- */

#include "fix_mesh_surface_stress.h"
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_gravity.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define EPSILON 0.0001

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStress::FixMeshSurfaceStress(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurface(lmp, narg, arg),
  //NP reference point is moved, rotated etc automatically
  //NP cannot do this in post_create since restart data would not have any influence
  p_ref_(*mesh()->prop().addGlobalProperty<VectorContainer<double,3> >("p_ref","comm_none","frame_general","restart_yes")),
  f_(0),
  sigma_n_(0),
  sigma_t_(0),
  wear_flag_(0),
  k_finnie_(0),
  wear_(0),
  wear_step_(0)
{
    vectorZeroize3D(f_total_);
    vectorZeroize3D(torque_total_);

    double zerovec[3] = {0.,0.,0.};
    mesh()->prop().setGlobalProperty<VectorContainer<double,3> >("p_ref",zerovec);

    // override default from base
    stress_flag_ = true;

    vector_flag = 1;
    size_vector = 9;
    global_freq = 1;
    extvector = 1;

    // parse further args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
      hasargs = false;
      if (strcmp(arg[iarg_],"reference_point") == 0) {
          if (narg < iarg_+4) error->fix_error(FLERR,this,"not enough arguments");
          if(manipulated())
            error->warning(FLERR,"Mesh for fix mesh/surface/stress has been scaled, moved, or rotated.\n"
                             "Please note that values for 'reference_point' refer to the scaled, moved, or rotated configuration");
          iarg_++;
          double _p_ref[3];
          _p_ref[0] = force->numeric(FLERR,arg[iarg_++]);
          _p_ref[1] = force->numeric(FLERR,arg[iarg_++]);
          _p_ref[2] = force->numeric(FLERR,arg[iarg_++]);
          mesh()->prop().setGlobalProperty<VectorContainer<double,3> >("p_ref",_p_ref);
          hasargs = true;
      } else if(strcmp(arg[iarg_],"stress") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          if(strcmp(arg[iarg_],"on") == 0) stress_flag_ = true;
          else if(strcmp(arg[iarg_],"off") == 0) stress_flag_ = false;
          else error->fix_error(FLERR,this,"expecting 'on' or 'off' as stress argument");
          iarg_++;
          hasargs = true;
      } else if(strcmp(arg[iarg_],"wear") == 0) {
          if (narg < iarg_+2) error->fix_error(FLERR,this,"not enough arguments");
          iarg_++;
          if(strcmp(arg[iarg_],"finnie") == 0) wear_flag_ = 1;
          else if(strcmp(arg[iarg_],"off") == 0) wear_flag_ = 0;
          else error->fix_error(FLERR,this,"expecting 'finnie' or 'off' as wear argument");
          iarg_++;
          hasargs = true;
      } else if(strcmp(style,"mesh/surface/stress") == 0) {
          char *errmsg = new char[strlen(arg[iarg_])+50];
          sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg_]);
          error->fix_error(FLERR,this,errmsg);
          delete []errmsg;
      }
    }
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStress::~FixMeshSurfaceStress()
{

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::post_create_pre_restart()
{
    //NP called after constructor and before restart
    //NP i.e. no restart: mesh elements present at this point
    //NP restart: no mesh elements present because created in restart()
    //Np -->register properties and set values for non-restart properties here

    if(stress_flag_)
        regStress();

    if(wear_flag_)
        regWear();
}

/* ---------------------------------------------------------------------- */

//NP in case of no restart, mem for properties is allocated directly
//NP by addElementProperty because mesh elements are already present

//NP restart allocates mem for all element properties (unititialized)
//NP (callback via addElement), but are unitilized, so have to init
//NP them in post_create

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::post_create()
{
    //NP called after constructor and restart
    //NP i.e. no restart: mesh elements present at this point
    //NP restart: no mesh elements present because created in restart()
    //Np --> set values for no-restart properties here

    if(stress_flag_)
        zeroizeStress();

    if(wear_flag_)
        zeroizeWear();

    FixMeshSurface::post_create();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::regStress()
{
    mesh()->prop().addElementProperty<VectorContainer<double,3> >("f","comm_reverse","frame_invariant","restart_no");
    mesh()->prop().addElementProperty<ScalarContainer<double> >("sigma_n","comm_none","frame_invariant","restart_no");
    mesh()->prop().addElementProperty<ScalarContainer<double> >("sigma_t","comm_none","frame_invariant","restart_no");
}
/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::zeroizeStress()
{
    mesh()->prop().getElementProperty<VectorContainer<double,3> >("f")->setAll(0.);
    mesh()->prop().getElementProperty<ScalarContainer<double> >("sigma_n")->setAll(0.);
    mesh()->prop().getElementProperty<ScalarContainer<double> >("sigma_t")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::regWear()
{
    mesh()->prop().addElementProperty<ScalarContainer<double> >("wear","comm_exchange_borders","frame_invariant","restart_yes");
    mesh()->prop().getElementProperty<ScalarContainer<double> >("wear")->setAll(0.);
    mesh()->prop().addElementProperty<ScalarContainer<double> >("wear_step","comm_reverse","frame_invariant","restart_no");
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::zeroizeWear()
{
    mesh()->prop().getElementProperty<ScalarContainer<double> >("wear_step")->setAll(0.);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::init()
{
    if(stress_flag_)
    {
        f_ = mesh()->prop().getElementProperty<VectorContainer<double,3> >("f");
        sigma_n_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("sigma_n");
        sigma_t_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("sigma_t");
        if(!f_ || !sigma_n_ || !sigma_t_)
            error->one(FLERR,"Internal error");
    }

    if(wear_flag_)
    {
        k_finnie_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("k_finnie","property/global","peratomtypepair",atom->ntypes,atom->ntypes,style))->get_array();
        wear_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("wear");
        wear_step_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("wear_step");
        if(!wear_ || ! wear_step_)
            error->one(FLERR,"Internal error");
    }
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceStress::setmask()
{
    int mask = FixMeshSurface::setmask();
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::setup_pre_force(int vflag)
{
    FixMeshSurface::setup_pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::pre_force(int vflag)
{
    FixMeshSurface::pre_force(vflag);

    if(trackStress())
    {
        vectorZeroize3D(f_total_);
        vectorZeroize3D(torque_total_);
    }

    /*NL*/ //if (screen) fprintf(screen,"force on tri0: %f %f %f\n",f(0)[0],f(0)[2],f(0)[2]);
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStress::final_integrate()
{
    //NP wear_step is communicated by reverse comm
    FixMeshSurface::final_integrate();

    calc_total_force();
}

/* ----------------------------------------------------------------------
   called during wall force calc
------------------------------------------------------------------------- */

void FixMeshSurfaceStress::add_particle_contribution(int ip,double *frc,
                                double *delta,int iTri,double *v_wall)
{
    double E,c[3],v_rel[3],v_rel_mag,cos_gamma,sin_gamma,sin_2gamma;
    double contactPoint[3]={},surfNorm[3], tmp[3], tmp2[3];

    // do not include if not in fix group
    if(!(atom->mask[ip] & groupbit)) return;

    double *x = atom->x[ip];
    double *v = atom->v[ip];

    vectorNegate3D(frc);

    vectorAdd3D(x,delta,contactPoint);

    if(trackStress())
    {
        /*NL*/ //if(screen && strcmp(id,"servo")==0) printVec3D(screen,"added force",frc);
        /*NL*/ //if(screen) fprintf(screen,"step " BIGINT_FORMAT ", added force %f %f %f\n",update->ntimestep,frc[0],frc[1],frc[2]);

        // add contribution to triangle force
        vectorAdd3D(f(iTri),frc,f(iTri));

        // add contribution to total body force and torque
        vectorAdd3D(f_total_,frc,f_total_);
        vectorSubtract3D(contactPoint,p_ref_(0),tmp);
        /*NL*/ //if (screen) fprintf(screen,"p_ref_ %f %f %f\n",p_ref_(0)[0],p_ref_(0)[1],p_ref_(0)[2]);
        /*NL*/ //if (screen) fprintf(screen,"tmp %f %f %f\n",tmp[0],tmp[1],tmp[2]);
        vectorCross3D(tmp,frc,tmp2); // tmp2 is torque contrib
        vectorAdd3D(torque_total_,tmp2,torque_total_);
    }

    // add wear if applicable
    if(trackWear())
    {
        //NP c pointing from particle center to contact point
        vectorSubtract3D(contactPoint,x,c);

        // calculate relative velocity
        vectorSubtract3D(v,v_wall,v_rel);

        //NP if v.c < 0 return
        //NP this ensures only the approach is taken into acctount
        //NP not the rebound
        if(vectorDot3D(c,v_rel) < 0.) return;
        v_rel_mag = vectorMag3D(v_rel);

        // get surface normal
        //NP surface normal is normalized
        triMesh()->surfaceNorm(iTri,surfNorm);

        // return if no relative velocity
        if(0.0000001 > v_rel_mag)
            return;

        //NP exchange sin and cos because
        //NP vectorDot/cross will give angle to surface normal
        //NP not to wall plane
        //NP cos(gamma) = sin(90-gamma)
        //NP sin(gamma) = cos(90-gamma)
        //NP need abs so that can cover surfNorm pointing in any direction

        sin_gamma = MathExtraLiggghts::abs(vectorDot3D(v_rel,surfNorm)) / (v_rel_mag);
        cos_gamma = MathExtraLiggghts::abs(vectorCrossMag3D(v_rel,surfNorm)) / (v_rel_mag);

        //NP both angles from 0..90 degree
        //NP both sin and cos will be > 0
        if(cos_gamma > 1.) cos_gamma = 1.;
        if(sin_gamma > 1.) sin_gamma = 1.;

        /*NL*/ //if (screen) fprintf(screen,"sin_gamma %e cos_gamma %e, stl normal mag %f ",sin_gamma,cos_gamma,vectorMag3D(surfNorm));

        //NP http://www1.ansys.com/customer/content/documentation/120/cfx/xthry.pdf
        //NP Eqn 5.70

        //NP cos_gamma cannot be < 0
        if(cos_gamma < EPSILON || 3.*sin_gamma > cos_gamma)
        {
            E = 0.33333 * cos_gamma * cos_gamma;
            /*NL*///if (screen) fprintf(screen," wear1 %e\n",E);
        }
        else
        {
            sin_2gamma = 2. * sin_gamma * cos_gamma;
            E = sin_2gamma - 3. * sin_gamma * sin_gamma;
            /*NL*///if (screen) fprintf(screen," wear2 %e\n",E);
        }
        E *= 2.*k_finnie_[atomTypeWall()-1][atom->type[ip]-1] * v_rel_mag * vectorMag3D(frc);
        //NP if (screen) fprintf(screen," k_finnie %f, vmag%f, frcmag %f ,wear %1.15f\n",k_finnie[atom_type_wall-1][atom->type[ip]-1],vmag,vectorMag3D(frc),E);
        //NP error->all("wear");

        wear_step(iTri) += E*update->dt / triMesh()->areaElem(iTri);
    }
}

/* ----------------------------------------------------------------------
   add external force (such as gravity)
   called by all procs, only proc0 adds
   has to be called before final_integrate()
------------------------------------------------------------------------- */

void FixMeshSurfaceStress::add_global_external_contribution(double *frc)
{
    if(0 == comm->me)
        vectorAdd3D(f_total_,frc,f_total_);
}

void FixMeshSurfaceStress::add_global_external_contribution(double *frc, double *trq)
{
    if(0 == comm->me)
    {
        vectorAdd3D(f_total_,frc,f_total_);
        vectorAdd3D(torque_total_,trq,torque_total_);
    }
}

/* ----------------------------------------------------------------------
   allreduce total force on tri
------------------------------------------------------------------------- */

void FixMeshSurfaceStress::calc_total_force()
{
    double surfNorm[3], invSurfArea, temp[3];
    int nTri = mesh()->size();

    // add wear from this step to total wear
    if(trackWear())
    {
        for(int i = 0; i < nTri; i++)
        {
            wear(i) += wear_step(i);
            wear_step(i) = 0.;
        }
    }

    // calculate normal and shear stress
    if(trackStress())
    {
        // total force and torque on mesh

        MPI_Sum_Vector(f_total_,3,world);
        MPI_Sum_Vector(torque_total_,3,world);
        /*NL*/ //if (screen) printVec3D(screen,"f_total_",f_total_);

        for(int i = 0; i < nTri; i++)
        {
            // get element surface norm and area
            triMesh()->surfaceNorm(i,surfNorm);
            invSurfArea = 1./triMesh()->areaElem(i);

            // calculate normal force
            sigma_n(i) = vectorDot3D(f(i),surfNorm);
            /*NL*/ //if (screen) printVec3D(screen,"force tri",f(i));

            // calculate tangential force
            vectorScalarMult3D(surfNorm,sigma_n(i),temp);
            vectorSubtract3D(f(i),temp,temp);
            sigma_t(i) = vectorMag3D(temp);

            // make both positive
            // necessary since orientation of surfNorm not known
            sigma_n(i) = MathExtraLiggghts::abs(sigma_n(i));
            sigma_t(i) = MathExtraLiggghts::abs(sigma_t(i));

            // divide by area so have stress
            sigma_n(i) *= invSurfArea;
            sigma_t(i) *= invSurfArea;
            /*NL*/ //if (screen) fprintf(screen,"sigma_n tri %f\n",sigma_n(i));
            /*NL*/ //if (screen) fprintf(screen,"sigma_t tri %f\n",sigma_t(i));
        }
    }
}

/* ----------------------------------------------------------------------
   return total force or torque component on body
------------------------------------------------------------------------- */

double FixMeshSurfaceStress::compute_vector(int n)
{
  if(n < 3) return f_total_[n];
  else if(n < 6)     return torque_total_[n-3];
  else return p_ref_(0)[n-6];
}
