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

#include "fix_mesh_surface_stress_contact.h"
#include <stdio.h>
#include "string.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"
#include "pair_gran.h"
#include "mech_param_gran.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "fix_gravity.h"

#define EPSILON 0.001

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressContact::FixMeshSurfaceStressContact(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurfaceStress(lmp, narg, arg),
  fix_wallcontacttime_(0),
  T_(0.),
  step_ave_start_(-1),
  area_correction_(false),
  deltan_ratio_(0)
{

  bool hasargs = true;
  while(iarg_ < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg_],"area_correction") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'area_correction'");
      if(strcmp(arg[iarg_+1],"yes") == 0)
        area_correction_ = true;
      else if(strcmp(arg[iarg_+1],"no") == 0)
        area_correction_ = false;
      else error->fix_error(FLERR,this,"");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(arg[iarg_],"start_step_ave") == 0) {
      if (iarg_+2 > narg) error->fix_error(FLERR,this,"not enough arguments for keyword 'start_step_ave'");
      step_ave_start_ = atoi(arg[iarg_+1]);
      if(step_ave_start_ < 0)
        error->fix_error(FLERR,this,"'start_step_ave' >= 0 required");
      iarg_ += 2;
      hasargs = true;
    } else if(strcmp(style,"mesh/surface/stress/contact") == 0)
        error->fix_error(FLERR,this,"unknown keyword");
  }
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressContact::~FixMeshSurfaceStressContact()
{

}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::post_create()
{
    FixMeshSurfaceStress::post_create();

    mesh()->prop().addElementProperty<ScalarContainer<double> >("contact_area",    "comm_exchange_borders","frame_invariant","restart_yes");
    mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area")->setAll(0.);
    mesh()->prop().addElementProperty<ScalarContainer<double> >("contact_area_abs","comm_exchange_borders","frame_invariant","restart_yes");
    mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_abs")->setAll(0.);

    mesh()->prop().addElementProperty<ScalarContainer<double> >("contact_area_step",    "comm_reverse","frame_invariant","restart_no");
    mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_step")->setAll(0.);
    mesh()->prop().addElementProperty<ScalarContainer<double> >("contact_area_step_abs","comm_reverse","frame_invariant","restart_no");
    mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_step_abs")->setAll(0.);

    char *wallctime_name = new char[strlen(style)+1+12];
    strcpy(wallctime_name,"contacttime_");
    strcat(wallctime_name,id);
    char **fixarg = new char*[11];
    fixarg[0] = wallctime_name;
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "property/atom";
    fixarg[3] = wallctime_name;
    fixarg[4] = (char *) "scalar";
    fixarg[5] = (char *) "yes";   // restart
    fixarg[6] = (char *) "no";    // communicate ghost
    fixarg[7] = (char *) "no";    // communicate rev
    fixarg[8] = (char *) "0.";
    modify->add_fix(9,fixarg);
    fix_wallcontacttime_ =
         static_cast<FixPropertyAtom*>(modify->find_fix_property(wallctime_name,"property/atom","scalar",0,0,style));
    delete []fixarg;
    delete []wallctime_name;
}

/* ---------------------------------------------------------------------- */

int FixMeshSurfaceStressContact::setmask()
{
    int mask = FixMeshSurfaceStress::setmask();
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::init()
{
    if(!trackStress())
        error->fix_error(FLERR,this,"need 'stress' = 'on'");

    FixMeshSurfaceStress::init();

    init_area_correction();

    contact_area_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area");
    contact_area_abs_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_abs");

    contact_area_step_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_step");
    contact_area_step_abs_ = mesh()->prop().getElementProperty<ScalarContainer<double> >("contact_area_step_abs");

    if(!contact_area_ || !contact_area_step_ || !contact_area_abs_ || !contact_area_step_abs_)
        error->one(FLERR,"internal error");

    if(force->cg() > 1.)
        error->fix_error(FLERR,this,"does not support coarse-graining");
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::init_area_correction()
{
    const double *Y, *nu, *Y_orig;
    double expo, Yeff_ij, Yeff_orig_ij, ratio;

    if(area_correction_)
    {
        PairGran *pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
        if(!pair_gran)
            error->fix_error(FLERR,this,"'area_correction' requires using a granular pair style");
        int max_type = pair_gran->mpg->max_type();

        if(force->pair_match("gran/hooke",0)) expo = 1.;
        else if(force->pair_match("gran/hertz",0)) expo = 2./3.;
        else error->fix_error(FLERR,this,"area correction could not identify the granular pair style you are using, supported are hooke and hertz types");

        Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
        nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
        Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_values();

        // allocate a new array within youngsModulusOriginal
        static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->new_array(max_type,max_type);

        // feed deltan_ratio into this array
        for(int i = 1; i < max_type+1; i++)
        {
          for(int j = 1; j < max_type+1; j++)
          {
            Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
            Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
            ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
            /*NL*/ //fprintf(screen,"ratio for type pair %d/%d is %f, Yeff_ij %f, Yeff_orig_ij %f, expo %f\n",i,j,ratio,Yeff_ij,Yeff_orig_ij,expo);
            static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->array_modify(i-1,j-1,ratio);
          }
        }

        // get reference to deltan_ratio
        deltan_ratio_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0,style))->get_array_modified();
    }
}

/* ----------------------------------------------------------------------
   re-set per-time-step contact area
------------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::pre_force(int vflag)
{
    FixMeshSurfaceStress::pre_force(vflag);

    int nTri = mesh()->size();

    for(int i = 0; i < nTri; i++)
    {
        contactAreaStep(i) = 0.;
        contactAreaStepAbs(i) = 0.;
    }
}

/* ----------------------------------------------------------------------
   called during wall force calc
------------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::add_particle_contribution(int ip,double *frc,
                                double *delta,int iTri,double *v_wall)
{
    if(!(atom->mask[ip] & groupbit)) return;

    FixMeshSurfaceStress::add_particle_contribution(ip,frc,delta,iTri,v_wall);

    double delta_n, r;
    double dt = update->dt;
    double rsq = vectorMag3DSquared(delta);
    double radius = atom->radius[ip];

    //NP adjust overlap that may be superficially large due to softening
    if(area_correction_)
    {
        delta_n = radius - sqrt(r);
        delta_n *= deltan_ratio_[atom->type[ip]-1][atomTypeWall()-1];
        r = radius - delta_n;
        rsq = r*r;
    }

    double Acont = (radius*radius-rsq)*M_PI; //contact area sphere-wall

    contactAreaStep(iTri) += Acont / triMesh()->areaElem(iTri);
    contactAreaStepAbs(iTri) += Acont;
    fix_wallcontacttime_->vector_atom[ip] += dt;
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::final_integrate()
{
    //NP wear_step is communicated by reverse comm
    FixMeshSurfaceStress::final_integrate();

    calc_contact_time_average();
}

/* ---------------------------------------------------------------------- */

void FixMeshSurfaceStressContact::calc_contact_time_average()
{
    if(step_ave_start_ > update->ntimestep)
        return;

    int nTri = mesh()->size();
    double dt = update->dt;

    for(int i = 0; i < nTri; i++)
    {
        contactArea(i) = T_ * contactArea(i) + dt * contactAreaStep(i);
        contactAreaAbs(i) = T_ * contactAreaAbs(i) + dt * contactAreaStepAbs(i);
        /*NL*/// if(contactArea(i) > 0.) fprintf(screen,"contactArea(i) %f\n",contactArea(i) );
        contactArea(i)    /= (T_+dt);
        contactAreaAbs(i) /= (T_+dt);
    }
    T_ += dt;
}

