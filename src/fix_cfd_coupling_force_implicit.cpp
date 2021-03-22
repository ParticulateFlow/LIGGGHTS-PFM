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

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_force_implicit.h"
#include "fix_property_atom.h"
#include "fix_nve_sphere.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::FixCfdCouplingForceImplicit(LAMMPS *lmp, int narg, char **arg) :
    FixCfdCouplingForce(lmp,narg,arg),
    useCN_(false),
    CNalpha_(0.0),
    useAM_(false), // superquadric
    CAddRhoFluid_(0.0), // superquadric
    onePlusCAddRhoFluid_(1.0), // superquadric
    implicitIntegration_(false),
    fix_Ksl_(0),
    fix_uf_(0),
    fix_KslRotation_(0), // superquadric
    fix_ex_(0), // superquadric
    fix_KslExtra_(0) // superquadric
{
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"CrankNicolson") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CrankNicholson'");
            iarg++;
            useCN_ = true;
            CNalpha_ = atof(arg[iarg]);
            if(CNalpha_<0 || CNalpha_>1)
                error->fix_error(FLERR,this,"incorrect choice for 'CrankNicholson': setting CNalpha_<0 or CNalpha_>1 is not appropriate");
            if (screen) fprintf(screen,"cfd_coupling_foce_implicit will use Crank-Nicholson scheme with %f\n", CNalpha_);
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"implicit_integration") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'implicit_integration'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                implicitIntegration_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                implicitIntegration_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'implicit_integration'");
            iarg++;
            hasargs = true;
        }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        else if (strcmp(arg[iarg],"CAddRhoFluid") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'CAddRhoFluid'");
            iarg++;
            useAM_ = true;
            CAddRhoFluid_        = atof(arg[iarg]);
            onePlusCAddRhoFluid_ = 1.0 + CAddRhoFluid_;
            if (screen) fprintf(screen,"cfd_coupling_force_implicit will consider added mass with CAddRhoFluid = %f\n", CAddRhoFluid_);
            iarg++;
            hasargs = true;
        }
#endif
    }

  nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingForceImplicit::~FixCfdCouplingForceImplicit()
{

}

/* ---------------------------------------------------------------------- */
int FixCfdCouplingForceImplicit::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_create()
{
    // do mother class init w/o dragforce
    FixCfdCouplingForce::post_create();

    // register Ksl
    if(!fix_Ksl_)
    {
        const char* fixarg[9];
        fixarg[0]="Ksl";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="Ksl";
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_Ksl_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    // register uf
    if(!fix_uf_)
    {
        const char* fixarg[11];
        fixarg[0]="uf";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="uf";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_uf_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(!fix_KslRotation_)
    {
      const char* fixarg[11];
      fixarg[0]="KslRotation";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="KslRotation";
      fixarg[4]="vector"; // 1 vector per particle to be registered
      fixarg[5]="yes";    // restart
      fixarg[6]="no";     // communicate ghost
      fixarg[7]="no";     // communicate rev
      fixarg[8]= "0.";
      fixarg[9]= "0.";
      fixarg[10]="0.";
      fix_KslRotation_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
    if(!fix_ex_)
    {
       const char* fixarg[11];
       fixarg[0]="ex";
       fixarg[1]="all";
       fixarg[2]="property/atom";
       fixarg[3]="ex";
       fixarg[4]="vector"; // 1 vector per particle to be registered
       fixarg[5]="yes";    // restart
       fixarg[6]="no";     // communicate ghost
       fixarg[7]="no";     // communicate rev
       fixarg[8]= "0.";
       fixarg[9]= "0.";
       fixarg[10]="0.";
       fix_ex_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
    if(!fix_KslExtra_)
    {
      const char* fixarg[11];
      fixarg[0]="KslExtra";
      fixarg[1]="all";
      fixarg[2]="property/atom";
      fixarg[3]="KslExtra";
      fixarg[4]="vector"; // 1 vector per particle to be registered
      fixarg[5]="yes";    // restart
      fixarg[6]="no";     // communicate ghost
      fixarg[7]="no";     // communicate rev
      fixarg[8]= "0.";
      fixarg[9]= "0.";
      fixarg[10]="0.";
      fix_KslExtra_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }
#endif
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_Ksl_) modify->delete_fix("Ksl");
    if(unfixflag && fix_uf_) modify->delete_fix("uf");
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(unfixflag && fix_KslRotation_) modify->delete_fix("KslRotation");
    if(unfixflag && fix_KslExtra_) modify->delete_fix("KslExtra");
    if(unfixflag && fix_ex_) modify->delete_fix("ex");
#endif
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::init()
{
    FixCfdCouplingForce::init();

    // values to come from OF
    fix_coupling_->add_pull_property("Ksl","scalar-atom");
    fix_coupling_->add_pull_property("uf","vector-atom");
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    fix_coupling_->add_pull_property("KslRotation","vector-atom");
    fix_coupling_->add_pull_property("KslExtra","vector-atom");
    fix_coupling_->add_pull_property("ex","vector-atom");
#endif

    deltaT_ = 0.5 * update->dt * force->ftm2v;

    if (implicitIntegration_)
    {
        fix_nve_sphere_ = static_cast<FixNVESphere*>(modify->find_fix_style("nve/sphere",0));
        if (!fix_nve_sphere_)
        {
            error->fix_error(FLERR,this,"Could not find fix ID 'nve/sphere' required for fully implicit treatment of the drag force.");
        }
        if (!fix_nve_sphere_->implicitIntegration())
        {
            error->fix_error(FLERR,this,"Fix 'nve/sphere' is not in mode for implicit integration.");
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::post_force(int)
{
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;
  double **dragforce = fix_dragforce_->array_atom;
  double frc[3];

  vectorZeroize3D(dragforce_total);
  vectorZeroize3D(hdtorque_total); // superquadric

  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        // calc force
        if(!useCN_)  //calculate drag force and add if not using Crank-Nicolson
        {
            if (implicitIntegration_)
            {
                frc[0] = uf[i][0];
                frc[1] = uf[i][1];
                frc[2] = uf[i][2];
            }
            else
            {
                vectorSubtract3D(uf[i],v[i],frc);
            }
            vectorScalarMult3D(frc,Ksl[i]);

            vectorAdd3D(f[i],frc,f[i]);
            vectorAdd3D(dragforce_total,frc,dragforce_total);
        }

        // add other force
        vectorAdd3D(f[i],dragforce[i],f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingForceImplicit::end_of_step()
{

  if(!useCN_) return; //return if CN not used

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ksl = fix_Ksl_->vector_atom;
  double **uf = fix_uf_->array_atom;
  //double **dragforce = fix_dragforce_->array_atom;
  double KslMDeltaT, deltaU;
  double vN32[3];
  double frc[3];

  vectorZeroize3D(dragforce_total);


  // add dragforce to force vector
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        if (rmass)  KslMDeltaT = Ksl[i]/(rmass[i]*onePlusCAddRhoFluid_)*deltaT_; // superquadric
        else        KslMDeltaT = Ksl[i]/(mass[type[i]]*onePlusCAddRhoFluid_)*deltaT_; // superquadric

        for(int dirI=0;dirI<3;dirI++)
        {
            //calculate new velocity
            vN32[dirI] = (  v[i][dirI]
                          + KslMDeltaT
                            *(   uf[i][dirI]
                              - (1.0-CNalpha_)*v[i][dirI]
                             )
                         )
                         /
                         (1.0+KslMDeltaT*CNalpha_);

            //calculate velocity difference and force
            deltaU    =  uf[i][dirI]
                           - (
                                (1.0-CNalpha_)*v[i][dirI]
                               +     CNalpha_ *vN32[dirI]
                             );
           frc[dirI] = Ksl[i] * deltaU;  //force required for the next time step

           //update the particle velocity
           v[i][dirI] = vN32[dirI];  //update velocity for a half step!
           //v[i][dirI] += KslMDeltaT/2.0 * deltaU;  //update velocity for a half step!
        }

         // add force
        vectorAdd3D(f[i],frc,f[i]);

        // add up forces for post-proc
        vectorAdd3D(dragforce_total,frc,dragforce_total);
    }
  }
}
