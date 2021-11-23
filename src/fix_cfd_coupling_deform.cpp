/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particulate Flow Modelling
   Copyright 2015-     JKU Linz

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
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include <math.h>
#include "vector_liggghts.h"
#include "math_const.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_deform.h"
#include "fix_property_atom.h"
#include "fix_property_atom_polydispparcel.h"
#include "fix_property_global.h"
#include "pair_gran.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define SMALL 1e-6

/* ---------------------------------------------------------------------- */

FixCfdCouplingDeform::FixCfdCouplingDeform(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    fix_coupling_(0),
    fix_partdeformations_(0),
    fix_effvolfactors_(0),
    verbose_(false),
    compress_flag_(1),
    mass_removed_(0.0),
    rmin_(0.0),
    fmax_(1.5), // 1/alpha_max \approx 1.5
    default_effvolfactors_(NULL),
    monitor_heat_(false),
    heat_removed_(0.),
    fix_temp_(NULL),
    fix_capacity_(NULL)
{

  int iarg = 3;

  while (iarg < narg) {
    if(strcmp(arg[iarg],"compress") == 0)
    {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"Illegal compress option");
      if(strcmp(arg[iarg+1],"yes") == 0)
        compress_flag_ = 1;
      else if(strcmp(arg[iarg+1],"no") == 0)
        compress_flag_ = 0;
      else
        error->fix_error(FLERR,this,"Illegal compress option, expecing 'yes' or 'no'");
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"verbose") == 0)
    {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'verbose'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        verbose_ = true;
      else if(strcmp(arg[iarg+1],"no"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'verbose'");
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"rmin") == 0)
    {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'rmin'");
      rmin_ = atof(arg[iarg+1]);
      if(rmin_ < 0.0)
        error->fix_error(FLERR,this,"rmin cannot be negative");
      iarg += 2;
    }
    else if(strcmp(arg[iarg],"monitor_heat") == 0) {
      if(narg < iarg+2)
        error->fix_error(FLERR,this,"not enough arguments for 'monitor_heat'");
      if(strcmp(arg[iarg+1],"yes") == 0)
        monitor_heat_ = true;
      else if(strcmp(arg[iarg+1],"no"))
        error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'monitor_heat_'");
      iarg += 2;
    }
    else error->fix_error(FLERR,this,"unknown keyword");
  }

  vector_flag = 1;
  if (monitor_heat_)
  {
    size_vector = 2;
  }
  else
  {
    size_vector = 1;
  }
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingDeform::~FixCfdCouplingDeform()
{
    delete [] default_effvolfactors_;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDeform::post_create()
{
//  register convective flux
    if(!fix_partdeformations_)
    {
        const char* fixarg[9];
        fixarg[0]="partDeformations";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="partDeformations";
        fixarg[4]="scalar";
        fixarg[5]="yes";    //NP restart
        fixarg[6]="no";    //NP communicate ghost
        fixarg[7]="no";    //NP communicate rev
        fixarg[8]="0.";
        fix_partdeformations_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    fix_effvolfactors_ = static_cast<FixPropertyAtomPolydispParcel*>(modify->find_fix_property("effvolfactor","property/atom/polydispparcel","scalar",0,0,style,false));
    if (!fix_effvolfactors_)
    {
        error->fix_error(FLERR,this,"Fix couple/cfd/deform needs a fix for type eff. vol. factor");
    }

    int num_defaultvalues = fix_effvolfactors_->num_defaultvalues();
    default_effvolfactors_ = new double[num_defaultvalues];
    for (int i=0;i<num_defaultvalues;i++)
    {
        default_effvolfactors_[i] = fix_effvolfactors_->defaultvalue(i);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDeform::pre_delete(bool unfixflag)
{
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingDeform::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDeform::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/deform needs a fix of type couple/cfd");

    // values to come from OF
    fix_coupling_->add_pull_property("partDeformations","scalar-atom");

    if (monitor_heat_)
    {
        fix_temp_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp", "property/atom", "scalar", 0, 0, style));

        PairGran* pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
        int max_type = pair_gran->get_properties()->max_type();
        fix_capacity_ = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalCapacity","property/global","peratomtype",max_type,0,style));
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingDeform::post_force(int)
{
    // deformation = 0 --> unchanged
    // deformation = 1 --> maximum deformation before removal
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int *type = atom->type;
    double *density = atom->density;
    double *radius = atom->radius;
    double *rmass = atom->rmass;

    double *partdeformations = fix_partdeformations_->vector_atom;
    double *effvolfactors_ = fix_effvolfactors_->vector_atom;
    double effvolfactor = 1.0;
    double neweffvolfactor = 1.0;
    double deformation = 0.0;

    int nremoved_this_me = 0;
    int nremoved_this = 0;
    double mass_removed_this_me = 0.0;
    double mass_removed_this = 0.0;

    double newradius = 0.0;

    double f0 = 0.0;

    double *T = NULL;
    if(monitor_heat_)
    {
        T = fix_temp_->vector_atom;
    }
    double heat_removed_this_me = 0., heat_removed_this = 0.;

    for (int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            deformation = partdeformations[i];
            if (deformation < 1.0 - SMALL && deformation > SMALL)
            {
                effvolfactor = effvolfactors_[i];
                f0 = default_effvolfactors_[type[i]];
                neweffvolfactor = f0 + deformation * (fmax_ - f0);
                // update properties only if new eff vol factor larger than old one (by e.g. 1%) 
                if (neweffvolfactor <= 1.01*effvolfactor) continue;

                // deform particle such that it keeps its volume
                newradius = pow(3.0*rmass[i]/(4.0*MY_PI*density[i]*neweffvolfactor),0.33333);
                // particles can gradually deform but won't recover state of higher sphericity
                if (newradius < radius[i] && newradius >= rmin_)
                {
                    radius[i] = newradius;
                    fix_effvolfactors_->set_vector(i,neweffvolfactor);
                }

            }
            else if (deformation >= 1.0 - SMALL)
            {
                nremoved_this_me++;
                mass_removed_this_me += rmass[i];

                if(monitor_heat_)
                {
                    double Cp = fix_capacity_->compute_vector(type[i]-1);
                    heat_removed_this_me += rmass[i]*T[i]*Cp;
                }

                delete_particle(i);
            }
        }
    }

    MPI_Allreduce(&mass_removed_this_me,&mass_removed_this,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&nremoved_this_me,&nremoved_this,1,MPI_INT,MPI_SUM,world);
    mass_removed_ += mass_removed_this;

    if(monitor_heat_)
    {
        MPI_Allreduce(&heat_removed_this_me,&heat_removed_this,1,MPI_DOUBLE,MPI_SUM,world);
        heat_removed_ += heat_removed_this;
    }


    if(comm->me == 0)
    {
        if(verbose_ && screen)
            fprintf(screen,"    Particle removal due to deformation:  n = %d, m = %f, m_tot = %f\n",
                    nremoved_this, mass_removed_this,mass_removed_);
    }

    if (atom->molecular == 0 && compress_flag_) {
      int *tag = atom->tag;
      for (int i = 0; i < atom->nlocal; i++) tag[i] = 0;
      atom->tag_extend();
    }

    if (atom->tag_enable) {
      if (atom->map_style) {
        atom->nghost = 0;
        atom->map_init();
        atom->map_set();
      }
    }
}

void FixCfdCouplingDeform::delete_particle(int i)
{
    atom->avec->copy(atom->nlocal-1,i,1);
    atom->nlocal--;
}


/* ----------------------------------------------------------------------
   provide accumulated removed mass and optionally heat
------------------------------------------------------------------------- */

double FixCfdCouplingDeform::compute_vector(int n)
{
  if (n==0)
  {
    return mass_removed_;
  }
  else if (n==1 && monitor_heat_)
  {
    return heat_removed_;
  }

  return -1.0;
}
