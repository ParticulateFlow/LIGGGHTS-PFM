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

/* ----------------------------------------------------------------------
   Contributing authors:
   Thomas Lichtenegger (JKU Linz)
------------------------------------------------------------------------- */


#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "comm.h"
#include "math.h"
#include "math_extra.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "fix_cfd_coupling_recurrence.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingRecurrence::FixCfdCouplingRecurrence(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg,arg),
    fix_coupling_(0),
    fix_vrec_(0),
    fix_vfluc_(0),
    fix_dragforce_(0),
    fix_volumeweight_(0),
    fix_tracerconcentration_(0),
    use_force_(false),
    use_fluc_(false),
    use_dens_(false),
    use_type_(false),
    use_tracer_(false),
    use_property_(false),
    limit_fluc(false),
    relative_limit(false),
    remove_vel_across_walls(false),
    iregion(-1)
{
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"transfer_density") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_density'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_dens_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_dens_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_density'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_type") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_type_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_type_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_type'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_force") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_force'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_force_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_force_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_force'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_fluctuations") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_fluctuations'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_fluc_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_fluc_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_fluctuations'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_tracer") == 0)
        {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_tracer'");
            iarg++;
            if(strcmp(arg[iarg],"yes") == 0)
                use_tracer_ = true;
            else if(strcmp(arg[iarg],"no") == 0)
                use_tracer_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'transfer_tracer'");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"transfer_property") == 0)
        {
            if(narg < iarg+5)
                error->fix_error(FLERR,this,"not enough arguments for 'transfer_type'");
            iarg++;
            use_property_ = true;
            if(strcmp(arg[iarg++],"name"))
                error->fix_error(FLERR,this,"expecting 'name' after 'transfer_property'");
            sprintf(property_name,"%s",arg[iarg++]);
            if(strcmp(arg[iarg++],"type"))
                error->fix_error(FLERR,this,"expecting 'type' after property name");
            sprintf(property_type,"%s",arg[iarg++]);
            iarg++;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"limit_fluctuations") == 0)
        {
            if (narg < iarg+2) error->all(FLERR,"Illegal fix couple/cfd/recurrence command");
            if(strcmp(arg[iarg+1],"yes") == 0)
                limit_fluc = true;
            else if(strcmp(arg[iarg+1],"no"))
                error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'limit_fluctuations'");
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"max_vfluc") == 0)
        {
            if (narg < iarg+2) error->all(FLERR,"Illegal fix couple/cfd/recurrence command");
            maxvfluc = atof(arg[iarg+1]);
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"relative_limit") == 0)
        {
            if (narg < iarg+2) error->all(FLERR,"Illegal fix couple/cfd/recurrence command");
            if(strcmp(arg[iarg+1],"yes") == 0)
                relative_limit = true;
            else if(strcmp(arg[iarg+1],"no"))
                error->fix_error(FLERR,this,"expecing 'yes' or 'no' for 'relative_limit'");
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"remove_vel_across_walls") == 0)
        {
            if (narg < iarg+3) error->all(FLERR,"Illegal fix couple/cfd/recurrence command");
            remove_vel_across_walls = true;
            int n = strlen(arg[iarg+1]) + 1 + 6;
            wallfixname = new char[n];
            strcpy(wallfixname,"force_");
            strcat(wallfixname,arg[iarg+1]);
            fwcrit = atof(arg[iarg+2]);
            iarg += 3;
            hasargs = true;
        }
        else if (strcmp(arg[iarg],"region") == 0)
        {
            if (narg < iarg+2) error->all(FLERR,"Illegal fix couple/cfd/recurrence command");
            iregion = domain->find_region(arg[iarg+1]);
            if (iregion == -1)
                error->all(FLERR,"Region ID for fix couple/cfd/recurrence does not exist");
            int n = strlen(arg[iarg+1]) + 1;
            idregion = new char[n];
            strcpy(idregion,arg[iarg+1]);
            iarg += 2;
            hasargs = true;
        }
        else if (strcmp(this->style,"couple/cfd/recurrence") == 0)
        {
            error->fix_error(FLERR,this,"unknown keyword");
        }
    }

    // flags for vector output
    vector_flag = 1;
    size_vector = 3;
    global_freq = 1;
    extvector = 1;
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingRecurrence::~FixCfdCouplingRecurrence()
{
  if(wallfixname) delete []wallfixname;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::post_create()
{
    // register vrec
    if(!fix_vrec_)
    {
        const char* fixarg[11];
        fixarg[0]="vrec";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="vrec";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_vrec_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    // register displacements
    if(!fix_vfluc_ && use_fluc_)
    {
        const char* fixarg[11];
        fixarg[0]="vfluc";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="vfluc";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_vfluc_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    // register dragforce
    if(!fix_dragforce_ && use_force_)
    {
        const char* fixarg[11];
        fixarg[0]="dragforce";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="dragforce";
        fixarg[4]="vector"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_dragforce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
    }

    // register volume weight for volume fraction calculation if not present
    // is 1 per default
    fix_volumeweight_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("volumeweight","property/atom","scalar",0,0,style,false));
    if(!fix_volumeweight_)
    {
        const char* fixarg[9];
        fixarg[0]="volumeweight";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="volumeweight";
        fixarg[4]="scalar"; // 1 scalar per particle to be registered
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="1.";
        fix_volumeweight_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

     // register tracer concentration
    if(!fix_tracerconcentration_ && use_tracer_)
    {
        const char* fixarg[11];
        fixarg[0]="tracerconcentration";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="tracerconcentration";
        fixarg[4]="scalar"; // 1 scalar per particle to be registered
        fixarg[5]="yes";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_tracerconcentration_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::pre_delete(bool unfixflag)
{
    if(unfixflag && fix_vrec_) modify->delete_fix("vrec");
    if(unfixflag && fix_vfluc_) modify->delete_fix("vfluc");
    if(unfixflag && fix_dragforce_) modify->delete_fix("dragforce");
    if(unfixflag && fix_volumeweight_) modify->delete_fix("volumeweight");
    if(unfixflag && fix_tracerconcentration_) modify->delete_fix("tracerconcentration");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingRecurrence::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"Fix couple/cfd/recurrence needs a fix of type couple/cfd");

    //  values to be transfered to OF

    fix_coupling_->add_push_property("x","vector-atom");
    fix_coupling_->add_push_property("v","vector-atom");
    fix_coupling_->add_push_property("radius","scalar-atom");
    if(use_type_) fix_coupling_->add_push_property("type","scalar-atom");
    if(use_dens_) fix_coupling_->add_push_property("density","scalar-atom");
    fix_coupling_->add_push_property("volumeweight","scalar-atom");

    if(use_property_) fix_coupling_->add_push_property(property_name,property_type);

    // values to come from OF
    fix_coupling_->add_pull_property("vrec","vector-atom");
    if(use_force_) fix_coupling_->add_pull_property("dragforce","vector-atom");
    if(use_fluc_) fix_coupling_->add_pull_property("vfluc","vector-atom");
    if(use_tracer_) fix_coupling_->add_pull_property("tracerconcentration","scalar-atom");


    vectorZeroize3D(dragforce_total);
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::initial_integrate(int)
{
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double **vrec = fix_vrec_->array_atom;
  double **fwall = NULL;

  double v_conv[3];
  double v_fluc[3];

  if (remove_vel_across_walls)
  {
    class FixPropertyAtom* fix_fwall = static_cast<FixPropertyAtom*>(modify->find_fix_property(wallfixname,"property/atom","vector",3,0,style,false));
    if (!fix_fwall)
        error->fix_error(FLERR,this,"Fix couple/cfd/recurrence could not find fix for particle-wall forces with provided name");
    fwall = fix_fwall->array_atom;
  }

  // displace particles if necessary
  if(use_fluc_)
  {
    double *rad = atom->radius;
    double **vfluc = fix_vfluc_->array_atom;
    double **x = atom->x;
    double x_new[3];
    double dt = update->dt;

    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        for (int j = 0; j < 3; j++)
        {
          v_fluc[j] = vfluc[i][j];
        }

        if (limit_fluc)
        {
          limit_vfluc(x[i],v_fluc,rad[i]/dt);
        }

        if (remove_vel_across_walls)
        {
          correct_vel_across_walls(v_fluc,fwall[i]);
        }

        for(int j=0; j<3; j++)
        {
          x_new[j] = x[i][j] + v_fluc[j] * dt;
          if (x_new[j] <= domain->boxlo[j] || x_new[j] >= domain->boxhi[j])
          {
            // if outside domain either a) reflect fluctuations or b) ignore them
            // if a), this assumes that fluctuations do not put particles outside domain
            // in two opposite directions at the same time

            //x_new[j] = x[i][j] - vfluc[i][j] * dt;
            x_new[j] = x[i][j];
          }
        }

        vectorCopy3D(x_new,x[i]);
      }
    }
  }

  // set particle velocity to that of recurrence field
  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
      for (int j = 0; j < 3; j++)
      {
        v_conv[j] = vrec[i][j];
      }

      if (remove_vel_across_walls)
      {
        correct_vel_across_walls(v_conv,fwall[i]);
      }

      vectorCopy3D(v_conv,v[i]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::post_force(int)
{
  if(use_force_)
  {
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double **dragforce = fix_dragforce_->array_atom;

    vectorZeroize3D(dragforce_total);

    // add dragforce to force vector
    //NP no allreduce here, is done in compute_vector
    for (int i = 0; i < nlocal; i++)
    {
      if (mask[i] & groupbit)
      {
        //vectorAdd3D(f[i],dragforce[i],f[i]);
        // I'm the only force!
        vectorCopy3D(dragforce[i],f[i]);
        vectorAdd3D(dragforce_total,dragforce[i],dragforce_total);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingRecurrence::limit_vfluc(double* pos, double* vfluc, double relativescale)
{
  double vmax = maxvfluc;

  if (iregion >= 0 && !domain->regions[iregion]->match(pos[0],pos[1],pos[2]))
  {
    return;
  }

  if (relative_limit)
  {
    vmax = maxvfluc * relativescale;
  }

  double absvfluc = MathExtra::len3(vfluc);
  if (absvfluc > vmax)
  {
    MathExtra::scale3(vmax/absvfluc,vfluc);
  }
}

void FixCfdCouplingRecurrence::correct_vel_across_walls(double* vel, double* fw)
{
  double fw2 = MathExtra::dot3(fw,fw);
  if (fw2 < fwcrit*fwcrit) return;

  double n[3];
  n[0] = fw[0];
  n[1] = fw[1];
  n[2] = fw[2];

  MathExtra::norm3(n);

  double vn = MathExtra::dot3(vel,n);

  // correct only if velocity has components towards the wall
  if (vn > 0) return;

  // reflect velocity component, hence factor 2
  vel[0] -= 2*vn*n[0];
  vel[1] -= 2*vn*n[1];
  vel[2] -= 2*vn*n[2];
}

