/* ----------------------------------------------------------------------
   LIGGGHTS academic

   Copyright 2005- JKU Linz

   LIGGGHTS academic is based on LIGGGHTS
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations
   Christoph Kloss, christoph.kloss@cfdem.com

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
   thomas.lichtenegger@jku.at
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "region.h"
#include "domain.h"
#include "neighbor.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_property_atom_cumulativetracer.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixPropertyAtomCumulativeTracer::FixPropertyAtomCumulativeTracer(LAMMPS *lmp, int narg, char **arg,bool parse) :
  FixPropertyAtom(lmp, narg, arg, false),
  iarg_(3),
  source_strength_(0),
  accumulated_source_strength_(0),
  absorbed_strength_(0),
  begin_time_(0),
  end_time_(0),
  tot_n_marked_(0),
  check_every_(1),
  iregion_(-1),
  idregion_(0)
{
    // do the base class stuff

    int n = strlen(id) + 1;
    tracer_name_ = new char[n];
    strcpy(tracer_name_,id);
    const char *baseargs[9];
    baseargs[0] = tracer_name_; //NP is actually not used
    baseargs[1] = "all";
    baseargs[2] = "property/atom/cumulativetracer";
    baseargs[3] = tracer_name_;
    baseargs[4] = "scalar"; //NP 1 scalar per particle to be registered
    baseargs[5] = "yes";    //NP restart yes
    baseargs[6] = "yes";    //NP communicate ghost yes
    baseargs[7] = "no";    //NP communicate rev no
    baseargs[8] = "0.";
    parse_args(9,(char**)baseargs);

    // do the derived class stuff

    // parse args

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg_],"region_mark") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'region_mark'");
            iarg_++;
            n = strlen(arg[iarg_]) + 1;
            idregion_ = new char[n];
            strcpy(idregion_,arg[iarg_]);
            iregion_ = domain->find_region(arg[iarg_++]);
            if (iregion_ == -1)
                error->fix_error(FLERR,this,"Region ID does not exist");
            hasargs = true;
        } else if(strcmp(arg[iarg_],"source_strength") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'source_strength'");
            iarg_++;
            source_strength_ = atof(arg[iarg_++]);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"begin_time") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'begin_time'");
            iarg_++;
            begin_time_ = atof(arg[iarg_++]);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"end_time") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'end_time'");
            iarg_++;
            end_time_ = atof(arg[iarg_++]);
            hasargs = true;
        } else if(strcmp(arg[iarg_],"check_mark_every") == 0) {
            if(narg < iarg_+2)
                error->fix_error(FLERR,this,"not enough arguments for 'check_mark_every'");
            iarg_++;
            check_every_ = atoi(arg[iarg_]);
            if(check_every_ < 0)
                error->fix_error(FLERR,this,"check_mark_every > 0 required");
            iarg_++;
            hasargs = true;
        } else if(strcmp(style,"property/atom/tracer") == 0 )
            error->fix_error(FLERR,this,"unknown keyword");
    }

    // error checks for this class

    if(strcmp(style,"property/atom/tracer") == 0)
    {
        if (iregion_ == -1)
            error->fix_error(FLERR,this,"expecting keyword 'region_mark'");
    }

    // settings
    nevery = check_every_;

    scalar_flag = 1;
    global_freq = 1;

    // distribute strength per unit time on time steps
    source_strength_*=update->dt*nevery;
}

/* ---------------------------------------------------------------------- */

FixPropertyAtomCumulativeTracer::~FixPropertyAtomCumulativeTracer()
{
    delete []tracer_name_;
    if(idregion_) delete []idregion_;
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixPropertyAtomCumulativeTracer::init()
{
    iregion_ = domain->find_region(idregion_);
    if (iregion_ == -1)
        error->fix_error(FLERR,this,"Region ID does not exist");
}

/* ---------------------------------------------------------------------- */

int FixPropertyAtomCumulativeTracer::setmask()
{
    int mask = FixPropertyAtom::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ----------------------------------------------------------------------
   mark particles
------------------------------------------------------------------------- */

void FixPropertyAtomCumulativeTracer::end_of_step()
{
    int ts = update->ntimestep;
    double dt = update->dt;
    double t=ts*dt;

    if(t < begin_time_ || t > end_time_)
        return;

    /*NL*/ //fprintf(screen,"FixPropertyAtomTracer::end_of_step(), proc %d, step " BIGINT_FORMAT "\n",comm->me,update->ntimestep);

    //NP mark all particles in region
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double *marker = this->vector_atom;
    Region *region = domain->regions[iregion_];

    int nmarked_this = 0;

    accumulated_source_strength_+=source_strength_;

    for(int i = 0; i < nlocal; i++)
    {
        if (region->match(x[i][0],x[i][1],x[i][2]))
        {
            nmarked_this++;
        }
    }

    MPI_Sum_Scalar(nmarked_this,world);

    if(nmarked_this<1)
        return;

    int newnmarked = 0;
    for(int i = 0; i < nlocal; i++)
    {
        if (region->match(x[i][0],x[i][1],x[i][2]))
        {
            if(marker[i] < 1e-7)
                newnmarked++;
            marker[i] += accumulated_source_strength_/nmarked_this;
        }
    }

    MPI_Sum_Scalar(newnmarked,world);
    tot_n_marked_+=newnmarked;
    absorbed_strength_ += accumulated_source_strength_;
    accumulated_source_strength_=0.0;
}

/* ----------------------------------------------------------------------
   return # marked particles
------------------------------------------------------------------------- */

double FixPropertyAtomCumulativeTracer::compute_scalar()
{
    return static_cast<double>(tot_n_marked_);
}
