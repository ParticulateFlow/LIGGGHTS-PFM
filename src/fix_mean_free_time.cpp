/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations
   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com
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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "error.h"
#include "group.h"
#include "region.h"
#include "domain.h"
#include "mpi_liggghts.h"
#include "fix_mean_free_time.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMeanFreeTime::FixMeanFreeTime(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  check_every_(1),
  t_start_(0.0)
{
    // do the derived class stuff

    // parse args

    bool hasargs = true;
    int iarg = 3;

    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"check_every") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'check_every'");
            iarg++;
            check_every_ = atoi(arg[iarg]);
            if(check_every_ < 0)
                error->fix_error(FLERR,this,"check_every > 0 required");
            iarg++;
            hasargs = true;
        }
        else if(strcmp(arg[iarg],"t_start") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 't_start'");
            iarg++;
            t_start_ = atof(arg[iarg]);
            iarg++;
            hasargs = true;
        } else if(strcmp(style,"MeanFreeTime") == 0 )
            error->fix_error(FLERR,this,"unknown keyword");
    }

    nevery = check_every_;
    scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

FixMeanFreeTime::~FixMeanFreeTime()
{
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixMeanFreeTime::init()
{
 if (force->pair == NULL)
    error->all(FLERR,"Fix mean_free_time requires a pair style be defined");

  // need an occasional neighbor list

  int irequest = neighbor->request((void *) this);
 // neighbor->requests[irequest]->half = 0;
 // neighbor->requests[irequest]->gran = 1;
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
 // neighbor->requests[irequest]->occasional = 1;

fprintf(screen, "irequest = %d\n",irequest);

  const char* fixarg[12];
  // register property/atom for time between collisions
  fixarg[0]="meanfreetime";            // fixid
  fixarg[1]="all";			//group;
  fixarg[2]="property/atom";
  fixarg[3]="meanfreetime";           // propertyid
  fixarg[4]="vector";          // vector with 3 values: free time so far, step of most recent collision, number of intervals so far, contact at previous step (0 no, 1 yes)
  fixarg[5]="yes";             // restart yes
  fixarg[6]="no";              // communicate ghost - yes
  fixarg[7]="no";              // communicate rev no
  fixarg[8]="0.0";             // tale 0 as default value
  fixarg[9]="0.0";
  fixarg[10]="0.0";
  fixarg[11]="1.0";
  fix_meanfreetime_ = modify->add_fix_property_atom(12,const_cast<char**>(fixarg),style);

}

/* ---------------------------------------------------------------------- */

void FixMeanFreeTime::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

int FixMeanFreeTime::setmask()
{
    int mask = 0;
    mask |= END_OF_STEP;
    return mask;
}

/*----------------------------------------------------------------------- */

void FixMeanFreeTime::end_of_step()
{
    int i,j,ii,jj,inum,jnum;
    double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
    double radi,radsum,radsumsq;
    int *ilist,*jlist,*numneigh,**firstneigh;

    double **meanfreetime = fix_meanfreetime_->array_atom;

    bigint currstep = update->ntimestep;

    if (update->dt * currstep < t_start_) return;
    // invoke neighbor list (will copy or build if necessary)
//    neighbor->build_one(list->index);

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // compute number of contacts for each atom in group
    // contact if distance <= sum of radii
    // tally for both I and J

    double **x = atom->x;
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      if (mask[i] & groupbit) {
        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        radi = radius[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          j &= NEIGHMASK;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radsum = radi + radius[j];
          radsumsq = radsum*radsum;
          if (rsq <= radsumsq) {
            if (meanfreetime[i][3] < 0.5) {
              meanfreetime[i][0] += currstep - meanfreetime[i][1];
              meanfreetime[i][2] += 1.0;
              meanfreetime[i][3] = 1.0;
            }
          }
          else {
            if (meanfreetime[i][3] > 0.5) {
              meanfreetime[i][1] = currstep;
              meanfreetime[i][3] = 0.0;
            }
          }
        }
      }
    }
}

/* ----------------------------------------------------------------------
   return average of mean free time
------------------------------------------------------------------------- */

double FixMeanFreeTime::compute_scalar()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double **meanfreetime = fix_meanfreetime_->array_atom;

  double value = 0.;

  for (int i = 0; i < nlocal; i++)
  {
      if ((mask[i] & groupbit) & (meanfreetime[i][2] > 0.5))
      {
          value += meanfreetime[i][0]/meanfreetime[i][2];
      }
  }

  value *= update->dt;

  MPI_Sum_Scalar(value,world);
  return value/atom->natoms;
}

