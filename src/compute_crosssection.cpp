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
#include "string.h"
#include "stdlib.h"
#include "compute_crosssection.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "modified_andrew.h"

using namespace LAMMPS_NS;
using MODIFIED_ANDREW_AUX::Circle;

/* ---------------------------------------------------------------------- */

ComputeCrosssection::ComputeCrosssection(LAMMPS *lmp, int narg, char **arg) :
  ComputeContactAtom(lmp, narg, arg),
  file_(0)
{
  if (narg != 15 && narg != 17)
    error->all(FLERR,"Illegal compute crosssection command, need exactly 15 or 17 args");

  int iarg = 5;

  if(strcmp(arg[iarg++],"dim"))
    error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'dim'");

  if(strcmp(arg[iarg],"x") == 0)
    dim_ = 0;
  else if(strcmp(arg[iarg],"y") == 0)
    dim_ = 1;
  else if(strcmp(arg[iarg],"z") == 0)
    dim_ = 2;
  else
    error->all(FLERR,"Illegal compute crosssection command, expecting 'x', 'y' or 'z' after 'dim'");
  iarg++;

  if(strcmp(arg[iarg++],"min"))
    error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'min'");
  min_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"max"))
    error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'max'");
  max_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"n_cuts"))
    error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'n_cuts'");
  n_cuts_ = atoi(arg[iarg++]);
  if(n_cuts_ < 2)
    error->all(FLERR,"Illegal compute crosssection command, 'n_cuts' >= 2 required");

  if(strcmp(arg[iarg++],"cut_thickness"))
    error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'cut_thickness'");
  cut_thickness_half_ = atof(arg[iarg++]) / 2.;

  if(narg == 17)
  {
      if(strcmp(arg[iarg++],"file"))
        error->all(FLERR,"Illegal compute crosssection command, expecting keyword 'file'");
      if(0 == comm->me)
      {
        file_ = fopen(arg[iarg++],"w");
        if(!file_)
            error->one(FLERR,"Illegal compute crosssection command, cannot open file");
      }
  }

  // calculate slices
  setup_cuts();

  vector_flag = 1;
  vector_flag = 1;
  size_vector = n_cuts_;

  vector = new double[n_cuts_];
}

/* ---------------------------------------------------------------------- */

ComputeCrosssection::~ComputeCrosssection()
{
    for(int i = 0; i < n_cuts_; i++)
        delete mod_andrew_[i];
    delete [] mod_andrew_;

    delete [] vector;

    if(file_) fclose(file_);
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::setup_cuts()
{
    cut_dist_ = (max_ - min_) / (n_cuts_-1);

    mod_andrew_ = new ModifiedAndrew*[n_cuts_];
    for(int i = 0; i < n_cuts_; i++)
        mod_andrew_[i] = new ModifiedAndrew(lmp);
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_vector()
{
  if(invoked_vector == update->ntimestep)
    return;

  invoked_vector = update->ntimestep;

  compute_peratom();

  for(int i = 0; i < n_cuts_; i++)
  {
      //NP this also clears list of contacts in mod_andrew[i]
      vector[i] = mod_andrew_[i]->area();
  }

  if(file_) write();
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  ComputeContactAtom::compute_peratom();

  compute_convex_hull();
}

/* ---------------------------------------------------------------------- */

inline int ComputeCrosssection::mod(double coo)
{
    double coo_shift = coo - min_;
    int result = round( coo_shift / cut_dist_);
    double remainder = coo_shift - static_cast<double>( result ) * cut_dist_;

    if(remainder > -cut_thickness_half_ && remainder < cut_thickness_half_)
        return result;
    return -1;
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::compute_convex_hull()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double *radius = atom->radius;
    double coo;
    int m;
    Circle circle;

    double mi = min_ - cut_thickness_half_;
    double ma = max_ + cut_thickness_half_;

    for(int i = 0; i < nlocal; i++)
    {
        coo = x[i][dim_];
        if(contact[i] >= 2 && coo > mi && coo < ma)
        {
            m = mod(coo);
            if(m >= 0)
            {
                circle.x = x[i][(1+dim_)%3];
                circle.y = x[i][(2+dim_)%3];
                circle.r = radius[i];
                mod_andrew_[m]->add_contact(circle);
            }
        }
    }
}


/* ---------------------------------------------------------------------- */

void ComputeCrosssection::write()
{
    for(int i = 0; i < n_cuts_; i++)
    {
        double coo = min_ + static_cast<double>(i)*2.*cut_thickness_half_;
        fprintf(file_,"%f %f %f\n",coo,vector[i],sqrt(vector[i]/M_PI));
    }
    fflush(file_);
}
