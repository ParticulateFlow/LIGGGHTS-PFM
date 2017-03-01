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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "compute_crosssection.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "domain.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "region.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "modified_andrew.h"
#include "math_extra_liggghts.h"

#if defined(_WIN32) || defined(_WIN64)
double inline round(double d);
#endif

using namespace LAMMPS_NS;
using MODIFIED_ANDREW_AUX::Circle;

/* ---------------------------------------------------------------------- */

ComputeCrosssection::ComputeCrosssection(LAMMPS *lmp, int narg, char **arg) :
  ComputeContactAtom(lmp, narg, arg),
  angle_(0),
  file_(0),
  iregion_(-1),
  idregion_(0)
{
  if (narg < 15)
    error->compute_error(FLERR,this,"need at least 15 args");

  int iarg = 5;

  if(strcmp(arg[iarg++],"dim"))
    error->compute_error(FLERR,this,"expecting keyword 'dim'");

  if(strcmp(arg[iarg],"x") == 0)
    dim_ = 0;
  else if(strcmp(arg[iarg],"y") == 0)
    dim_ = 1;
  else if(strcmp(arg[iarg],"z") == 0)
    dim_ = 2;
  else
    error->compute_error(FLERR,this,"expecting 'x', 'y' or 'z' after 'dim'");
  iarg++;

  if(strcmp(arg[iarg++],"min"))
    error->compute_error(FLERR,this,"expecting keyword 'min'");
  min_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"max"))
    error->compute_error(FLERR,this,"expecting keyword 'max'");
  max_ = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"n_cuts"))
    error->compute_error(FLERR,this,"expecting keyword 'n_cuts'");
  n_cuts_ = atoi(arg[iarg++]);
  if(n_cuts_ < 2)
    error->compute_error(FLERR,this,"'n_cuts' >= 2 required");

  if(strcmp(arg[iarg++],"cut_thickness"))
    error->compute_error(FLERR,this,"expecting keyword 'cut_thickness'");
  cut_thickness_half_ = atof(arg[iarg++]) / 2.;

  // parse args
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
      hasargs = false;
      if(strcmp(arg[iarg],"file") == 0) {
          if(narg < iarg+2)
              error->compute_error(FLERR,this,"not enough arguments for 'file'");
          if(0 == comm->me)
          {
            file_ = fopen(arg[iarg+1],"w");
            if(!file_)
                error->one(FLERR,"Illegal compute crosssection command, cannot open file");
          }
          iarg += 2;
          hasargs = true;
      } else if(strcmp(arg[iarg],"region") == 0) {
            if(narg < iarg+2)
                error->compute_error(FLERR,this,"not enough arguments for 'region'");
            int iregion = domain->find_region(arg[iarg+1]);
            if (iregion == -1)
                error->compute_error(FLERR,this,"region ID does not exist");
            int n = strlen(arg[iarg+1]) + 1;
            idregion_ = new char[n];
            strcpy(idregion_,arg[iarg+1]);
            iarg += 2;
            hasargs = true;
      }
      else error->compute_error(FLERR,this,"unknown keyword");
  }

  // calculate slices
  setup_cuts();

  vector_flag = 1;

  if(file_)
    size_vector = n_cuts_+1;
  else
    size_vector = n_cuts_;

  vector = new double[n_cuts_+1];
}

/* ---------------------------------------------------------------------- */

ComputeCrosssection::~ComputeCrosssection()
{
    for(int i = 0; i < n_cuts_; i++)
        delete mod_andrew_[i];
    delete [] mod_andrew_;

    delete [] vector;

    if(file_) fclose(file_);

    if(idregion_) delete []idregion_;
}

/* ---------------------------------------------------------------------- */

void ComputeCrosssection::init()
{
  ComputeContactAtom::init();

  // set index and check validity of region

  if (idregion_) {
    iregion_ = domain->find_region(idregion_);
    if (iregion_ == -1)
      error->compute_error(FLERR,this,"region ID does not exist");
  }
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

  if(file_)
  {
      vector[n_cuts_] = calc_ang();
      //fprintf(file_,"%f\n",ang);
      write();
  }
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
           if (iregion_ >= 0 && !domain->regions[iregion_]->match(x[i][0],x[i][1],x[i][2]))
             continue;

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
    //NP write to text file
    //NP format: coordinate area r*
    //NP where r* is the radius of area-equivalent cricle

    for(int i = 0; i < n_cuts_; i++)
    {
        double coo = min_ + static_cast<double>(i)*cut_dist_;
        fprintf(file_,"%f %f %f\n",coo,vector[i],sqrt(vector[i]/M_PI));
    }
    fflush(file_);
}

/* ---------------------------------------------------------------------- */

double ComputeCrosssection::calc_ang()
{
    //NP calculates angle assuming areas are forming a cone

    int ilo = 0;
    int ihi = n_cuts_-1;
    int imid, ilomid, ihimid, idelta;
    double rhimid, rlomid, del, ang;

    //NP calc ilo,ihi

    for(int i = 0; i < n_cuts_; i++)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,1e-6))
            ilo += 1;
        else
            break;
    }

    for(int i = n_cuts_-1; i >= 0; i--)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,1e-6))
            ihi -= 1;
        else
            break;
    }

    //NP error checks ilo,ihi

    if(ilo == ihi)
        error->one(FLERR,"Compute crossection could not calculate angle (1)");

    for(int i = ilo; i <= ihi; i++)
    {
        if(MathExtraLiggghts::compDouble(vector[i],0.,1e-6))
            error->one(FLERR,"Compute crossection could not calculate angle - internal error");
    }

    //NP calc ilomid,ihimid


    imid = static_cast<int>(0.5*static_cast<double>(ilo+ihi));
    idelta = static_cast<int>(0.25*static_cast<double>(ihi-ilo));

    ilomid = imid-idelta;
    ihimid = imid+idelta;

    //NP check ilomid,ihimid
    //NP in case of flat heap, just take ilo and ihi

    if(ilomid == ihimid || ilomid == ilo || ihimid == ihi)
    {
        //ilomid = ilo;
        //ihimid = ihi;
        /*NL*/ if(screen) fprintf(screen,"ilo %d ilomid %d ihimid %d ihi %d\n",ilo,ilomid,ihimid,ihi);
        error->one(FLERR,"Compute crossection could not calculate angle (2)");
    }


    //NP calc angle

    rlomid = sqrt(vector[ilomid]/M_PI);
    rhimid = sqrt(vector[ihimid]/M_PI);
    del = (ihimid-ilomid)*cut_dist_;

    if(rhimid > rlomid)
        ang = 90. + 180./M_PI*atan((rhimid-rlomid)/del);
    else
        ang = 180./M_PI*atan(del/(rlomid-rhimid));

    /*NL*/ if(screen) fprintf(screen,"rlomid %f rhimid %f del %f ang %f\n",rlomid,rhimid,del,ang);

    return ang;
}
