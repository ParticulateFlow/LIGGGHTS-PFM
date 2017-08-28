/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Stefan Radl, radl@tugraz.at
   Copyright 2013-  Graz University of Technology (TU Graz)

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

#include "fix_roughness.h"

#include "atom.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "compute_pair_gran_local.h"
#include "neigh_list.h"
#include "force.h"
#include "group.h"
#include "math_extra.h"
#include "math_const.h"
#include "modify.h"
#include "pair_gran.h"
#include <stdlib.h>
#include "random_park.h"
#include "comm.h"


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRoughness::FixRoughness(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg){

  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all(FLERR,"Fix liquidTrack needs per particle radius and mass");

  if (narg < 7)
    error->fix_error(FLERR,this,"not enough arguments");

  int iarg = 3;

  if(strcmp(arg[iarg++],"roughMean"))
    error->fix_error(FLERR,this,"roughMean");
  roughMean = atof(arg[iarg++]);

  if(strcmp(arg[iarg++],"deltaGammaMean"))
    error->fix_error(FLERR,this,"deltaGammaMean");
  deltaGammaMean = atof(arg[iarg++]) * MathConst::MY_PI / 180.0;
  if(deltaGammaMean<=0.0)
     error->fix_error(FLERR,this,"deltaGammaMean is too small. deltaGammaMean must be larger than 0°.");
  if(deltaGammaMean>20*MathConst::MY_PI / 180.0)
     error->warning(FLERR,"deltaGammaMean is too large. deltaGammaMean should not exceed 20°.");

  if(strcmp(arg[iarg++],"seed"))
    error->fix_error(FLERR,this,"seed'");
  seed = atoi(arg[iarg++]);


  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  time_depend = 1;

  scalar_flag = 1;
  global_freq = 1;

  haveNonZeroDeltaGamma_ = false;
  if(deltaGammaMean>0.0) haveNonZeroDeltaGamma_= true;

  cpl = NULL;

  FHG_init_flag = false;

  // random number generator, different for CPUs
  random_equal = new RanPark(lmp,seed + 3000 * comm->me);

}
FixRoughness::~FixRoughness()
{
  delete random_equal;
}
/* ---------------------------------------------------------------------- */

void FixRoughness::post_create()
{

}

/* ---------------------------------------------------------------------- */
void FixRoughness::pre_delete(bool unfixflag)
{

  // tell cpl that this fix is deleted
  if(cpl && unfixflag) cpl->reference_deleted();

}

/* ---------------------------------------------------------------------- */

void FixRoughness::updatePtrs(){

  vector_atom = liqOnParticle;


}

/* ---------------------------------------------------------------------- */

void FixRoughness::init(){

  if (!atom->radius_flag || !atom->rmass_flag) error->all(FLERR,"Please use a granular atom style for fix liq/gran");

    // check if a fix of this style already exists
  if(modify->n_fixes_style(style) > 1)
    error->fix_error(FLERR,this,"cannot have more than one fix of this style");

  if(!force->pair_match("gran", 0)) error->all(FLERR,"Please use a granular pair style for fix liq/gran");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();


  dnum = pair_gran->dnum_pair();
  dnum_mine = pair_gran->fix_extra_dnum_index(this);

  updatePtrs();

  FHG_init_flag = true;
}

/* ---------------------------------------------------------------------- */

double FixRoughness::compute_scalar()
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */

int FixRoughness::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixRoughness::n_history_extra() const
{
    return 2+3+3; //2 values for orientation, 3 for contact point, 3 for first normal vector of plane
}

/* ---------------------------------------------------------------------- */

bool FixRoughness::history_args(char** args) const
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value is same
    args[0] = (char *) "deltaGamma";
    args[1] = (char *) "0";
    args[2] = (char *) "psi";
    args[3] = (char *) "0";
    args[4] = (char *) "contactPX";
    args[5] = (char *) "0";
    args[6] = (char *) "contactPY";
    args[7] = (char *) "0";
    args[8] = (char *) "contactPZ";
    args[9] = (char *) "0";
    args[10] = (char *) "contactN1X";
    args[11] = (char *) "0";
    args[12] = (char *) "contactN1Y";
    args[13] = (char *) "0";
    args[14] = (char *) "contactN1Z";
    args[15] = (char *) "0";

    return true;
}

/* ---------------------------------------------------------------------- */
double FixRoughness::generateDeltaGamma(double gamma)
{
    double deltaGamma = 0.0;
    do
    {
        deltaGamma = random_equal->gaussian() * deltaGammaMean;
    } while( (gamma/2.0+deltaGamma) < 0.0);

    if(deltaGamma>MathConst::MY_PI/2.0)
        deltaGamma=MathConst::MY_PI/2.0 - gamma;

    return deltaGamma;
}

/* ---------------------------------------------------------------------- */
double FixRoughness::generatePsi()
{
    double psi;
    psi = MathConst::MY_PI * random_equal->uniform() ;

    return psi;
}

/* ---------------------------------------------------------------------- */

void FixRoughness::initial_integrate(int vflag)
{
  int i,j,ii,jj,inum,jnum;

  //reset heat flux
  //sources are not reset
  int *mask = atom->mask;


  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touchpair,**firsttouch;
  double *allhist, **firsthist;
  double *deltaGamma, *psi;

  updatePtrs();

  //Pull the contact history
  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  firsttouch = pair_gran->listgranhistory->firstneigh;
  firsthist = pair_gran->listgranhistory->firstdouble;

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    jnum = numneigh[i];
    jlist = firstneigh[i];
    touchpair = firsttouch[i];

    allhist = firsthist[i];
    for (jj = 0; jj < jnum; jj++)
    {
      j = jlist[jj];
      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      //reset the extra values if particles not longer touch
      deltaGamma = &allhist[(dnum+n_history_extra())*jj+3]; //the deltaGamma of the contact with jj
      psi        = &allhist[(dnum+n_history_extra())*jj+4];   //the psi value of contact with jj
      if (!touchpair[jj] )
      {
            deltaGamma[0] = 0.0;
            psi[0]        = 0.0;
            psi[1]        = 0.0; //contact point coordinates
            psi[2]        = 0.0; //contact point coordinates
            psi[3]        = 0.0; //contact point coordinates
      }
      //DEBUG: Save history to file
#if 0
	  history = fopen("roughnessHistory.txt","a+");
	  fprintf(history,"i/j: %d %d,dnum/n_history_extra: %d %d, touchpair[jj] %d, allhistory %f %f %f , deltaGamma/psi %f %f \n",
				   i, j,
				   dnum, n_history_extra(),
				   touchpair[jj],
				   allhist[(dnum+n_history_extra())*jj], allhist[(dnum+n_history_extra())*jj+1], allhist[(dnum+n_history_extra())*jj+2], deltaGamma[0], psi[0]);
	  fclose(history);
#endif
    }

  }



}


/* ---------------------------------------------------------------------- */

void FixRoughness::post_force(int vflag){

  //template function for using touchflag or not
//  if(history_flag == 0) post_force_eval<0>(vflag,0);
//  if(history_flag == 1) post_force_eval<1>(vflag,0);

}




/* ---------------------------------------------------------------------- */

void FixRoughness::register_compute_pair_local(class ComputePairGranLocal *ptr)
{

   if(cpl != NULL)
      error->all(FLERR,"Fix liquidTracking/instant allows only one compute of type pair/local");
   cpl = ptr;

}

/* ---------------------------------------------------------------------- */

void FixRoughness::unregister_compute_pair_local(class ComputePairGranLocal *ptr)
{

   if(cpl != ptr)
       error->all(FLERR,"Illegal situation in FixLiquidTrackInst::unregister_compute_pair_local");
   cpl = NULL;

}
