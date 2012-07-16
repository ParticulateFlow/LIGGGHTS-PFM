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

#include "fix_heat_gran_radiation.h"

#include <cstdlib>
#include <cstring>
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

//NP assumptions:
//NP if background_temperature is not passed it is assumed that
//NP    initial_temperature is the background temperature
FixHeatGranRad::FixHeatGranRad(class LAMMPS *lmp, int narg, char **arg) : FixHeatGran(lmp, narg, arg){
	int iarg = 5;

	//NP use initial temperature as background temperature if
	//NP TB is not passed as argument background_temperature
	TB = T0;

	bool hasargs = true;
	while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"background_temperature") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'background_temperature'");
      TB = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if(strcmp(style,"heat/gran/radiation") == 0)
      	error->fix_error(FLERR,this,"unknown keyword");
  }
}

/* ---------------------------------------------------------------------- */

FixHeatGranRad::~FixHeatGranRad(){}

/* ---------------------------------------------------------------------- */

int FixHeatGranRad::setmask(){
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::init(){

	if (FHG_init_flag == false){
		FixHeatGran::init();
	}

  //NP Get pointer to all the fixes (also those that have the material properties)
  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::post_force(int vflag){

}

/* ---------------------------------------------------------------------- */





















