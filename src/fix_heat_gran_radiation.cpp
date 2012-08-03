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
  //NP this is handled in init()
	TB = -1.0;

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

  if (TB == -1.0){
    TB = T0;
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::post_force(int vflag){

}

/* ---------------------------------------------------------------------- */
/*
finds the intersection point of a ray and a sphere, if it exists.
o ......... origin of ray
d ......... direction of ray
center .... center of sphere
radius .... radius of sphere
t ......... return value - intersection point (of parameterized ray intc = o + t*d)
buffer3 ... undefined return value - array of length three that intersectRaySphere() can use for calculations
*/
//NP unit tested
bool FixHeatGranRad::intersectRaySphere(double *o, double *d, double *center, double radius, double &t, double *buffer3){

  double A, B, C;
  double discr, q;
  double t0, t1, ttemp;

  t = 0.0;

  A = lensq3(d);
  sub3(o, center, buffer3);
  B = 2.0 * dot3(buffer3, d);
  C = lensq3(buffer3) - radius*radius;

  discr = B*B - 4.0*A*C;

  // miss
  if (discr < 0.0){
    t = 0.0;
    return false;
  }

  // hit
  q = 0.0;
  if (B < 0.0){
    q = (-B - sqrt(discr)) * 0.5;
  } else {
    q = (-B + sqrt(discr)) * 0.5;
  }

  // calculate ts
  t0 = q / A;
  t1 = C / q;

  // sort ts
  if (t0 > t1){
    ttemp = t0;
    t0 = t1;
    t1 = ttemp;
  }

  // if t1 is less than zero, the object is in the ray's negative direction
  // and consequently the ray misses the sphere
  if (t1 <= 0.0){
    t = 0.0;
    return false;
  }

  // if t0 is less than zero, the intersection point is at t1
  if (t0 < 0.0){
    t = t1;
  } else {
    t = t0;
  }

  return true;
}

/* ---------------------------------------------------------------------- */
/*
generates a random point on the surface of a sphere
c ...... center of sphere
r ...... radius of sphere
ansP ... return value - random point on sphere
ansD ... return value - direction vector that points out of the sphere
*/
//NP unit tested
void FixHeatGranRad::randOnSphere(double *c, double r, double *ansP, double *ansD){

  // generate random direction
  ansD[0] = 2.0*rGen.uniform() - 1.0;
  ansD[1] = 2.0*rGen.uniform() - 1.0;
  ansD[2] = 2.0*rGen.uniform() - 1.0;

  normalize3(ansD, ansD);

  // calculate corresponding point on surface
  snormalize3(r, ansD, ansP);

  ansP[0] += c[0];
  ansP[1] += c[1];
  ansP[2] += c[2];

}

/* ---------------------------------------------------------------------- */
/*
  generates a random direction in direction of the vector n
  n ... normal vector
  o ... return value - random direction
*/
//NP unit tested
void FixHeatGranRad::randDir(double *n, double *o){
  double side;

  // generate random direction
  o[0] = 2.0*rGen.uniform() - 1.0;
  o[1] = 2.0*rGen.uniform() - 1.0;
  o[2] = 2.0*rGen.uniform() - 1.0;
  normalize3(o, o);

  // adjust direction, if it points to the wrong side
  // and if n != 0
  side = dot3(o, n);
  if (side < 0.0 && (n[0] != 0.0 || n[1] != 0.0 || n[2] != 0.0)){
    o[0] *= -1.0;
    o[1] *= -1.0;
    o[2] *= -1.0;
  }
}
