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

#include "atom.h"
#include "error.h"
#include "fix_property_global.h"
#include "math_const.h"
#include "math_extra.h"
#include "mech_param_gran.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair_gran.h"
#include "random_mars.h"
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

using MathExtra::dot3;
using MathExtra::lensq3;
using MathExtra::sub3;
using MathExtra::normalize3;
using MathExtra::snormalize3;
using MathConst::MY_4PI;

/* ---------------------------------------------------------------------- */

//NP assumptions:
//NP if background_temperature is not passed it is assumed that
//NP    initial_temperature is the background temperature
FixHeatGranRad::FixHeatGranRad(class LAMMPS *lmp, int narg, char **arg) : FixHeatGran(lmp, narg, arg),
                                                                          RGen(lmp, 1){
	int iarg = 5;

	//NP use initial temperature as background temperature if
	//NP TB is not passed as argument background_temperature
  //NP background temperature is further handled in init()

  Qr = 0.0;
  avgNRays = 10; //NP TODO optimize this parameter
  maxBounces = 10;

  TB  = -1.0;
  Qtot = 0.0;

  fix_emissivity = NULL;
  emissivity = NULL;

  Sigma = 5.670373E-8;

  // parse input arguments:
  //  - backgroundTemperature
  //  - averageNumberOfRaysPerParticle
	bool hasargs = true;
	while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"backgroundTemperature") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'backgroundTemperature'");
      TB = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"averageNumberOfRaysPerParticle") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'backgroundTemperature'");
      avgNRays = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"maxNumberOfBounces") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'backgroundTemperature'");
      maxBounces = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(style,"heat/gran/radiation") == 0)
      	error->fix_error(FLERR,this,"unknown keyword");
  }
}

/* ---------------------------------------------------------------------- */

FixHeatGranRad::~FixHeatGranRad(){
  if (emissivity){
    delete [] emissivity;
  }
}

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

  // set background temperature to T0 if it was not provided as an argument
  if (TB == -1.0){
    TB = T0;
  }

  // find mechanical-parameter-granular "emissivity" and fetch data
  // if "emissivity" not found: error handling happens inside find_fix_property().
  int max_type = pair_gran->mpg->max_type();
  if (emissivity){
    delete [] emissivity;
  }
  emissivity = new double[max_type];
  fix_emissivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalEmissivity","property/global","peratomtype",max_type,0,style));
  for (int i = 0; i < max_type; i++)
  {
    emissivity[i] = fix_emissivity->compute_vector(i);
    if (emissivity[i] < 0.0 || emissivity[i] > 2.0){
      error->all(FLERR,"Fix heat/gran/radiation: Thermal emissivity must not be < 0.0 or > 1.0.");
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::post_force(int vflag){
  //NP neighborlist data
  int *ilist;
  int *numneigh;
  int i;
  int inum;

  //NP particle data
  double **x;
  double *radius;
  int *mask;
  int *type;
  int nlocal;

  //NP ray data
  double *buffer = new double[3]; //NP buffer for computations later on
  double *d = new double[3]; //NP direction of ray
  double *o = new double[3]; //NP origin of ray
  double t; //NP parameter of ray

  //NP individual particle data
  double areai; //NP area
  double emisi; //NP emissivity
  double radi;  //NP radius
  double tempi; //NP temperature

  //NP others
  int NRaysTot; //NP total number of rays per timestep

  //NP fetch neighborlist data
  ilist = pair_gran->listfull->ilist;
  inum = pair_gran->listfull->inum;
  numneigh = pair_gran->listfull->numneigh;

  //NP fetch particle data
  mask = atom->mask;
  nlocal = atom->nlocal;
  radius = atom->radius;
  type = atom->type;
  x = atom->x;

  //NP TODO should I do this here? What's the purpose?
  updatePtrs();

  // calculate total heat of all particles to update energy of one ray
  if (Qtot == 0.0 || neighbor->decide()){
    for (int ii = 0; ii < inum; ii++){
      i = ilist[ii];
      radi = radius[i];

      areai = MY_4PI * radi * radi;
      emisi = emissivity[type[i]-1];
      tempi = Temp[i];

      Qtot += areai * emisi * Sigma * tempi * tempi * tempi * tempi;
    }
    NRaysTot = avgNRays * inum;
    Qr = Qtot / (double)NRaysTot;
  }

  // all owned particles radiate
  for (int ii = 0; ii < inum; ii++)
  {
    i = ilist[ii];
    //NP TODO
  }

  // reverse communication of heatflux (ghost particles send flux to owner)
  fix_heatFlux->do_reverse_comm();
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::radiate(int i, )


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
  ansD[0] = 2.0*RGen.uniform() - 1.0;
  ansD[1] = 2.0*RGen.uniform() - 1.0;
  ansD[2] = 2.0*RGen.uniform() - 1.0;

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
void FixHeatGranRad::randDir(double *n, double *d){
  double side;

  // generate random direction
  d[0] = 2.0*RGen.uniform() - 1.0;
  d[1] = 2.0*RGen.uniform() - 1.0;
  d[2] = 2.0*RGen.uniform() - 1.0;
  normalize3(d, d);

  // adjust direction, if it points to the wrong side
  // and if n != 0
  side = dot3(d, n);
  if (side < 0.0 && (n[0] != 0.0 || n[1] != 0.0 || n[2] != 0.0)){
    d[0] *= -1.0;
    d[1] *= -1.0;
    d[2] *= -1.0;
  }
}
