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
#include "fix_property_atom.h"
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
//DEBUG
#include "update.h"
#include <iostream>
#include <fstream>
//ENDDEBUG

using namespace LAMMPS_NS;
using namespace FixConst;

using MathExtra::add3;
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

  //DEBUG
  bool debugFlg = false;
  //ENDDEBUG

	//NP use initial temperature as background temperature if
	//NP TB is not passed as argument background_temperature
  //NP background temperature is further handled in init()

  Qr         = 0.0;
  avgNRays   = 1; //NP TODO optimize these parameters
  maxBounces = 3;

  TB   = -1.0;
  Qtot = 0.0;

  fix_emissivity = NULL;
  emissivity     = NULL;

  Sigma = 5.670373E-8;

  // assert that group for 'heat/radiation' is 'all'
  if (strcmp(arg[1],"all")){
    error->fix_error(FLERR, this, "Group for heat/gran/radiation needs to be 'all'.");
  }

  // parse input arguments:
  //  - backgroundTemperature
  //  - averageNumberOfRaysPerParticle
  //  - maxNumberOfBounces
	bool hasargs = true;
	while(iarg < narg && hasargs)
  {
    //DEBUG
    if (debugFlg) printf("DEBUG: SETTING PARAMETERS\n");
    if (debugFlg) printf("DEBUG: arg[%d]: %s\n", iarg, arg[iarg]);
    //ENDDEBUG
    hasargs = false;
    if(strcmp(arg[iarg],"backgroundTemperature") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'backgroundTemperature'");
      TB = atof(arg[iarg+1]);
      //DEBUG
      if (debugFlg) printf("DEBUG: TB: %f\n", TB);
      //ENDDEBUG
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"averageNumberOfRaysPerParticle") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'averageNumberOfRaysPerParticle'");
      avgNRays = atoi(arg[iarg+1]);
      //DEBUG
      if (debugFlg) printf("DEBUG: avgNRays: %d\n", avgNRays);
      //ENDDEBUG
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"maxNumberOfBounces") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'maxNumberOfBounces'");
      maxBounces = atoi(arg[iarg+1]);
      //DEBUG
      if (debugFlg) printf("DEBUG: maxBounces: %d\n", maxBounces);
      //ENDDEBUG
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
    if (emissivity[i] < 0.0 || emissivity[i] > 1.0){
      error->all(FLERR,"Fix heat/gran/radiation: Thermal emissivity must not be < 0.0 or > 1.0.");
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::post_force(int vflag){

  //DEBUG
  bool debugFlg = false;
  if (debugFlg) printf("DEBUG: in post_force()\n");
  if (debugFlg) printf("DEBUG: TB: %f\n",TB);
  //ENDDEBUG
  //NP neighborlist data
  int i, j;
  int ibin;

  //NP particle data
  double **x;
  double *radius;
  int *mask;
  int *type;
  int nlocal, nghost;

  //NP ray data
  bool hit;                  //NP if a ray hit one particular particle
  double *hitp = new double[3]; //NP the point where a ray hit a particle
  double hitEmis;            //NP emissivity of hit particle
  int hitId;                 //NP index of particle that was hit by ray
  int hitBin;                //NP no of bin where particle has been hit by ray
  double *nextO = new double[3];
  double *nextNormal = new double[3];

  double *buffer3 = new double[3]; //NP buffer for computations in intersectRaySphere
  double *ci;                //NP center of one particle
  double *cj;                //NP center of another particle
  double *d = new double[3]; //NP direction of ray
  double *o = new double[3]; //NP origin of ray
  double flux;               //NP energy of ray
  double sendflux;           //NP flux that is sent in a particular ray
  double t;                  //NP parameter of ray
  int m;                     //NP number of rays for one specific particle
  int NRaysTot;              //NP total number of rays per timestep

  //NP individual particle data
  double areai; //NP area
  double areaj; //NP area
  double emisi; //NP emissivity
  double emisj; //NP emissivity
  double radi;  //NP radius
  double radj;  //NP radius
  double tempi; //NP temperature

  //NP fetch particle data
  mask   = atom->mask;
  nlocal = atom->nlocal;
  nghost = atom->nghost;
  radius = atom->radius;
  type   = atom->type;
  x      = atom->x;

  //NP TODO should I do this here? What's the purpose?
  updatePtrs();

  // calculate total heat of all particles to update energy of one ray
  if (Qtot == 0.0 || neighbor->decide()){
    updateQr();

  //DEBUG
  if (debugFlg) printf("DEBUG: updated Qr\n");
  if (debugFlg) printf("DEBUG: Qr:%f\n", Qr);
  //ENDDEBUG
  }

  // all owned particles radiate
  for (i = 0; i < nlocal + nghost; i++)
  {
    // get basic data of atom
    radi  = radius[i];
    ci    = x[i];
    emisi = emissivity[type[i]-1];
    tempi = Temp[i];
    //DEBUG
    if (debugFlg) printf("------------------------------------------------------\n");
    if (debugFlg) printf("DEBUG: radiating from atom %d with center: [%f %f %f] and radius %f\n", i, ci[0], ci[1], ci[2], radi);
    //ENDDEBUG

    // check in which box we are in
    ibin = neighbor->coord2bin(ci);
    //DEBUG
    // if (debugFlg) printf("DEBUG: ibin: %d\n", ibin);
    //ENDDEBUG

    // calculate heat flux from this particle
    areai = MY_4PI * radi * radi;
    flux  = areai * emisi * Sigma * tempi * tempi * tempi * tempi;
    //DEBUG
    // if (debugFlg) printf("DEBUG: areai: %f\n", areai);
    // if (debugFlg) printf("DEBUG: emisi: %f\n", emisi);
    // if (debugFlg) printf("DEBUG: Sigma: %.10f\n", Sigma);
    if (debugFlg) printf("DEBUG: tempi: %f\n", tempi);
    if (debugFlg) printf("DEBUG: heatFlux[%d] -= %f\n", i, flux);
    //ENDDEBUG

    // calculate number of rays
    m = (int)(flux/Qr) + 1;

    // calculate effective energy of each ray
    sendflux = flux / (double) m;

    // let this particle radiate (flux is reduced)
    heatFlux[i] -= flux;

    // shoot rays
    for (int r = 0; r < m; r++){
      //DEBUG
      if (debugFlg) printf("DEBUG: handling ray %d of %d\n", r+1, m);
      //ENDDEBUG

      // generate random point and direction
      randOnSphere(ci, radi, o, buffer3);
      randDir(buffer3, d);
      //DEBUG
      // if (debugFlg) printf("DEBUG: o: [%f %f %f]\n", o[0], o[1], o[2]);
      // if (debugFlg) printf("DEBUG: d: [%f %f %f]\n", d[0], d[1], d[2]);
      //ENDDEBUG

      // start radiating and tracing
      hitId = trace(i, ibin, o, d, sendflux, buffer3, hitp);

      if (hitId != -1){ // the ray hit a particle: reflect from hit particle

        // update heatflux of particle j
        hitEmis = emissivity[type[hitId]-1];
        heatFlux[hitId] += hitEmis * sendflux;
        //DEBUG
        if (debugFlg) printf("DEBUG: a hit -> heatFlux[%d] += %f\n", hitId, hitEmis*sendflux);
        //ENDDEBUG

        hitBin   = neighbor->coord2bin(x[hitId]);
        nextO[0] = hitp[0];
        nextO[1] = hitp[1];
        nextO[2] = hitp[2];

        // calculate new normal vector of reflection for reflected ray
        sub3(hitp, x[hitId], buffer3);
        normalize3(buffer3, nextNormal);

        // reflect ray at the hitpoint.
        reflect(hitId, hitBin, nextO, nextNormal, (1.0-hitEmis) * sendflux, maxBounces, buffer3);

      } else { // if a ray does not hit a particle we assume radiation from the background
        heatFlux[i] += (areai * emisi * Sigma * TB * TB * TB * TB) / ((double) m);
        //DEBUG
        if (debugFlg) printf("DEBUG: nohit -> heatFlux[%d] += %f\n", i, (areai * emisi * Sigma * TB * TB * TB * TB) / ((double) m));
        // if (debugFlg) error->all(FLERR, "muhaha");
        //ENDDEBUG
      }
    }
  }

  //DEBUG
  if (debugFlg) printf("\nDEBUG: === summary of FLUX ===\n");
  for (int i = 0; i < nlocal; i++){
    if (debugFlg) printf("DEBUG: flux[%d]: %f\n", i, heatFlux[i]);
  }
  if (debugFlg) printf("\n");
  //ENDDEBUG

  // reverse communication of heatflux (ghost particles send flux to owner)
  // fix_heatFlux->do_reverse_comm();
  //DEBUG
  // if (debugFlg) error->all(FLERR, "end this, now!!");
  //ENDDEBUG
  delete [] buffer3;
  delete [] d;
  delete [] hitp;
  delete [] nextNormal;
  delete [] nextO;
  delete [] o;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::reflect(int orig_id, int ibin, double *o, double *d,
  double influx, int n, double *buffer3){

  //DEBUG
  bool debugFlg = false;
  if (debugFlg) printf("DEBUG: reflect()\n");
  if (debugFlg) printf("DEBUG: n: %d\n", n);
  //ENDDEBUG

  double ratio;
  double sendflux;
  double hitEmis;
  int hitId, hitBin;
  int m;

  // base case
  if (n == 0){
    return;
  }

  ratio = influx / Qr;
  m = (int) ratio + 1;

  //DEBUG
  if (debugFlg) printf("DEBUG: ratio: %f\n", ratio);
  if (debugFlg) printf("DEBUG: m: %d\n", m);
  //ENDDEBUG

  // if energy of one ray would be too small -> stop.
  //NP TODO: optimize this parameter
  if (ratio < 0.01){
    return;
  }

  sendflux = influx / (double) m;

  // shoot rays.
  double **x         = atom->x;
  double *dd         = new double[3];
  double *hitp       = new double[3];
  double *nextNormal = new double[3];
  double *nextO      = new double[3];
  int *type          = atom->type;

  for (int r = 0; r < m; r++){

    // generate random (diffuse) direction
    randDir(d, dd);

    hitId = trace(orig_id, ibin, o, dd, sendflux, buffer3, hitp);

    if (hitId != -1){ // ray hit a particle

      // update heat flux for particle
      hitEmis = emissivity[type[hitId]-1];
      heatFlux[hitId] += hitEmis * sendflux;
      //DEBUG
      if (debugFlg) printf("DEBUG: hitId: %d\n", hitId);
      if (debugFlg) printf("DEBUG: hitEmis: %f\n", hitEmis);
      if (debugFlg) printf("DEBUG: sendflux: %f\n", sendflux);
      if (debugFlg) printf("DEBUG: a hit -> heatFlux[%d] += %f\n", hitId, hitEmis*sendflux);
      //ENDDEBUG

      // calculate new starting point for reflected ray
      hitBin = neighbor->coord2bin(x[hitId]);
      nextO[0] = hitp[0];
      nextO[1] = hitp[1];
      nextO[2] = hitp[2];

      // calculate new normal vector of reflection for reflected ray
      sub3(hitp, x[hitId], buffer3);
      normalize3(buffer3, nextNormal);

      // reflect ray at the hitpoint.
      reflect(hitId, hitBin, nextO, nextNormal, (1.0-hitEmis) * sendflux, n-1, buffer3);
    }
  }

  delete [] dd;
  delete [] hitp;
  delete [] nextNormal;
  delete [] nextO;

}

/* ---------------------------------------------------------------------- */
/*
  orig_id .. input - id of particle that radiated
  ibin .. input - bin within which tracing starts
  o ... input - origin of ray
  d ... input - direction of ray (length 1)
  flux ... input - heatflux, energy of rays.
  n ... input - recursion parameter (if 0 recursion stops)
  buffer3 ... neutral - non-const, no meaning
  hitp ... output - point where ray hit a particle

  return value:
  particle id that has been hit, or -1
*/
int FixHeatGranRad::trace(int orig_id, int ibin, double *o, double *d,
  double flux, double *buffer3, double *hitp){
  //DEBUG
  bool debugFlg = false;
  //ENDDEBUG

  //DEBUG
  if (debugFlg) printf("DEBUG: in trace()\n");
  // printf("DEBUG: radiating atom: %d\n", orig_id);
  //ENDDEBUG

  double *raypoint = new double[3];

  int *binhead = neighbor->binhead;
  int *bins = neighbor->bins;

  // individual particle data
  double *c;  // center of particle i
  double rad; // radius of particle i

  // individual ray data
  bool hit;
  double t;

  // global ray data
  bool hitFlag = false;
  double hitT  = 1.0e300;
  int hitId    = -1;

  //NP fetch particle data
  double **x     = atom->x;
  double *radius = atom->radius;
  int *mask      = atom->mask; //NP? needed ?
  int *type      = atom->type; //NP? needed ?

  int i, j;
  int currentBin = neighbor->coord2bin(o);
  int dx, dy, dz;
  int *binStencil = new int[27]; // binStencil[x*9 + y*3 + z]

  // initialize binStencil
  for (i = 0; i < 27; i++)
    binStencil[i] = -1;
  j = 0;
  for (int ix = -1; ix <= 1; ix++){
    for (int iy = -1; iy <= 1; iy++){
      for (int iz = -1; iz <= 1; iz++){
        binStencil[j++] = neighbor->binHop(currentBin, ix, iy, iz);
      }
    }
  }

  while ((currentBin != -1) && (hitFlag == false)){

    //DEBUG
    // printf("DEBUG: currentBin: %d\n", currentBin);
    // for (int debug000 = 0; debug000 < 27; debug000++){
      // printf("DEBUG: binStencil[%d]: %d\n", debug000, binStencil[debug000]);
    // }
    //ENDDEBUG

    // walk the stencil of bins related to this bin, check all of their atoms
    for (int k = 0; k < 27; k++){
      if (binStencil[k] == -1){
        continue;
      }
      //DEBUG
      if (debugFlg) printf("DEBUG: checking bin (from stencil) %d\n", binStencil[k]);
      //ENDDEBUG
      // first atom in this bin
      i = binhead[binStencil[k]];

      // walk all atoms in this bin
      while (i != -1){
        //DEBUG
        if (debugFlg) printf("DEBUG: checking atom: %d\n", i);
        //ENDDEBUG
        // do not intersect with reflecting particle
        if (i == orig_id){
          i = bins[i];
          continue;
        }

        // center and radius
        c = x[i];
        rad = radius[i];

        // check if atom intersects ray
        hit = intersectRaySphere(o, d, c, rad, t, buffer3);
        if (hit){
          hitFlag = true;
          if (t < hitT){
            hitT = t;
            hitId = i;
          }
        }

        // next atom
        i = bins[i];
      }
    }

    if (hitFlag){
      // calculate hit point 'hitp'
      snormalize3(hitT, d, buffer3);
      add3(o, buffer3, hitp);
      //DEBUG
      if (debugFlg) printf("DEBUG: hit in stencil around bin %d\n", currentBin);
      //ENDDEBUG
      delete [] raypoint;
      delete [] binStencil;
      return hitId;
    }
    //DEBUG
    if (debugFlg) printf("DEBUG: no hit in stencil around bin %d\n", currentBin);
    //ENDDEBUG
    // find the next bin, since in this bin there was no intersection found
    ibin       = currentBin;
    currentBin = nextBin(ibin, o, d, raypoint, dx, dy, dz);
    o[0] = raypoint[0];
    o[1] = raypoint[1];
    o[2] = raypoint[2];

    // update binStencil, only bins that are in the direction currentBin was
    // moved to are set to non-negative, old ones are not checked again.
    j = 0;
    for (int ix = -1; ix <= 1; ix++){
      for (int iy = -1; iy <= 1; iy++){
        for (int iz = -1; iz <= 1; iz++){
          if (
            (dx != 0 && ix == dx) ||
            (dy != 0 && iy == dy) ||
            (dz != 0 && iz == dz)
            ){ // set the stencil to non-negative only for new bins
            binStencil[j] = neighbor->binHop(ibin, ix, iy, iz);
          }
          else
            binStencil[j] = -1;
          j++;
        }
      }
    }
  }

  delete [] raypoint;
  delete [] binStencil;
  return -1;
}

/* ---------------------------------------------------------------------- */
/*
returns -1 if border of processor is hit
returns id of new bin otherwise

- input:
ibin ... input - local index of current bin
o ... origin of ray
d ... direction of ray

- output:
p ... intersection point of ray and border of bin
dx ... bin hopped in x direction? (-1/0/1)
dy ... bin hopped in y direction? (-1/0/1)
dz ... bin hopped in z direction? (-1/0/1)
*/
int FixHeatGranRad::nextBin(int ibin, double *o, double *d, double *p, int &dx, int &dy, int &dz){
  //DEBUG
  // printf("DEBUG: in nextBin()\n");
  //ENDDEBUG
  double s;
  double smax = 0.0;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  int nextBinId;
  //DEBUG
  bool debugFlg = false;
  //ENDDEBUG

  dx = dy = dz = 0;

  //DEBUG
  // d[0] = 1.0;
  // d[1] = 0.0;
  // d[2] = 0.0;
  //ENDDEBUG

  // calculate borders of this bin
  neighbor->binBorders(ibin, xlo, xhi, ylo, yhi, zlo, zhi);

  //DEBUG
  if (debugFlg) printf("DEBUG: bin: %d\n", ibin);
  if (debugFlg) printf("DEBUG: binBorders: xlo: %f, xhi: %f, ylo: %f, yhi: %f, zlo: %f, zhi: %f\n", xlo, xhi, ylo, yhi, zlo, zhi);
  if (debugFlg) printf("DEBUG: o: %f %f %f\n",o[0],o[1],o[2]);
  if (debugFlg) printf("DEBUG: d: %f %f %f\n",d[0],d[1],d[2]);
  //ENDDEBUG

  // xlo, xhi
  if (d[0] != 0.0){
    //DEBUG
    if (debugFlg) printf("\nDEBUG: case xlo, xhi\n");
    //ENDDEBUG
    s = (xlo - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG
      smax = s;
      dy = dz = 0;
      dx = -1;
    }

    s = (xhi - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG

    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG

      smax = s;
      dy = dz = 0;
      dx = 1;
    }
  }

  // ylo, yhi
  if (d[1] != 0.0){
    //DEBUG
    if (debugFlg) printf("\nDEBUG: case ylo, yhi\n");
    //ENDDEBUG
    s = (ylo - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG


    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG

      smax = s;
      dx = dz = 0;
      dy = -1;
    }

    s = (yhi - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG

      smax = s;
      dx = dz = 0;
      dy = 1;
    }
  }

  // zlo, zhi
  if (d[2] != 0.0){
    //DEBUG
    if (debugFlg) printf("\nDEBUG: case zlo, zhi\n");
    //ENDDEBUG
    s = (zlo - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG

      smax = s;
      dx = dy = 0;
      dz = -1;
    }

    s = (zhi - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: s: %f\n", s);
    if (debugFlg) printf("DEBUG: p: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){
      //DEBUG
      if (debugFlg) printf("DEBUG: inside conditional\n");
      //ENDDEBUG

      smax = s;
      dx = dy = 0;
      dz = 1;
    }
  }

  p[0] = o[0] + smax*d[0];
  p[1] = o[1] + smax*d[1];
  p[2] = o[2] + smax*d[2];
    //DEBUG
    if (debugFlg) printf("DEBUG: sfin: %f\n", smax);
    if (debugFlg) printf("DEBUG: pfin: [%f %f %f]\n", p[0], p[1], p[2]);
    //ENDDEBUG

  nextBinId = neighbor->binHop(ibin,dx,dy,dz);
  if (nextBinId != ibin){
    //DEBUG
    if (debugFlg) printf("DEBUG: returning next bin: %d\n", nextBinId);
    //ENDDEBUG
    return nextBinId;
  }
  else {
    //DEBUG
    if (debugFlg) printf("DEBUG: WARNING: nextBin() did not find suitable next bin.\n");
    //ENDDEBUG
    error->all(FLERR, "FixHeatGranRad::nextBin() did not find a suitable next bin.");
    return -1;
  }

}

/* ---------------------------------------------------------------------- */
/* calculates energy per ray from total radiative energy in the system */

void FixHeatGranRad::updateQr(){

  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  double areai, emisi, tempi;
  double Qtot = 0.0;
  double radi;
  int i;
  unsigned long NRaysTot;

  for (int i = 0; i < nlocal + nghost; i++){
    radi = radius[i];

    areai = MY_4PI * radi * radi;
    emisi = emissivity[type[i]-1];
    tempi = Temp[i];

    Qtot += areai * emisi * Sigma * tempi * tempi * tempi * tempi;
  }

  NRaysTot = avgNRays * nlocal;
  Qr = Qtot / (double)NRaysTot;
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

note:
if a ray originates from within a sphere its intersection point is the point on
the sphere, where the ray's parameter is negative!
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
  // if (t0 < 0.0){
  //   t = t1;
  // } else {
  //   t = t0;
  // }

  t = t0;
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
  n ... input - normal vector
  o ... output - random direction
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
