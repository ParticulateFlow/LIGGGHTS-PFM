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
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "math_const.h"
#include "math_extra.h"
#include "properties.h"
#include "modify.h"
#include "neigh_list.h"
#include "force.h"
#include "neighbor.h"
#include "pair_gran.h"
#include "random_mars.h"
#include <cstdlib>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

using MathConst::MY_4PI;
using namespace MathExtra;

/* ---------------------------------------------------------------------- */

// assumptions:
// if background_temperature is not passed it is assumed that
//    initial_temperature is the background temperature
FixHeatGranRad::FixHeatGranRad(class LAMMPS *lmp, int narg, char **arg) : FixHeatGran(lmp, narg, arg)
{
  int iarg = 5;

  Qr            = 0.0;
  maxBounces    = 100;
  nTimesteps    = 1000;
  updateCounter = 0;
  cutGhost      = 0.0;

  TB   = -1.0;
  Qtot = 0.0;
  seed = 0;

  fix_emissivity = NULL;
  emissivity     = NULL;

  Sigma = 5.670373E-8;

  // assert that group for 'heat/radiation' is 'all'
  if (strcmp(arg[1],"all")){
    error->fix_error(FLERR, this, "Group for heat/gran/radiation needs to be 'all'.");
  }

  // parse input arguments:
  //  - background_temperature
  //  - max_bounces
  //  - cutoff
  //  - seed
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"background_temperature") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'background_temperature'");
      TB = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"max_bounces") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'max_bounces'");
      maxBounces = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'cutoff'");
      cutGhost = atof(arg[iarg+1]);
      cutGhostsq = cutGhost*cutGhost;
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR, this,"not enough arguments for keyword 'seed'");
      seed = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    }
    else if(strcmp(style,"heat/gran/radiation") == 0)
        error->fix_error(FLERR,this,"unknown keyword");
  }
  if (seed == 0){
    error->fix_error(FLERR,this,"expecting keyword 'seed'");
  }
  if (cutGhost == 0){
    error->fix_error(FLERR, this, "expecting keyword cutoff");
  }

  RGen = new RanMars(lmp, seed + comm->me);

  // for optimization of trace() preallocate these
  raypoint = new double[3];

  stencilLength = NULL;
  binStencildx = NULL;
  binStencilmdx = NULL;
  binStencildy = NULL;
  binStencilmdy = NULL;
  binStencildz = NULL;
  binStencilmdz = NULL;
}

/* ---------------------------------------------------------------------- */

FixHeatGranRad::~FixHeatGranRad()
{
  delete [] emissivity;
  delete [] raypoint;
  delete [] stencilLength;
  delete [] binStencildx;
  delete [] binStencilmdx;
  delete [] binStencildy;
  delete [] binStencilmdy;
  delete [] binStencildz;
  delete [] binStencilmdz;
  delete RGen;
}

/* ---------------------------------------------------------------------- */

int FixHeatGranRad::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

double FixHeatGranRad::extend_cut_ghost() const
{
  return cutGhost;
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::init()
{
  // init base class
  FixHeatGran::init();

  //NP Get pointer to all the fixes (also those that have the material properties)
  updatePtrs();

  // set background temperature to T0 if it was not provided as an argument
  if (TB == -1.0){
    TB = T0;
  }

  // find mechanical-parameter-granular "emissivity" and fetch data
  // if "emissivity" not found: error handling happens inside find_fix_property().
  int max_type = pair_gran->get_properties()->max_type();

  delete [] emissivity;
  emissivity = new double[max_type];

  fix_emissivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalEmissivity","property/global","peratomtype",max_type,0,style));
  for (int i = 0; i < max_type; i++)
  {
    emissivity[i] = fix_emissivity->compute_vector(i);
    if (emissivity[i] < 0.0 || emissivity[i] > 1.0){
      error->all(FLERR,"Fix heat/gran/radiation: Thermal emissivity must not be < 0.0 or > 1.0.");
    }
  }

  // error checks on coarsegraining
  if(force->cg_active())
    error->cg(FLERR,this->style);
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::setup(int i)
{
  // forces algorithm to call updateQr() in post_force()
  Qr = 0.0;

  // create the 2D binStencils needed for ray tracing
  createStencils();
}

/* ---------------------------------------------------------------------- */

void FixHeatGranRad::post_force(int vflag)
{
  //NP neighborlist data
  int i;
  int ibin;

  //NP particle data
  double **x;
  double *radius;
  int *type;
  int nlocal, nghost;

  //NP ray data
  double hitEmis;            //NP emissivity of hit particle
  int hitId;                 //NP index of particle that was hit by ray
  int hitBin;                //NP no of bin where particle has been hit by ray
  double hitp[3]; //NP the point where a ray hit a particle
  double nextNormal[3];

  double *ci;                //NP center of one particle
  double flux;               //NP energy of ray
  double buffer3[3]; //NP buffer for computations in intersectRaySphere
  double d[3]; //NP direction of ray
  double o[3]; //NP origin of ray

  //NP individual particle data
  double areai; //NP area
  double emisi; //NP emissivity
  double radi;  //NP radius
  double tempi; //NP temperature

  //NP fetch particle data
  nlocal = atom->nlocal;
  nghost = atom->nghost;
  radius = atom->radius;
  type   = atom->type;
  x      = atom->x;

  updatePtrs();

  // calculate total heat of all particles to update energy of one ray
  // if (Qr == 0.0 || neighbor->ago() == 0){
  //   updateQr();
  // }

  // all owned particles radiate
  for (i = 0; i < nlocal + nghost; i++)
  {
    // get basic data of atom
    radi  = radius[i];
    ci    = x[i];
    emisi = emissivity[type[i]-1];
    tempi = Temp[i];

    // check in which box we are in
    ibin = neighbor->coord2bin(ci);

    // calculate heat flux from this particle
    areai = MY_4PI * radi * radi;
    flux  = areai * emisi * Sigma * tempi * tempi * tempi * tempi;

    // let this particle radiate (flux is reduced)
    heatFlux[i] -= flux;

    // generate random point and direction
    randOnSphere(ci, radi, o, buffer3);
    randDir(buffer3, d);

    // start radiating and tracing
    // /*NL*/ printf("before trace (post_force)\n");
    hitId = trace(i, ibin, o, d, buffer3, hitp);
    // /*NL*/ printf("after trace (post_force)\n");

    if (hitId != -1){ // the ray hit a particle: reflect from hit particle

      // update heatflux of particle j
      const double& sendflux = flux;
      hitEmis = emissivity[type[hitId]-1];
      heatFlux[hitId] += hitEmis * sendflux;

      hitBin   = neighbor->coord2bin(x[hitId]);

      // calculate new normal vector of reflection for reflected ray
      sub3(hitp, x[hitId], buffer3);
      normalize3(buffer3, nextNormal);

      // reflect ray at the hitpoint.
      reflect(i, hitId, hitBin, hitp, nextNormal, sendflux, 1.0-hitEmis, maxBounces, buffer3);

    } else { // if a ray does not hit a particle we assume radiation from the background
      //NP TODO DEBUG ERROR
      heatFlux[i] += (areai * emisi * Sigma * TB * TB * TB * TB);
    }
  }
}

/* ----------------------------------------------------------------------

// recursive method!
radID ... id of particle that originally radiated (source of energy)
orig_id ... id of particle whereupon the ray was reflected (source of ray)

---------------------------------------------------------------------- */
void FixHeatGranRad::reflect(int radID, int orig_id, int ibin, const double *o, const double *d,
  double flux, double accum_eps, int n, double *buffer3)
{
  const double influx = flux * accum_eps;

  // base case
  if (n == 0){
    heatFlux[radID] += influx;
    return;
  }

  // if energy of one ray would be too small -> stop. //NP TODO: optimize this.
  if (accum_eps < 0.001){
    heatFlux[radID] += influx;
    return;
  }

  // shoot rays.
  double hitp[3];
  int hitId;

  {
    double dd[3];
    // generate random (diffuse) direction
    randDir(d, dd);

    hitId = trace(orig_id, ibin, o, dd, buffer3, hitp);
  }

  int *type = atom->type;

  if (hitId != -1){ // ray hit a particle

    double **x = atom->x;
    const double& sendflux = influx;
    const double hitEmis = emissivity[type[hitId]-1];

    // update heat flux for particle
    heatFlux[hitId] += hitEmis * sendflux;

    // calculate new starting point for reflected ray
    int hitBin = neighbor->coord2bin(x[hitId]);

    double nextNormal[3];

    // calculate new normal vector of reflection for reflected ray
    sub3(hitp, x[hitId], nextNormal);
    norm3(nextNormal);

    // reflect ray at the hitpoint.
    reflect(radID, hitId, hitBin, hitp, nextNormal, flux, (1.0-hitEmis) * accum_eps, n-1, buffer3);
  }
  else {
    const double *radius = atom->radius;
    const double radRad  = radius[radID];
    const double radArea = MY_4PI * radRad * radRad;
    const double radEmis = emissivity[type[radID]-1];
    //NP TODO ERROR DEBUG
    heatFlux[radID] += (radArea * radEmis * accum_eps * Sigma * TB * TB * TB * TB);
  }
}

/* ---------------------------------------------------------------------- */

/*
builds the 2d stencils that are stencils of the borders
of the regular stencil neighbor->stencil when one moves
in directions dx, -dx, dy, -dy, dz or -dz
*/
void FixHeatGranRad::createStencils()
{
  // stencil might need to be re-allocated, if bin-size
  // has been updated
  if(stencilLength == NULL)
    stencilLength = new int[6];

  int mbinx = neighbor->mbinx;
  int mbiny = neighbor->mbiny;
  //int mbinz = neighbor->mbinz;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  double cutneighmaxsq = neighbor->cutneighmaxsq;

  int n;
  int i,j,k;
  int nstencil;

  bool hit;

  // binStencildx
  n = (2*sy+1)*(2*sz+1);
  delete [] binStencildx;
  binStencildx = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (j = -sy; j <= sy; j++){
      hit = false;
      i = sx;
      while (!hit && i >= -sz){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencildx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          i--;
      }
    }
  }
  stencilLength[0] = nstencil;

  // binStencildy
  n = (2*sx+1)*(2*sz+1);
  delete [] binStencildy;
  binStencildy = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      j = sy;
      while (!hit && j >= -sy){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencildy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          j--;
      }
    }
  }
  stencilLength[1] = nstencil;

  // binStencildz
  n = (2*sx+1)*(2*sy+1);
  delete [] binStencildz;
  binStencildz = new int[n];
  nstencil = 0;
  for (j = -sy; j <= sy; j++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      k = sz;
      while (!hit && k >= -sz){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencildz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          k--;
      }
    }
  }
  stencilLength[2] = nstencil;

  // binStencilmdx
  n = (2*sy+1)*(2*sz+1);
  delete [] binStencilmdx;
  binStencilmdx = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (j = -sy; j <= sy; j++){
      hit = false;
      i = -sz;
      while (!hit && i <= sz){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdx[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          i++;
      }
    }
  }
  stencilLength[3] = nstencil;

  // binStencilmdy
  n = (2*sx+1)*(2*sz+1);
  delete [] binStencilmdy;
  binStencilmdy = new int[n];
  nstencil = 0;
  for (k = -sz; k <= sz; k++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      j = -sy;
      while (!hit && j <= sy){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdy[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          j++;
      }
    }
  }
  stencilLength[4] = nstencil;

  // binStencilmdz
  n = (2*sx+1)*(2*sy+1);
  delete [] binStencilmdz;
  binStencilmdz = new int[n];
  nstencil = 0;
  for (j = -sy; j <= sy; j++){
    for (i = -sx; i <= sx; i++){
      hit = false;
      k = -sz;
      while (!hit && k <= sz){
        if (neighbor->bin_distance(i,j,k) < cutneighmaxsq){
          binStencilmdz[nstencil++] = k*mbiny*mbinx + j*mbinx + i;
          hit = true;
        }
        else
          k++;
      }
    }
  }
  stencilLength[5] = nstencil;
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
int FixHeatGranRad::trace(int orig_id, int ibin, const double *o, const double *d, double *buffer3, double *hitp)
{
  int *binhead = neighbor->binhead;
  int *bins = neighbor->bins;
  int stencilbin, stbX, stbY, stbZ;
  int currX, currY, currZ;
  int sx = neighbor->sx;
  int sy = neighbor->sy;
  int sz = neighbor->sz;
  int mbins = neighbor->mbins;

  // individual ray data
  double distsq, distx, disty, distz;
  bool hit;
  double t;

  // global ray data
  bool hitFlag = false;
  double hitT  = 1.0e300;
  int hitId    = -1;

  //NP fetch particle data
  double **x     = atom->x;
  double *radius = atom->radius;

  int i;
  int currentBin = neighbor->coord2bin(o);
  int dx = 0;
  int dy = 0;
  int dz = 0;

  int var_nstencil;
  int *stencil = pair_gran->list->stencil;
  int nstencil = pair_gran->list->nstencil;
  int nstencil2D = 0;

  int *currentStencil = stencil;

  // at first iteration all bins of the full stencil need to be checked
  bool check_boundary_only = false;

  while ((currentBin != -1) && (hitFlag == false)){

    // number of new bins in direction of bin-hop
    var_nstencil = check_boundary_only ? nstencil2D : nstencil;

    // walk the stencil of bins related to this bin, check all of their atoms
    for (int k = 0; k < var_nstencil; k++){

      // check if bin is inside domain of comm->me
      stencilbin = currentBin + currentStencil[k];
      neighbor->bin2XYZ(stencilbin, stbX, stbY, stbZ);
      neighbor->bin2XYZ(currentBin, currX, currY, currZ);
      if (stencilbin < 0 || stencilbin >= mbins || abs(stbX - currX) > sx || abs(stbY - currY) > sy || abs(stbZ - currZ) > sz)
        continue;

      // walk all atoms in this bin
      for (i = binhead[stencilbin]; i >= 0; i = bins[i]){

        // do not intersect with reflecting particle
        if (i == orig_id){
          continue;
        }

        // check if atom intersects ray
        hit = intersectRaySphere(o, d, x[i], radius[i], t, buffer3);
        if (hit){
          hitFlag = true;
          if (t < hitT){
            hitT = t;
            hitId = i;
          }
        }
      }
    }

    if (hitFlag){
      // calculate hit point 'hitp'
      addscaled3(o,d,hitT,hitp);
      return hitId;
    }

    // find the next bin, since in this bin there was no intersection found
    ibin       = currentBin;
    currentBin = nextBin(ibin, o, d, raypoint, dx, dy, dz);

    // from now on only check the boundary any more
    check_boundary_only = true;
    if (dx == 1){
      currentStencil = binStencildx;
      nstencil2D = stencilLength[0];
    }
    else if(dx == -1){
      currentStencil = binStencilmdx;
      nstencil2D = stencilLength[1];
    }
    else if(dy == 1){
      currentStencil = binStencildy;
      nstencil2D = stencilLength[2];
    }
    else if(dy == -1){
      currentStencil = binStencilmdy;
      nstencil2D = stencilLength[3];
    }
    else if(dz == 1){
      currentStencil = binStencildz;
      nstencil2D = stencilLength[4];
    }
    else{
      currentStencil = binStencilmdz;
      nstencil2D = stencilLength[5];
    }

    // maximum radiation distance
    distx = raypoint[0] - o[0];
    disty = raypoint[1] - o[1];
    distz = raypoint[2] - o[2];
    distsq = distx*distx + disty*disty + distz*distz;
    if (distsq >= cutGhostsq){
      return -1;
    }
  }

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
int FixHeatGranRad::nextBin(int ibin, const double *o, const double *d, double *p, int &dx, int &dy, int &dz)
{
  double s;
  double smax = 0.0;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  int nextBinId;

  dx = dy = dz = 0;

  // calculate borders of this bin
  neighbor->binBorders(ibin, xlo, xhi, ylo, yhi, zlo, zhi);

  // xlo, xhi
  if (d[0] != 0.0){
    s = (xlo - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];
    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){
      smax = s;
      dy = dz = 0;
      dx = -1;
    }

    s = (xhi - o[0]) / d[0];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];

    if ((s > smax) && (ylo <= p[1] && p[1] <= yhi) && (zlo <= p[2] && p[2] <= zhi)){

      smax = s;
      dy = dz = 0;
      dx = 1;
    }
  }

  // ylo, yhi
  if (d[1] != 0.0){
    s = (ylo - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){

      smax = s;
      dx = dz = 0;
      dy = -1;
    }

    s = (yhi - o[1]) / d[1];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (zlo <= p[2] && p[2] <= zhi)){

      smax = s;
      dx = dz = 0;
      dy = 1;
    }
  }

  // zlo, zhi
  if (d[2] != 0.0){
    s = (zlo - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){

      smax = s;
      dx = dy = 0;
      dz = -1;
    }

    s = (zhi - o[2]) / d[2];
    p[0] = o[0] + s*d[0];
    p[1] = o[1] + s*d[1];
    p[2] = o[2] + s*d[2];

    if ((s > smax) && (xlo <= p[0] && p[0] <= xhi) && (ylo <= p[1] && p[1] <= yhi)){

      smax = s;
      dx = dy = 0;
      dz = 1;
    }
  }

  p[0] = o[0] + smax*d[0];
  p[1] = o[1] + smax*d[1];
  p[2] = o[2] + smax*d[2];

  nextBinId = neighbor->binHop(ibin,dx,dy,dz);

  if (nextBinId != ibin){
    return nextBinId;
  }
  else {
    error->one(FLERR, "FixHeatGranRad::nextBin() did not find a suitable next bin.");
    return -1;
  }

}

/* ---------------------------------------------------------------------- */
/* calculates energy per ray from total radiative energy in the system */

void FixHeatGranRad::updateQr()
{
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  double areai, emisi, tempi;
  double radi;
  unsigned long NRaysTot;

  Qtot = 0.0;
  for (int i = 0; i < nlocal + nghost; i++){
    radi = radius[i];

    areai = MY_4PI * radi * radi;
    emisi = emissivity[type[i]-1];
    tempi = Temp[i];

    Qtot += areai * emisi * Sigma * tempi * tempi * tempi * tempi;
  }

  NRaysTot = nlocal + nghost;
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
bool FixHeatGranRad::intersectRaySphere(const double *o, const double *d, const double *center, double radius, double &t, double *buffer3)
{
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

  // since the variale "raypoint" in nextBin() could possibly be inside a sphere
  // we leave this part.
  // // if t0 is less than zero, the intersection point is at t1
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

 * This is a variant of the algorithm for computing a random point
 * on the unit sphere; the algorithm is suggested in Knuth, v2,
 * 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
 * 326.
 * see http://fossies.org/dox/gsl-2.6/sphere_8c_source.html#l00066
*/
void FixHeatGranRad::randOnSphere(const double *c, double r, double *ansP, double *ansD)
{
  double s;

  // generate random direction
  do
  {
    ansD[0] = 2.0*RGen->uniform() - 1.0;
    ansD[1] = 2.0*RGen->uniform() - 1.0;
    s = ansD[0]*ansD[0] + ansD[1]*ansD[1];

  } while (s > 1.0);

  ansD[2] = 2.0*s - 1.0;
  const double a = 2.0 * sqrt(1.0 - s);

  ansD[0] *= a;
  ansD[1] *= a;

  // calculate corresponding point on surface
  addscaled3(c,ansD,r,ansP);

}

/* ---------------------------------------------------------------------- */
/*
  generates a random direction in direction of the vector n
  n ... input - normal vector
  o ... output - random direction

* This is a variant of the algorithm for computing a random point
* on the unit sphere; the algorithm is suggested in Knuth, v2,
* 3rd ed, p136; and attributed to Robert E Knop, CACM, 13 (1970),
* 326.
* see http://fossies.org/dox/gsl-2.6/sphere_8c_source.html#l00066
*/
void FixHeatGranRad::randDir(const double *n, double *d)
{
  double s;

  do
  {
    d[0] = 2.0*RGen->uniform() - 1.0;
    d[1] = 2.0*RGen->uniform() - 1.0;
    s = d[0]*d[0] + d[1]*d[1];

  } while (s > 1.0);

  d[2] = 2.0*s - 1.0;
  const double a = 2.0 * sqrt(1.0 - s);

  d[0] *= a;
  d[1] *= a;

  // adjust direction, if it points to the wrong side
  const double side = dot3(d, n);
  if (side < 0.0){
    d[0] = -d[0];
    d[1] = -d[1];
    d[2] = -d[2];
  }
}
