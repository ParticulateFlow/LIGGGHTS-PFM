#include "lmptype.h"
#include "mpi.h"
#include "region_neighbor_list.h"
#include <limits>
#include "assert.h"

#define SMALL 1.0e-6
#define BIG 1.0e20
#define CUT2BIN_RATIO 100

using namespace LAMMPS_NS;
using namespace std;

RegionNeighborList::RegionNeighborList() {
}


bool RegionNeighborList::hasOverlap(double * x, double radius) const {

  /*for(vector<Particle>::const_iterator it = xnear.begin(); it != xnear.end(); ++it) {
    const Particle & p = *it;
    double del[3];
    vectorSubtract3D(x, p.x, del);
    const double rsq = vectorMag3DSquared(del);
    const double radsum = radius + p.radius;

    if (rsq <= radsum*radsum) return true;
  }*/

  int ibin = coord2bin(x);

  for(std::vector<int>::const_iterator it = stencil.begin(); it != stencil.end(); ++it) {
    const int offset = *it;
    assert(ibin+offset >= 0 || ibin+offset < bins.size());
    const ParticleBin & bin = bins[ibin+offset];

    for(ParticleBin::const_iterator pit = bin.begin(); pit != bin.end(); ++pit) {
      const Particle & p = *pit;
      double del[3];
      vectorSubtract3D(x, p.x, del);
      const double rsq = vectorMag3DSquared(del);
      const double radsum = radius + p.radius;
      if (rsq <= radsum*radsum) return true;
    }
  }

  return false;
}


void RegionNeighborList::insert(double * x, double radius) {
  //xnear.push_back(Particle(x, radius));

  int ibin = coord2bin(x);
  bins[ibin].push_back(Particle(x, radius));
  ++ncount;
}

void RegionNeighborList::reset() {
  xnear.clear();
  bins.clear();
  stencil.clear();
  ncount = 0;
}

size_t RegionNeighborList::count() const {
  //return xnear.size();
  return ncount;
}

bool RegionNeighborList::setBoundingBox(InsertBoundingBox & bb, double maxrad) {

  printf("setting insertion bounding box: [%g, %g] x [%g, %g] x [%g, %g]\n", bb.xlo, bb.xhi, bb.ylo, bb.yhi, bb.zlo, bb.zhi);


  double bbox[3];

  bbox[0] = bb.xhi - bb.xlo;
  bbox[1] = bb.yhi - bb.ylo;
  bbox[2] = bb.zhi - bb.zlo;

  if(bbox[0]*bbox[1]*bbox[2] < 0) {
    // insertion area not in subbox
    bins.clear();
    stencil.clear();
    return false;
  }

  bboxlo[0] = bb.xlo;
  bboxlo[1] = bb.ylo;
  bboxlo[2] = bb.zlo;

  bboxhi[0] = bb.xhi;
  bboxhi[1] = bb.yhi;
  bboxhi[2] = bb.zhi;


  // testing code
  double binsize_optimal = 4*maxrad;
  double binsizeinv = 1.0/binsize_optimal;

  // test for too many global bins in any dimension due to huge global domain
  const int max_small_int = std::numeric_limits<int>::max();

  if (bbox[0]*binsizeinv > max_small_int || bbox[1]*binsizeinv > max_small_int ||
      bbox[2]*binsizeinv > max_small_int) {
    printf("ERROR\n");
    return false;
  }

  // create actual bins
  nbinx = static_cast<int>(bbox[0]*binsizeinv);
  nbiny = static_cast<int>(bbox[1]*binsizeinv);
  nbinz = static_cast<int>(bbox[2]*binsizeinv);

  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;

  printf("nbinx %d nbiny %d, nbinz %d\n",nbinx,nbiny,nbinz);

  binsizex = bbox[0]/nbinx;
  binsizey = bbox[1]/nbiny;
  binsizez = bbox[2]/nbinz;

  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  /*if (binsize_optimal*bininvx > CUT2BIN_RATIO ||
      binsize_optimal*bininvy > CUT2BIN_RATIO ||
      binsize_optimal*bininvz > CUT2BIN_RATIO)
    error->all(FLERR,"Cannot use neighbor bins - box size << cutoff");*/

  // mbinlo/hi = lowest and highest global bins my ghost atoms could be in
  // coord = lowest and highest values of coords for my ghost atoms
  // static_cast(-1.5) = -1, so subract additional -1
  // add in SMALL for round-off safety

  double bsubboxlo[3], bsubboxhi[3];

  bsubboxlo[0] = bb.xlo;
  bsubboxlo[1] = bb.ylo;
  bsubboxlo[2] = bb.zlo;
  bsubboxhi[0] = bb.xhi;
  bsubboxhi[1] = bb.yhi;
  bsubboxhi[2] = bb.zhi;

  double coord = bsubboxlo[0] - SMALL*bbox[0];
  mbinxlo = static_cast<int> ((coord-bboxlo[0])*bininvx);
  if (coord < bboxlo[0]) mbinxlo = mbinxlo - 1;
  coord = bsubboxhi[0] + SMALL*bbox[0];
  int mbinxhi = static_cast<int> ((coord-bboxlo[0])*bininvx);

  coord = bsubboxlo[1] - SMALL*bbox[1];
  mbinylo = static_cast<int> ((coord-bboxlo[1])*bininvy);
  if (coord < bboxlo[1]) mbinylo = mbinylo - 1;
  coord = bsubboxhi[1] + SMALL*bbox[1];
  int mbinyhi = static_cast<int> ((coord-bboxlo[1])*bininvy);

  coord = bsubboxlo[2] - SMALL*bbox[2];
  mbinzlo = static_cast<int> ((coord-bboxlo[2])*bininvz);
  if (coord < bboxlo[2]) mbinzlo = mbinzlo - 1;
  coord = bsubboxhi[2] + SMALL*bbox[2];
  int mbinzhi = static_cast<int> ((coord-bboxlo[2])*bininvz);

  // extend bins by 1 to insure stencil extent is included

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;

  printf("mbinxlo: %d, mbinxhi: %d\n", mbinxlo, mbinxhi);

  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  printf("mbinylo: %d, mbinyhi: %d\n", mbinylo, mbinyhi);
  mbiny = mbinyhi - mbinylo + 1;

  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  printf("mbinzlo: %d, mbinzhi: %d\n", mbinzlo, mbinzhi);
  mbinz = mbinzhi - mbinzlo + 1;

  printf("mbinx %d mbiny %d, mbinz %d\n",mbinx,mbiny,mbinz);

  double x0[] = {0.001, 0.001, 0.001 };
  int ibinzero = coord2bin(x0);
  printf("ibinzero %d\n", ibinzero);

  // memory for bin ptrs

  bigint bbin = ((bigint) mbinx) * ((bigint) mbiny) * ((bigint) mbinz);
  if (bbin > max_small_int) {
    printf("Too many neighbor bins\n");
    return false;
  }
  bins.resize(bbin);


  for (int k = -1; k <= 1; k++)
    for (int j = -1; j <= 1; j++)
      for (int i = -1; i <= 1; i++)
        stencil.push_back(k*mbiny*mbinx + j*mbinx + i);

  return true;
}

/* ----------------------------------------------------------------------
   compute closest distance between central bin (0,0,0) and bin (i,j,k)
------------------------------------------------------------------------- */

double RegionNeighborList::bin_distance(int i, int j, int k)
{
  double delx,dely,delz;

  if (i > 0) delx = (i-1)*binsizex;
  else if (i == 0) delx = 0.0;
  else delx = (i+1)*binsizex;

  if (j > 0) dely = (j-1)*binsizey;
  else if (j == 0) dely = 0.0;
  else dely = (j+1)*binsizey;

  if (k > 0) delz = (k-1)*binsizez;
  else if (k == 0) delz = 0.0;
  else delz = (k+1)*binsizez;

  return (delx*delx + dely*dely + delz*delz);
}


int RegionNeighborList::coord2bin(double *x) const
{
  int ix,iy,iz;

  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = std::min(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = std::min(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = std::min(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/* ----------------------------------------------------------------------
   same as coord2bin, but also return ix,iy,iz offsets in each dim
------------------------------------------------------------------------- */

int RegionNeighborList::coord2bin(double *x, int &ix, int &iy, int &iz) const
{
  if (x[0] >= bboxhi[0])
    ix = static_cast<int> ((x[0]-bboxhi[0])*bininvx) + nbinx;
  else if (x[0] >= bboxlo[0]) {
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx);
    ix = std::min(ix,nbinx-1);
  } else
    ix = static_cast<int> ((x[0]-bboxlo[0])*bininvx) - 1;

  if (x[1] >= bboxhi[1])
    iy = static_cast<int> ((x[1]-bboxhi[1])*bininvy) + nbiny;
  else if (x[1] >= bboxlo[1]) {
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy);
    iy = std::min(iy,nbiny-1);
  } else
    iy = static_cast<int> ((x[1]-bboxlo[1])*bininvy) - 1;

  if (x[2] >= bboxhi[2])
    iz = static_cast<int> ((x[2]-bboxhi[2])*bininvz) + nbinz;
  else if (x[2] >= bboxlo[2]) {
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz);
    iz = std::min(iz,nbinz-1);
  } else
    iz = static_cast<int> ((x[2]-bboxlo[2])*bininvz) - 1;

  ix -= mbinxlo;
  iy -= mbinylo;
  iz -= mbinzlo;
  return iz*mbiny*mbinx + iy*mbinx + ix;
}
