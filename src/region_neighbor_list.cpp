#include "mpi.h"
#include "region_neighbor_list.h"

using namespace LAMMPS_NS;
using namespace std;

RegionNeighborList::RegionNeighborList() {
}


bool RegionNeighborList::hasOverlap(double * x, double radius) const {

  for(vector<Particle>::const_iterator it = xnear.begin(); it != xnear.end(); ++it) {
    const Particle & p = *it;
    double del[3];
    vectorSubtract3D(x, p.x, del);
    const double rsq = vectorMag3DSquared(del);
    const double radsum = radius + p.radius;

    if (rsq <= radsum*radsum) return true;
  }

  return false;
}


void RegionNeighborList::insert(double * x, double radius) {
  xnear.push_back(Particle(x, radius));
}

void RegionNeighborList::reset() {
  xnear.clear();
}

int RegionNeighborList::count() const {
  return xnear.size();
}
