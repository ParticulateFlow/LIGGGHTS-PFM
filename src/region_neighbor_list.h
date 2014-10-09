#ifndef REGION_NEIGHBOR_LIST_H
#define REGION_NEIGHBOR_LIST_H

#include <vector>
#include "vector_liggghts.h"

struct Particle {
  double x[3];
  double radius;

  Particle(double * pos, double rad) {
    LAMMPS_NS::vectorCopy3D(pos, x);
    radius = rad;
  }
};

class RegionNeighborList
{
  std::vector<Particle> xnear;

public:
    RegionNeighborList();

    bool hasOverlap(double * x, double radius) const;
    void insert(double * x, double radius);
    int count() const;
    void reset();
};

#endif // REGION_NEIGHBOR_LIST_H
