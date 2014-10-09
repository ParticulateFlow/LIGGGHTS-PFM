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

struct InsertBoundingBox {
  double xlo;
  double xhi;
  double ylo;
  double yhi;
  double zlo;
  double zhi;

  InsertBoundingBox() :
    xlo(0.0),
    xhi(0.0),
    ylo(0.0),
    yhi(0.0),
    zlo(0.0),
    zhi(0.0)
  {}
};


class RegionNeighborList
{
  typedef std::vector<Particle> ParticleBin;

  std::vector<Particle> xnear;
  std::vector<ParticleBin> bins;
  std::vector<int> stencil;
  size_t ncount;

  double bboxlo[3];
  double bboxhi[3];

  int nbinx,nbiny,nbinz;           // # of global bins
  int mbinx,mbiny,mbinz;
  int mbinxlo,mbinylo,mbinzlo;

  double binsizex,binsizey,binsizez;  // actual bin sizes and inverse sizes
  double bininvx,bininvy,bininvz;

  double bin_distance(int i, int j, int k);
  int coord2bin(double *x) const;
  int coord2bin(double *x, int &ix, int &iy, int &iz) const;

public:
    RegionNeighborList();

    bool hasOverlap(double * x, double radius) const;
    void insert(double * x, double radius);
    size_t count() const;
    void reset();
    void setBoundingBox(InsertBoundingBox & bb, double maxrad);

};

#endif // REGION_NEIGHBOR_LIST_H
