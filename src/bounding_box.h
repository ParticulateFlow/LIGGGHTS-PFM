/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2015-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_BOUNDING_BOX
#define LMP_BOUNDING_BOX

#include "mpi_liggghts.h"

namespace LAMMPS_NS
{

class BoundingBox
{
  friend class FixNeighlistMesh;

  public:

    BoundingBox();
    explicit BoundingBox(double bounds[6]);
    BoundingBox(double xLo, double xHi, double yLo, double yHi, double zLo, double zHi);
    virtual ~BoundingBox();

    void reset();

    void extendToContain(double const *pt)
    {
        if(initGiven){
            if(pt[0] < xLo) xLo = pt[0];
            else if(pt[0] > xHi) xHi = pt[0];

            if(pt[1] < yLo) yLo = pt[1];
            else if(pt[1] > yHi) yHi = pt[1];

            if(pt[2] < zLo) zLo = pt[2];
            else if(pt[2] > zHi) zHi = pt[2];
        } else{
          xLo = pt[0]; xHi = pt[0];
          yLo = pt[1]; yHi = pt[1];
          zLo = pt[2]; zHi = pt[2];
          initGiven = true;
        }
    }

    void extendToContain(BoundingBox const &other)
    {
      if(initGiven){
        if(other.xLo < xLo) xLo = other.xLo;
        if(other.xHi > xHi) xHi = other.xHi;

        if(other.yLo < yLo) yLo = other.yLo;
        if(other.yHi > yHi) yHi = other.yHi;

        if(other.zLo < zLo) zLo = other.zLo;
        if(other.zHi > zHi) zHi = other.zHi;
      } else{
        xLo = other.xLo; xHi = other.xHi;
        yLo = other.xLo; yHi = other.xHi;
        zLo = other.xLo; zHi = other.xHi;
      }
    }

    void extendToParallel(MPI_Comm comm)
    {
      double limit[6];
      limit[0] = -xLo;
      limit[1] =  xHi;
      limit[2] = -yLo;
      limit[3] =  yHi;
      limit[4] = -zLo;
      limit[5] =  zHi;

      MPI_Max_Vector(limit,6,comm);

      xLo = -limit[0];
      xHi =  limit[1];
      yLo = -limit[2];
      yHi =  limit[3];
      zLo = -limit[4];
      zHi =  limit[5];
    }

    void extrude(double length, const double * vec);

    void getBoxBounds(double *lo,double *hi)
    {
      lo[0] = xLo;
      lo[1] = yLo;
      lo[2] = zLo;
      hi[0] = xHi;
      hi[1] = yHi;
      hi[2] = zHi;
    }

    void getExtent(double extent[3]) const {
      extent[0] = xHi - xLo;
      extent[1] = yHi - yLo;
      extent[2] = zHi - zLo;
    }

    void getCenter(double center[3]) const {
      center[0] = 0.5 * (xHi + xLo);
      center[1] = 0.5 * (yHi + yLo);
      center[2] = 0.5 * (zHi + zLo);
    }

    bool hasVolume() const {
      return (xHi > xLo && yHi > yLo && zHi > zLo);
    }

    void extendByDelta(double delta);

    void getBoxBoundsExtendedByDelta(double *lo,double *hi,double delta)
    {
      lo[0] = xLo-delta;
      lo[1] = yLo-delta;
      lo[2] = zLo-delta;
      hi[0] = xHi+delta;
      hi[1] = yHi+delta;
      hi[2] = zHi+delta;
    }

    void shrinkToSubbox(double *sublo,double *subhi)
    {
      if(xLo < sublo[0])
          xLo = sublo[0];
      if(xHi > subhi[0])
          xHi = subhi[0];

      if(yLo < sublo[1])
          yLo = sublo[1];
      if(yHi > subhi[1])
          yHi = subhi[1];

      if(zLo < sublo[2])
          zLo = sublo[2];
      if(zHi > subhi[2])
          zHi = subhi[2];
    }

    bool isInitialized()
    { return initGiven; }

    bool isInside(const double *p)
    {
       // check bbox
       // test for >= and < as in Domain class
       return (p[0] >= xLo && p[0] < xHi &&
           p[1] >= yLo && p[1] < yHi &&
           p[2] >= zLo && p[2] < zHi);
    }

    bool isDirty() const {
        return dirty;
    }

    void setDirty(bool value) {
        dirty = value;
    }

  private:
    double xLo, xHi, yLo, yHi, zLo, zHi;
    bool initGiven;
    bool dirty;
};

} /* LAMMPS_NS */
#endif /* BOUNDINGBOX_H_ */
