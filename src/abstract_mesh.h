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

#ifndef LMP_ABSTRACT_MESH_H
#define LMP_ABSTRACT_MESH_H

#include "pointers.h"

//NP this is a pure interface class
//NP reason is to be able to interface to derived classes
//NP without having to specifcy their concrete template
//NP args, so they can vary at runtime

namespace LAMMPS_NS
{
  class AbstractMesh : protected Pointers
  {
      friend class FixMoveMesh;
      friend class MeshMover;

      public:

        //NP ***************************************
        //NP interface to MultiNodeMesh
        //NP ***************************************

        virtual void setMeshID(const char *_mesh_id) = 0;

        virtual void setPrecision(double _precision) = 0;

        virtual void autoRemoveDuplicates() = 0;

        // scale mesh
        virtual void scale(double factor) = 0;

        // linear move w/ total and incremental displacement
        virtual void move(double *vecTotal, double *vecIncremental) = 0;

        // linear move w/ incremental displacement
        virtual void move(double *vecIncremental) = 0;

        // rotation w/ total and incremental displacement
        //   calls rotate(double *totalQuat,double *dQuat,double *displacement)
        virtual void rotate(double totalAngle, double dAngle, double *axis, double *p) = 0;

        // rotation w/ incremental displacement
        //   calls rotate(double *dQuat,double *displacement)
        virtual void rotate(double dAngle, double *axis, double *p) = 0;

        // rotation using quaternions
        virtual void rotate(double *totalQ, double *dQ,double *origin) = 0;

        // initialize movement
        virtual bool registerMove(bool _scale, bool _translate, bool _rotate) = 0;
        virtual void unregisterMove(bool _scale, bool _translate, bool _rotate) = 0;

        virtual bool isMoving() = 0;

        // get node j of element i
        //NP !!warning!!
        //NP this function is slow due to inevitable vtable lookup
        virtual void node_slow(int i,int j,double *node) = 0;

        // neigh list stuff for moving mesh
        virtual bool decideRebuild() = 0;

        //NP ***************************************
        //NP interface to MultiNodeMeshParallel
        //NP ***************************************

        virtual void initalSetup() = 0;
        virtual void pbcExchangeBorders(int setupFlag) = 0;
        virtual void clearReverse() = 0;
        virtual void forwardComm() = 0;
        virtual void reverseComm() = 0;

        virtual void writeRestart(FILE *fp) = 0;
        virtual void restart(double *list) = 0;

        virtual bool allNodesInsideSimulationBox() = 0;

        virtual int numNodes() = 0;

        virtual inline class CustomValueTracker& prop() = 0;

        //NP ***************************************
        //NP interface to TrackingMesh
        //NP ***************************************
        /*
        virtual ContainerBase* container(double type,int lenVec) = 0;
        virtual ContainerBase* container(int type,int lenVec) = 0;
        virtual ContainerBase* container(bool type,int lenVec) = 0;
        */
        virtual int id_slow(int i) = 0;

        virtual void setVerbose() = 0;

        //NP ***************************************
        //NP interface to SurfaceMesh
        //NP ***************************************

        virtual int nBelowAngle() = 0;
        virtual double angleLimit() = 0;
        virtual int nTooManyNeighs() = 0;

        //NP ***************************************
        //NP inline access functions for size
        //NP ***************************************

        // size includes owned and ghost elements
        inline int size()
        { return sizeLocal()+sizeGhost(); }

        // virtual functions for size
        // parallelism implemented in derived class
        virtual int sizeLocal() = 0;
        virtual int sizeGhost() = 0;
        virtual int sizeGlobal() = 0;


        //NP virtual destructor
        virtual ~AbstractMesh()
        {}

      protected:

          AbstractMesh(LAMMPS *lmp)
          : Pointers(lmp)
          {}

      private:

         virtual double*** nodePtr() = 0;
  };


} /* LAMMPS_NS */
#endif
