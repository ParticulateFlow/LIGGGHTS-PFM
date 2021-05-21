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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_TRACKING_MESH_H
#define LMP_TRACKING_MESH_H

#include "multi_node_mesh_parallel.h"
#include "custom_value_tracker.h"
#include "container.h"
#include "memory.h"

namespace LAMMPS_NS{

  template<int NUM_NODES>
  class TrackingMesh : public MultiNodeMeshParallel<NUM_NODES>
  {
      public:

        void clearReverse();

        void setVerbose();

        virtual void buildNeighbours() = 0;

        virtual void move(const double *vecTotal, const double *vecIncremental);
        virtual void move(const double *vecIncremental);

        virtual void scale(double factor);

        virtual int generateRandomOwnedGhost(double *pos) = 0;
        virtual int generateRandomOwnedGhostWithin(double *pos,double delta) = 0;
        virtual int generateRandomSubbox(double *pos) = 0;

        // inline access

        inline CustomValueTracker& prop()
        {return customValues_;}

        // global to local ID mapping
        inline int map(int global)
        { return (mapArray_ ? mapArray_[global] : -1);}

        inline int tag_max()
        { return mapTagMax_; }

        int id_slow(int i)
        { return id_(i); }

        inline int id(int i)
        { return id_(i); }

        inline int lineNo(int i)
        { return (lineNo_?(*lineNo_)(i):-1); }

        inline bool verbose()
        { return verbose_; }

        inline void check_element_property_consistency()
        { customValues_.check_element_property_consistency(this->sizeLocal()+this->sizeGhost()); }

      protected:

        TrackingMesh(LAMMPS *lmp);
        virtual ~TrackingMesh();

        virtual bool addElement(double **nodeToAdd,int lineNumb);
        virtual void deleteElement(int n);

        void clearGhostForward(bool scale,bool translate,bool rotate);

        // called via mesh move functions
        bool resetToOrig();

        void postInitialSetup();

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        void clearMap();
        void generateMap();

        virtual void moveElement(int i, const double *vecIncremental);

        virtual void rotate(const double *totalQ, const double *dQ, const double *origin);
        virtual void rotate(const double *dQ, const double *origin);

        // buffer operations

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate,bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate,bool rotate);

        int meshPropsBufSize(int operation,bool scale,bool translate,bool rotate);
        int pushMeshPropsToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        int popMeshPropsFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

      private:

        // class holding fields
        CustomValueTracker &customValues_;

        // ID of element
        ScalarContainer<int> &id_;

        // line where element was read
        //NP ATTENTON is free'd after initialSetup
        ScalarContainer<int> *lineNo_;

        // global-local lookup
        int mapTagMax_;
        int *mapArray_;

        bool verbose_;

  };

  // *************************************
  #include "tracking_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* TRACKINGMESH_H_ */
