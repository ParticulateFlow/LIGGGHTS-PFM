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

#ifndef LMP_MULTI_NODE_MESH_H
#define LMP_MULTI_NODE_MESH_H

#include "abstract_mesh.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "math_extra_liggghts.h"
#include "update.h"
#include "container.h"
#include "bounding_box.h"
#include "random_park.h"

#define EPSILON_PRECISION 1e-8

namespace LAMMPS_NS
{
  template<int NUM_NODES>
  class MultiNodeMesh : public AbstractMesh
  {
      friend class FixMoveMesh;
      friend class MeshMover;

      public:

        void setMeshID(const char *_mesh_id);

        void setPrecision(double _precision);

        void setElementExclusionList(FILE *_file);

        void autoRemoveDuplicates();

        // scale mesh
        virtual void scale(double factor);

        // linear move w/ total and incremental displacement
        virtual void move(const double *vecTotal, const double *vecIncremental);

        // linear move w/ incremental displacement
        virtual void move(const double *vecIncremental);

        // rotation w/ total and incremental displacement
        //   calls rotate(double *totalQuat,double *dQuat,double *displacement)
        void rotate(double totalAngle, double dAngle, const double *axis, const double *p);

        // rotation w/ incremental displacement
        //   calls rotate(double *dQuat,double *displacement)
        void rotate(double dAngle, const double *axis, const double *p);

        //NP must be called after nodes are manipulated from outside
        //NP eg deforming mesh
        void updateCenterRbound(int ilo,int ihi);

        // initialize movement
        bool registerMove(bool _scale, bool _translate, bool _rotate);
        void unregisterMove(bool _scale, bool _translate, bool _rotate);

        // bbox stuff
        BoundingBox getGlobalBoundingBox();
        BoundingBox getElementBoundingBoxOnSubdomain(int const n);
        void updateGlobalBoundingBox();

        // neigh list stuff for moving mesh
        bool decideRebuild();
        void storeNodePosRebuild();

        // inline access

        inline bool isMoving()
        { return nMove_ > 0; }

        inline int nMove()
        { return nMove_; }

        inline bool isScaling()
        { return nScale_ > 0; }

        inline bool isTranslating()
        { return nTranslate_ > 0; }

        inline bool isRotating()
        { return nRotate_ > 0; }

        inline void node(int i,int j,double *node)
        { vectorCopy3D(node_(i)[j],node);}

        void node_slow(int i,int j,double *node)
        { vectorCopy3D(node_(i)[j],node);}

        inline void center(int i,double *center)
        { vectorCopy3D(center_(i),center);}

        inline int numNodes()
        { return NUM_NODES; }

        inline char* mesh_id()
        { return mesh_id_; }

        inline bool removeDuplicates()
        { return autoRemoveDuplicates_; }

        // virtual functions for size
        // parallelism implemented in derived class
        virtual int sizeLocal() = 0;
        virtual int sizeGhost() = 0;
        virtual int sizeGlobal() = 0;

        virtual bool isDeforming()
        { return false; }

      protected:
        MultiNodeMesh(LAMMPS *lmp);
        virtual ~MultiNodeMesh();


        virtual bool addElement(double **nodeToAdd);
        virtual void deleteElement(int n);

        virtual void refreshOwned(int setupFlag);
        virtual void refreshGhosts(int setupFlag);

        bool nodesAreEqual(int iSurf, int iNode, int jSurf, int jNode);
        bool nodesAreEqual(double *nodeToCheck1,double *nodeToCheck2);

        // returns true if surfaces share 2 nodes
        // called with local element indices
        // returns indices of nodes with iNode1 < iNode2
        bool share2Nodes(int iElem, int jElem, int &iNode1, int &jNode1, int &iNode2, int &jNode2);

        // returns number of shared nodes
        // called with local index
        int nSharedNodes(int iElem, int jElem);

        // returns node index if iElem contains nodeToCheck
        int containsNode(int iElem, double *nodeToCheck);

        void extendToElem(int const nElem) const;

        // linear move of single element w/ incremental displacement
        virtual void moveElement(int i, const double *vecIncremental);

        // rotation using quaternions
        //NP called by rotation functions above
        virtual void rotate(const double *totalQ, const double *dQ, const double *origin);
        virtual void rotate(const double *dQ, const double *origin);

        // mesh nodes
        MultiVectorContainer<double,NUM_NODES,3> node_;

        // original mesh node_ position, used for moving mesh
        MultiVectorContainer<double,NUM_NODES,3> *node_orig_;
        double** node_orig(int i) {return (*node_orig_)(i);}

        // node pos stored by external trigger at neigh build
        MultiVectorContainer<double,NUM_NODES,3> nodesLastRe_;

        // mesh center
        VectorContainer<double,3> center_; //NP used as reference for parallelization
        ScalarContainer<double> rBound_; //NP bounding radius

        // global bounding box for mesh across all processors
        BoundingBox bbox_;

        // random generator
        RanPark *random_;

        // mesh ID - same as fix mesh ID
        char *mesh_id_;

        //NP
        inline void reset_stepLastReset()
        { stepLastReset_ = -1; }

        // reset mesh nodes to original position
        // called via mesh move functions
        //NP resets node pos and global Props
        virtual bool resetToOrig();

        inline double precision()
        { return precision_; }

        inline FILE* elementExclusionList()
        { return element_exclusion_list_; }

      private:

        // mesh precision
        double precision_;

        FILE *element_exclusion_list_;

        // state if elements should be automatically removed if duplicate
        bool autoRemoveDuplicates_;

        // flags stating how many move operations are performed on the mesh
        int nMove_;
        int nScale_,nTranslate_,nRotate_;

        // store current node position for use by moving mesh
        void storeNodePosOrig(int ilo, int ihi);

        // step when nodes have been reset the last time
        // only relevant for moving mesh
        int stepLastReset_;

        // extends a given bbox to include element number nElem
        void extendToElem(BoundingBox &box, int const nElem);
        void extendToElem(int const nElem);

        inline double*** nodePtr()
        { return node_.begin(); }
  };

  // *************************************
  #include "multi_node_mesh_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* MULTINODEMESH_H_ */
