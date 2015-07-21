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

#ifndef LMP_SURFACE_MESH_I_H
#define LMP_SURFACE_MESH_I_H

#define NTRY_MC_SURFACE_MESH_I_H 30000
#define NITER_MC_SURFACE_MESH_I_H 5
#define TOLERANCE_MC_SURFACE_MESH_I_H 0.05

/*NL*/ #define DEBUGMODE_SURFACE_MESH false //true

/* ----------------------------------------------------------------------
   constructors, destructor
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::SurfaceMesh(LAMMPS *lmp)
:   TrackingMesh<NUM_NODES>(lmp),
    curvature_(1.-EPSILON_CURVATURE),
    minAngle_(cos(MIN_ANGLE_MESH*M_PI/180.)),

    // TODO should keep areaMeshSubdomain up-to-date more often for insertion faces
    areaMesh_     (*this->prop().template addGlobalProperty   < ScalarContainer<double> >                 ("areaMesh",     "comm_none","frame_trans_rot_invariant","restart_no",2)),

    nBelowAngle_(0),
    nTooManyNeighs_(0),
    nOverlapping_(0),

    //NP neigh topology is communicated at exchange and borders
    //NP  (neigh topology is created once and never changed)

    //NP no forward communication at all
    //NP if scale,move,translate: properties are manipulated/updated via CustomValueTracker

    area_         (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("area",         "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    areaAcc_      (*this->prop().template addElementProperty< ScalarContainer<double> >                   ("areaAcc",      "comm_none","frame_trans_rot_invariant", "restart_no",2)),
    edgeLen_      (*this->prop().template addElementProperty< VectorContainer<double,NUM_NODES> >         ("edgeLen",      "comm_none","frame_trans_rot_invariant", "restart_no")),
    edgeVec_      (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeVec",      "comm_none","frame_scale_trans_invariant","restart_no")),
    edgeNorm_     (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_NODES,3> >  ("edgeNorm",     "comm_none","frame_scale_trans_invariant","restart_no")),
    surfaceNorm_  (*this->prop().template addElementProperty< VectorContainer<double,3> >                 ("surfaceNorm",  "comm_none","frame_scale_trans_invariant","restart_no")),
    obtuseAngleIndex_   (*this->prop().template addElementProperty< ScalarContainer<int> >                ("obtuseAngleIndex","comm_exchange_borders","frame_invariant","restart_no")),
    nNeighs_      (*this->prop().template addElementProperty< ScalarContainer<int> >                      ("nNeighs",      "comm_exchange_borders","frame_invariant","restart_no")),
    neighFaces_   (*this->prop().template addElementProperty< VectorContainer<int,NUM_NEIGH_MAX> >        ("neighFaces",   "comm_exchange_borders","frame_invariant","restart_no")),
    hasNonCoplanarSharedNode_(*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >("hasNonCoplanarSharedNode","comm_exchange_borders","frame_invariant", "restart_no")),
    edgeActive_   (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("edgeActive",   "comm_exchange_borders","frame_invariant","restart_no")),
    cornerActive_ (*this->prop().template addElementProperty< VectorContainer<bool,NUM_NODES> >           ("cornerActive", "comm_exchange_borders","frame_invariant","restart_no"))

{
    //NP allocate 3 scalar spaces
    //NP add directly instead of setGlobalProperty
    //NP so no resetToOrig functionality (but not needed here)
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    areaMesh_.add(0.);
    /*NL*///this->error->all(FLERR,"check: use ID instead of index for neigh list, areCoplanar etc");
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::~SurfaceMesh()
{}

/* ----------------------------------------------------------------------
   set mesh curvature, used for mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::setCurvature(double _curvature)
{
    curvature_ = _curvature;
}

/* ----------------------------------------------------------------------
   add and delete an element
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::addElement(double **nodeToAdd,int lineNumb)
{
    if(TrackingMesh<NUM_NODES>::addElement(nodeToAdd,lineNumb))
    {

        //NP need to do this because some classes may access data before
        //NP setup() is called
        calcSurfPropertiesOfNewElement();
        return true;
    }
    return false;
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::refreshOwned(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshOwned(setupFlag);
    // (re)calculate all properties for owned elements
    //NP calculates properties for newly arrived elements
    //NP also removes round-off isues for moving mesh
    recalcLocalSurfProperties();
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::refreshGhosts(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshGhosts(setupFlag);

    recalcGhostSurfProperties();
}

/* ----------------------------------------------------------------------
   recalculate properties of local elements
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::recalcLocalSurfProperties()
{
    //NP could use this function instead of rotating
    //NP all properties
    //NP execute this after re-neighboring via
    //NP refresh() so no round-off issues occur

    // areaMeshGlobal [areaMesh_(0)] and areaMeshOwned [areaMesh_(1)]
    // calculated here

    areaMesh_(0) = 0.;
    areaMesh_(1) = 0.;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      for(int j=0;j<NUM_NODES;j++)
      {
          double dot;
          calcObtuseAngleIndex(i,j,dot);
      }

      area(i) = calcArea(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);

      // add to local area
      areaMesh_(1) += area(i);
      /*NL*///fprintf(this->screen,"triangle %d: area %f, areaacc %f, mesharea %f\n",i,area_(i),areaAcc_(i),areaMesh_);
    }

    // mesh area must be summed up
    MPI_Sum_Scalar(areaMesh_(1),areaMesh_(0),this->world);

    /*NL*/// fprintf(this->screen,"proc %d, areaMeshGlobal() %f,areaMeshOwned() %f,areaMeshGhost() %f\n",
    /*NL*///         this->comm->me,areaMeshGlobal(),areaMeshOwned(),areaMeshGhost());
    /*NL*/// this->error->all(FLERR,"check this");
}

/* ----------------------------------------------------------------------
   recalculate properties of ghost elements
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::recalcGhostSurfProperties()
{
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();

    // areaMeshGhost [areaMesh_(2)] calculated here

    // accumulated area includes owned and ghosts
    areaMesh_(2) = 0.;
    for(int i = nlocal; i < nall; i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));

      for(int j=0;j<NUM_NODES;j++)
      {
          double dot;
          calcObtuseAngleIndex(i,j,dot);
      }

      area(i) = calcArea(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);

      // add to ghost area
      areaMesh_(2) += area(i);
    }

    /*NL*/// fprintf(this->screen,"proc %d, areaMeshGlobal() %f,areaMeshOwned() %f,areaMeshGhost() %f\n",
    /*NL*///         this->comm->me,areaMeshGlobal(),areaMeshOwned(),areaMeshGhost());
    /*NL*/// this->error->all(FLERR,"check this");

    /*NL*///fprintf(this->screen,"proc %d: areaMeshOwned+Ghost %f areaAcc(lastGhost) %f SHOULD BE EQUAL\n",
    /*NL*///        this->comm->me,areaMeshOwned()+areaMeshGhost(),areaAcc(nall-1));
    /*NL*/// this->error->all(FLERR,"CHECK this");

    /*NL*///fprintf(this->screen,"proc %d: isInsertionMesh_ %s\n",this->comm->me,isInsertionMesh_?"true":"false");
    /*NL*/// this->error->all(FLERR,"CHECK this");

    /*NL*/ //if(this->map(21) >= 0) fprintf(this->screen,"proc %d has ID 21 and edgeActive(21)[1] is %s\n",this->comm->me,edgeActive(this->map(21))[1]?"y":"n");
}

/* ----------------------------------------------------------------------
   generate a random Element by areaAcc
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
inline int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::randomOwnedGhostElement()
{
    //NP disallow to use this unless this is an insertion mesh
    if(!this->isInsertionMesh())
        this->error->one(FLERR,"Illegal call for non-insertion mesh");

    //NP isInsertionMesh: parallel
    double area = areaMeshOwned()+areaMeshGhost();

    double r = this->random_->uniform() * area;
    /*NL*/ //fprintf(this->screen,"area %f\n",areaMeshGlobal());
    /*NL*/ //fprintf(this->screen,"areaMeshOwned()+areaMeshGhost() %f\n",areaMeshOwned()+areaMeshGhost());

    int first = 0;
    int last = this->sizeLocal()+this->sizeGhost()-1;

    return searchElementByAreaAcc(r,first,last);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
inline int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::searchElementByAreaAcc(double area,int lo, int hi)
{
    /*NL*/ //fprintf(this->screen,"areaAcc(lo) %f areaAcc(hi) %f\n",areaAcc(lo),areaAcc(hi));

    if( (lo < 1 || area > areaAcc(lo-1)) && (area <= areaAcc(lo)) )
        return lo;
    if( (hi < 1 || area > areaAcc(hi-1)) && (area <= areaAcc(hi)) )
        return hi;

    int mid = static_cast<int>((lo+hi)/2);
    if(area > areaAcc(mid))
        return searchElementByAreaAcc(area,mid,hi);
    else
        return searchElementByAreaAcc(area,lo,mid);
}

/* ----------------------------------------------------------------------
   calculate surface properties of new element
   only called once on import
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcSurfPropertiesOfNewElement()
{
    //NP IMPORTANT: do not use add() functions here
    //NP rather use set
    //NP this is b/c elements have been added already
    //NP in TrackingMesh::addElement()
    //NP via customValues_.grow()

    //NP always called in serial, so use sizeLocal()

    int n = this->sizeLocal()-1;

    double *vecTmp3,*vecTmpNumNodes,**nodeTmp;
    create<double>(vecTmp3,3);
    create<double>(vecTmpNumNodes,NUM_NODES);
    create<double>(nodeTmp,NUM_NODES,3);

    // calculate edge vectors and lengths
    calcEdgeVecLen(n,vecTmpNumNodes,nodeTmp);
    edgeLen_.set(n,vecTmpNumNodes);
    edgeVec_.set(n,nodeTmp);

    /*NP
    printVec3D(this->screen,"edge len", vecTmpNumNodes);
    printVec3D(this->screen,"edge vec0", nodeTmp[0]);
    printVec3D(this->screen,"edge vec1", nodeTmp[1]);
    printVec3D(this->screen,"edge vec2", nodeTmp[2]);

    printVec3D(this->screen,"edge len", edgeLen(n));
    printVec3D(this->screen,"edge vec0", edgeVec(n)[0]);
    printVec3D(this->screen,"edge vec1", edgeVec(n)[1]);
    printVec3D(this->screen,"edge vec2", edgeVec(n)[2]);

    printVec3D(this->screen,"edge len", edgeLen(n+1));
    printVec3D(this->screen,"edge vec0", edgeVec(n+1)[0]);
    printVec3D(this->screen,"edge vec1", edgeVec(n+1)[1]);
    printVec3D(this->screen,"edge vec2", edgeVec(n+1)[2]);*/

    // calc surface normal
    calcSurfaceNorm(n,vecTmp3);
    surfaceNorm_.set(n,vecTmp3);

    // calc edge normal in plane pointing outwards of area_
    // should be (edgeVec_ cross surfaceNormal)
    calcEdgeNormals(n,nodeTmp);
    edgeNorm_.set(n,nodeTmp);

    obtuseAngleIndex_.set(n,NO_OBTUSE_ANGLE);

    //NP check if triangle has an obtuse angle
    //NP and save its index
    //NP this is needed for contact detection
    //NP also check for highly skewed triangles
    bool hasSmallAngle = false;

    for(int i=0;i<NUM_NODES;i++){
      double dot;
      calcObtuseAngleIndex(n,i,dot);
      if(-dot > minAngle_)
        hasSmallAngle = true;
    }

    if(hasSmallAngle)
    {
        if(TrackingMesh<NUM_NODES>::verbose() && 0 == this->comm->me)
            fprintf(this->screen,"Mesh %s: elements %d (line %d) has high aspect ratio (angle < %f °) \n",
                    this->mesh_id_,n,TrackingMesh<NUM_NODES>::lineNo(n),this->angleLimit());
        nBelowAngle_++;
    }

    // calc area_ from previously obtained values and add to container
    // calcArea is pure virtual and implemented in derived class(es)
    //NP need not parallelize, every provess calculates areaMesh_ and areaAcc_
    double area_elem = calcArea(n);
    areaMesh_(0) += area_elem;
    area_(n) = area_elem;
    areaAcc_(n) = area_elem;
    if(n > 0) areaAcc_(n) += areaAcc_(n-1);

    // cannot calc areaMesh_(1), areaMesh_(2), areaMesh_(3) here since
    // not parallelized at this point
    //NP but is calculated anyway via refresh() and refreshGhosts()
    //NP in MultiNodeMeshParallel::setup() and initialSetup()

    /*NL*///fprintf(this->screen,"triangle %d: id %d,area %f, areaacc %f, mesharea %f\n",n,this->id_(n),area_(n),areaAcc_(n),areaMesh_(0));

    destroy<double>(nodeTmp);
    destroy<double>(vecTmpNumNodes);
    destroy<double>(vecTmp3);

    //NP neigh topology is initialized later - need not do this here
}

/* ----------------------------------------------------------------------
   sub-functions needed to calculate mesh properties
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcEdgeVecLen(int nElem, double *len, double **vec)
{
    for(int i=0;i<NUM_NODES;i++)
    {
      vectorSubtract3D(
        MultiNodeMesh<NUM_NODES>::node_(nElem)[(i+1)%NUM_NODES],
        MultiNodeMesh<NUM_NODES>::node_(nElem)[i],vec[i]);
      len[i] = vectorMag3D(vec[i]);
      vectorScalarDiv3D(vec[i],len[i]);
    }
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcSurfaceNorm(int nElem, double *surfNorm)
{
    vectorCross3D(edgeVec(nElem)[0],edgeVec(nElem)[1],surfNorm);
    vectorScalarDiv3D(surfNorm, vectorMag3D(surfNorm));
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcEdgeNormals(int nElem, double **edgeNorm)
{
    for(int i=0;i<NUM_NODES;i++){
      vectorCross3D(edgeVec(nElem)[i],surfaceNorm(nElem),edgeNorm[i]);
      vectorScalarDiv3D(edgeNorm[i],vectorMag3D(edgeNorm[i]));
    }
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::calcObtuseAngleIndex(int nElem, int iNode, double &dot)
{
    //NP this is not AB dot AC but AB dot CA --> swap sign
    dot = vectorDot3D(edgeVec_(nElem)[iNode],edgeVec_(nElem)[(iNode-1+NUM_NODES)%NUM_NODES]);
    /*NL*/ //fprintf(this->screen,"angle %f, minAngle_ %f\n",i,angle,minAngle_);
    if(dot > 0.)
    {
        /*NL*/ //if(n==12880) fprintf(this->screen,"TRI 12880 is obtuse at node %d, angle %f, minAngle_ %f\n",i,angle,minAngle_);
        obtuseAngleIndex_.set(nElem,iNode);
    }
    else
        obtuseAngleIndex_.set(nElem,NO_OBTUSE_ANGLE);
}

/* ----------------------------------------------------------------------
   build neighlist, generate mesh topology, check (in)active edges and nodes
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::buildNeighbours()
{
    /*NL*/ //fprintf(this->screen,"building neigh topology\n");

    // iterate over all surfaces, over ghosts as well
    //NP this is important for parallel correction!!
    int nall = this->sizeLocal()+this->sizeGhost();

    // inititalize neigh topology - reset to default, ~n
    bool t[NUM_NODES], f[NUM_NODES];
    int neighs[NUM_NEIGH_MAX];

    for(int i = 0; i < NUM_NODES; i++)
    {
        t[i] = true;
        f[i] = false;
    }
    for(int i = 0; i < NUM_NEIGH_MAX; i++)
        neighs[i] = -1;

    for(int i = 0; i < nall; i++)
    {
        nNeighs_.set(i,0);
        neighFaces_.set(i,neighs);
        edgeActive_.set(i,t);
        cornerActive_.set(i,t);
        hasNonCoplanarSharedNode_.set(i,f);
    }

    // build neigh topology and edge activity, ~n*n/2
    for(int i = 0; i < nall; i++)
    {
      for(int j = i+1; j < nall; j++)
      {
        //NP continue of do not share any node so can not share an edge
        int iEdge(0), jEdge(0);

        //NP assumption: 2 surface elements only share 1 edge at maximum
        //NP so for duplicate elements, only 1 edge is handled here!!
        if(shareEdge(i,j,iEdge,jEdge))
          handleSharedEdge(i,iEdge,j,jEdge, areCoplanar(TrackingMesh<NUM_NODES>::id(i),TrackingMesh<NUM_NODES>::id(j)));
      }
    }

    int *idListVisited = new int[nall];
    int *idListHasNode = new int[nall];
    double **edgeList,**edgeEndPoint;
    this->memory->create(edgeList,2*nall,3,"SurfaceMesh:edgeList");
    this->memory->create(edgeEndPoint,2*nall,3,"SurfaceMesh:edgeEndPoint");

    // recursively handle corner activity, ~n
    for(int i = 0; i < nall; i++)
    {
        for(int iNode = 0; iNode < NUM_NODES; iNode++)
            handleCorner(i,iNode,idListVisited,idListHasNode,edgeList,edgeEndPoint);
    }

    delete []idListVisited;
    delete []idListHasNode;
    this->memory->destroy(edgeList);
    this->memory->destroy(edgeEndPoint);
    /*NL*/ //fprintf(this->screen,"neigh end");

    // correct edge and corner activation/deactivation in parallel
    //NP this is to avoid false positives for cases where edges and corners are
    //NP not located at the proc where the element is owned
    parallelCorrection();
}

/* ----------------------------------------------------------------------
   quality check for surface mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::qualityCheck()
{
    // iterate over surfaces
    //NP do not need ghosts here
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();
    int me = this->comm->me;

    // check duplicate elements, n^2/2 operation
    //NP doing here makes it a local n^2/2 operation
    //NP checking local elements only is ok, since if they are
    //NP duplicate, they must be owned by same proc
    for(int i = 0; i < nlocal; i++)
    {
        for(int j = i+1; j < nall; j++)
        {
            if(this->nSharedNodes(i,j) == NUM_NODES)
            {
                fprintf(this->screen,"ERROR: Mesh %s: elements %d and %d (lines %d and %d) are duplicate\n",
                        this->mesh_id_,TrackingMesh<NUM_NODES>::id(i),TrackingMesh<NUM_NODES>::id(j),
                        TrackingMesh<NUM_NODES>::lineNo(i),TrackingMesh<NUM_NODES>::lineNo(j));
                if(!this->removeDuplicates())
                    this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. You can try re-running with 'heal auto_remove_duplicates'");
                else
                    this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. The mesh probably reached the precision you defined. "
                                           "You can try re-running with a lower value for 'precision'");
            }
        }
    }

    //NP check curvature for local elements
    for(int i = 0; i < nlocal; i++)
    {
      for(int iNode = 0; iNode < NUM_NODES; iNode++)
      {
        double dot;
        calcObtuseAngleIndex(i,iNode,dot);
        if(-dot > curvature_)
        {
            fprintf(this->screen,"ERROR: Mesh %s: The minumum angle of mesh element %d (line %d) is lower than the specified curvature. "
                                   "Increase mesh quality or decrease curvature (currently %f°)\n",
                                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(i),TrackingMesh<NUM_NODES>::lineNo(i),acos(curvature_)*180./M_PI);
            this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. You can try setting 'curvature' to 1e-5 or lower");
        }
      }
    }

    if(this->nBelowAngle() > 0 && 0 == me)
    {
        fprintf(this->screen,"Mesh %s: %d elements have high aspect ratio (angle < %f °)\n",
                this->mesh_id_,this->nBelowAngle(),this->angleLimit());
        this->error->warning(FLERR,"Fix mesh: Mesh contains highly skewed element, moving mesh (if used) will not parallelize well");
    }

    if(this->nTooManyNeighs() > 0 && 0 == me)
    {
        //NP warning message already printed to screen in handleEdge
        //NP no error there because could be caused by duplicate faces
        fprintf(this->screen,"Mesh %s: %d mesh elements have more than %d neighbors \n",
                this->mesh_id_,this->nTooManyNeighs(),NUM_NEIGH_MAX);
        this->error->one(FLERR,"Fix mesh: Bad mesh, cannot continue. Possibly corrupt elements with too many neighbors.\n"
                                "If you know what you're doing, you can try to change the definition of SurfaceMeshBase in tri_mesh.h and recompile");
    }

    if(nOverlapping() > 0)
    {
        fprintf(this->screen,"WARNING: Mesh %s: proc %d has %d element pairs that are coplanar, "
                "share an edge and overlap (but are not duplicate)\n",
                this->mesh_id_,me,nOverlapping());
    }
}

/* ----------------------------------------------------------------------
   correct edge and corner activation/deactivation in parallel
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::parallelCorrection()
{
    //NP global arrays for allreduce
    int iGlobal,iLocal;
    int mine = this->sizeLocal()+this->sizeGhost();
    int sizeGlob = this->sizeGlobal();
    int len = NUM_NODES*sizeGlob;

    /*NL*/ //fprintf(this->screen,"tag max %d, sizeGlob %d\n",this->tag_max(),sizeGlob);

    int *edgea = new int[len];
    int *cornera = new int[len];
    vectorInitializeN(edgea,len,2);
    vectorInitializeN(cornera,len,2);

    for(int i = 0; i < mine; i++)
    {
        iGlobal = TrackingMesh<NUM_NODES>::id(i);

        for(int j = 0; j < NUM_NODES; j++)
        {
            edgea[iGlobal*NUM_NODES+j] = edgeActive(i)[j]?1:0;
            cornera[iGlobal*NUM_NODES+j] = cornerActive(i)[j]?1:0;
        }
    }

    //NP allreduce to get a global minimum
    MPI_Min_Vector(edgea,len,this->world);
    MPI_Min_Vector(cornera,len,this->world);

    for(int i = 0; i < sizeGlob; i++)
    {
        iLocal = this->map(i);
        if(iLocal >= 0)
        {
            for(int j = 0; j < NUM_NODES; j++)
            {
                if(edgea[i*NUM_NODES+j] == 0)
                    edgeActive(iLocal)[j] = false;
                else if(edgea[i*NUM_NODES+j] == 1)
                    edgeActive(iLocal)[j] = true;
                else
                    this->error->one(FLERR,"Illegal situation in SurfaceMesh::parallelCorrection()");
                if(cornera[i*NUM_NODES+j] == 0)
                    cornerActive(iLocal)[j] = false;
                else if(cornera[i*NUM_NODES+j] == 1)
                    cornerActive(iLocal)[j] = true;
                else
                    this->error->one(FLERR,"Illegal situation in SurfaceMesh::parallelCorrection()");
            }
        }
    }

    delete []edgea;
    delete []cornera;
}

/* ----------------------------------------------------------------------
   functions to generate mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanar(int tag_a, int tag_b)
{
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanar()");

    // check if two faces are coplanar
    //NP used for building neigh topology, ONLY CALLED for neigh faces

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    /*NL*/ //if(127 == tag_a || 127 == tag_b){
    /*NL*/ // fprintf(this->screen,"a %d b %d  dot %f curvature_ %f\n",tag_a,tag_b, dot,curvature_);
    /*NL*/ // printVec3D(this->screen,"surfaceNorm(a)",surfaceNorm(a));
    /*NL*/ // printVec3D(this->screen,"surfaceNorm(b)",surfaceNorm(b));
    /*NL*/ // if(fabs(dot) > curvature_) fprintf(this->screen,"tag_a %d tag_b %d  are coplanar \n",tag_a,tag_b);
    /*NL*/ //}

    // need fabs in case surface normal is other direction
    if(fabs(dot) >= curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanarNeighs(int tag_a, int tag_b)
{
    bool areNeighs = false;
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanarNeighs()");

    // check if two faces are coplanar
    //NP used to transfer shear history btw planar faces

    // must be neighs, otherwise not considered coplanar
    for(int i = 0; i < nNeighs_(a); i++)
        if(neighFaces_(a)[i] == tag_b)
            areNeighs = true;

    if(!areNeighs) return false;

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    /*NL*/ //fprintf(this->screen,"a %d b %d  dot %f\n",a,b, dot);
    /*NL*/ //printVec3D(this->screen,"surfaceNorm(a)",surfaceNorm(a));
    /*NL*/ //printVec3D(this->screen,"surfaceNorm(b)",surfaceNorm(b));
    /*NL*/ //if(fabs(dot) > curvature_) fprintf(this->screen,"a %d b %d  are coplanar \n",a,b);
    /*NL*/ //else  fprintf(this->screen,"a %d b %d  are NOT coplanar \n",a,b);

    // need fabs in case surface normal is other direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::areCoplanarNodeNeighs(int tag_a, int tag_b)
{
    bool areNeighs = false;
    int a = this->map(tag_a);
    int b = this->map(tag_b);

    if(a < 0 || b < 0)
        this->error->one(FLERR,"Internal error: Illegal call to SurfaceMesh::areCoplanarNeighs()");

    // check if two faces are coplanar
    //NP used to transfer shear history btw planar faces

    // must be neighs, otherwise not considered coplanar
    for(int i = 0; i < nNeighs_(a); i++)
        if(neighFaces_(a)[i] == tag_b)
            areNeighs = true;

    //NP explicity calculate nSharedNodes here
    //NP could be pre-calculated, but is a rare case
    if(!areNeighs && MultiNodeMesh<NUM_NODES>::nSharedNodes(a,b) == 0) return false;

    double dot = vectorDot3D(surfaceNorm(a),surfaceNorm(b));
    /*NL*/ //fprintf(this->screen,"a %d b %d  dot %f\n",a,b, dot);
    /*NL*/// printVec3D(this->screen,"surfaceNorm(a)",surfaceNorm(a));
    /*NL*/// printVec3D(this->screen,"surfaceNorm(b)",surfaceNorm(b));
    /*NL*/ //if(fabs(dot) > curvature_) fprintf(this->screen,"a %d b %d  are coplanar \n",a,b);

    // need fabs in case surface normal is other direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::coplanarNeighsOverlap(int iSrf,int iEdge,int jSrf,int jEdge)
{
    //NP take node index iEdge of iSrf as reference point
    //NP one node of iSrf that is not in iEdge is node with index iEdge+2
    //NP one node of jSrf that is not in jEdge is node with index jEdge+2

    double vecI[3],vecJ[3], pRef[3], edgeN[3], dot1, dot2;

    vectorCopy3D(MultiNodeMesh<NUM_NODES>::node_(iSrf)[iEdge],pRef);
    vectorCopy3D(edgeNorm(iSrf)[iEdge],edgeN);

    vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node_(iSrf)[(iEdge+2)%NUM_NODES],pRef,vecI);
    vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node_(jSrf)[(jEdge+2)%NUM_NODES],pRef,vecJ);

    dot1 = vectorDot3D(vecI,edgeN);
    dot2 = vectorDot3D(vecJ,edgeN);

    /*NL*///fprintf(this->screen,"surfaces %d %d dot1*dot2 %f\n",iSrf,jSrf,dot1*dot2);

    if(dot1*dot2 > 0.)
    {
        if(TrackingMesh<NUM_NODES>::verbose())
        {
            //NP have info about lineNo only for local elements
            //NP in parallel, pair is handled by two procs so all info is written out
            int nlocal = this->sizeLocal();
            fprintf(this->screen,"WARNING: Mesh %s: elements %d and %d are coplanar, "
                    "share an edge and overlap (but are not duplicate)\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
            if(iSrf < nlocal)
                fprintf(this->screen,"INFO: Mesh %s: element %d corresponds to line # %d\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::lineNo(iSrf));
            if(jSrf < nlocal)
                fprintf(this->screen,"INFO: Mesh %s: element %d corresponds to line # %d\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(jSrf),TrackingMesh<NUM_NODES>::lineNo(jSrf));
        }

        nOverlapping_++;
        /*NL*///fprintf(this->screen,"i %d j %d share node: %s\n",iSrf,jSrf,MultiNodeMesh<NUM_NODES>::shareNode(iSrf,jSrf,a,b)?"yes":"no");
        //this->error->warning(FLERR,"Fix mesh: Check overlapping mesh elements");
        return true;
    }
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeVecsColinear(double *v,double *w)
{
    // need normalized vectors
    double dot = vectorDot3D(v,w);
    // need fabs in case vectors are in different direction
    if(fabs(dot) > curvature_) return true;
    else return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::growSurface(int iSrf, double by)
{
    double *tmp = new double[3];
    for(int i=0;i<NUM_NODES;i++)
    {
      vectorSubtract3D(MultiNodeMesh<NUM_NODES>::node(iSrf)[i],this->center_(iSrf),tmp);
      vectorScalarMult3D(tmp,by);
      vectorAdd3D(MultiNodeMesh<NUM_NODES>::node(iSrf)[i],
                      tmp,MultiNodeMesh<NUM_NODES>::node(iSrf)[i]);
    }
    delete[] tmp;
    return;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::shareEdge(int iSrf, int jSrf, int &iEdge, int &jEdge)
{
    int iNode1=0,jNode1=0,iNode2,jNode2;
    if(this->share2Nodes(iSrf,jSrf,iNode1,jNode1,iNode2,jNode2)){
      // following implementation of shareNode(), the only remaining option to
      // share an edge is that the next node of iSrf is equal to the next or previous
      // node if jSrf
      /*NL*/ //fprintf(this->screen,"surfaces %d %d do share nodes %d %d\n",iSrf,jSrf,i,j);

      //NP if nodes are 0 and 2, then edge index is 2
      if(2 == iNode1+iNode2)
        iEdge = 2;
      //NP otherwise them minimum node index is the edge index
      else
        iEdge = MathExtraLiggghts::min(iNode1,iNode2);

      //NP same for j
      if(2 == jNode1+jNode2)
        jEdge = 2;
      else
        jEdge = MathExtraLiggghts::min(jNode1,jNode2);

      return true;
    }
    /*NL*/ //fprintf(this->screen,"surfaces %d %d do not share edge\n",iSrf,jSrf);
    iEdge = -1; jEdge = -1;
    return false;
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::handleSharedEdge(int iSrf, int iEdge, int jSrf, int jEdge,
                                            bool coplanar, bool neighflag)
{
    //NP neighflag = true: build both neighbors and edge activity
    //NP neighflag = false: build edge activity only

    if(neighflag)
    {
        /*NL*/ //if(1971 == TrackingMesh<NUM_NODES>::id(iSrf))
        /*NL*/ //     fprintf(this->screen,"Mesh %s: element id %d (line %d) has element id %d as neigh\n",
        /*NL*/ //             this->mesh_id_,TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::lineNo(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
        /*NL*/ //if(1971 == TrackingMesh<NUM_NODES>::id(jSrf))
        /*NL*/ //     fprintf(this->screen,"Mesh %s: element id %d (line %d) has element id %d as neigh\n",
        /*NL*/ //             this->mesh_id_,TrackingMesh<NUM_NODES>::id(jSrf),TrackingMesh<NUM_NODES>::lineNo(jSrf),TrackingMesh<NUM_NODES>::id(iSrf));


        if(nNeighs_(iSrf) == NUM_NEIGH_MAX || nNeighs_(jSrf) == NUM_NEIGH_MAX)
        {
            int ii;
            if(nNeighs_(iSrf) == NUM_NEIGH_MAX)
                ii = iSrf;
            else
                ii = jSrf;

            nTooManyNeighs_++;
            fprintf(this->screen,"Mesh %s: element id %d (line %d) has %d neighs, but only %d expected\n",
                    this->mesh_id_,TrackingMesh<NUM_NODES>::id(ii),TrackingMesh<NUM_NODES>::lineNo(ii),nNeighs_(ii)+1,NUM_NEIGH_MAX);
            //NP no error here because could be caused by duplicate faces
            //NP so should throw error on duplicate faces first
        }

        // set neighbor topology
        if(nNeighs_(iSrf) < NUM_NEIGH_MAX)
            neighFaces_(iSrf)[nNeighs_(iSrf)] = TrackingMesh<NUM_NODES>::id(jSrf);
        if(nNeighs_(jSrf) < NUM_NEIGH_MAX)
            neighFaces_(jSrf)[nNeighs_(jSrf)] = TrackingMesh<NUM_NODES>::id(iSrf);
        nNeighs_(iSrf)++;
        nNeighs_(jSrf)++;
    }

    // deactivate one egde
    // other as well if coplanar
    //NP dont do this if mesh elements overlap
    //NP IMPORTANT have to use ID criterion in parallel because local i/j are different
    if(!coplanar || coplanarNeighsOverlap(iSrf,iEdge,jSrf,jEdge))
    {
        if(TrackingMesh<NUM_NODES>::id(iSrf) < TrackingMesh<NUM_NODES>::id(jSrf))
        {
            /*NL*/ //if(127 == TrackingMesh<NUM_NODES>::id(iSrf) || 127 == TrackingMesh<NUM_NODES>::id(jSrf))fprintf(this->screen,"non-coplanar %d and %d\n",TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
            edgeActive(iSrf)[iEdge] = false;
            edgeActive(jSrf)[jEdge] = true;
        }
        else
        {
            /*NL*/ //if(127 == TrackingMesh<NUM_NODES>::id(iSrf) || 127 == TrackingMesh<NUM_NODES>::id(jSrf))fprintf(this->screen,"non-coplanar %d and %d\n",TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
            edgeActive(iSrf)[iEdge] = true;
            edgeActive(jSrf)[jEdge] = false;
        }
    }
    else // coplanar
    {
        if(!coplanar) this->error->one(FLERR,"internal error");
        /*NL*/ //if(127 == TrackingMesh<NUM_NODES>::id(iSrf) || 127 == TrackingMesh<NUM_NODES>::id(jSrf))fprintf(this->screen,"Coplanar %d and %d\n",TrackingMesh<NUM_NODES>::id(iSrf),TrackingMesh<NUM_NODES>::id(jSrf));
        edgeActive(iSrf)[iEdge] = false;
        edgeActive(jSrf)[jEdge] = false;
    }
}

/* ---------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::handleCorner(int iSrf, int iNode,
        int *idListVisited,int *idListHasNode,double **edgeList,double **edgeEndPoint)
{
    double nodeToCheck[3];
    bool hasTwoColinearEdges, anyActiveEdge;
    int nIdListVisited = 0, nIdListHasNode = 0, maxId = -1, nEdgeList;

    this->node(iSrf,iNode,nodeToCheck);
    anyActiveEdge = false;
    checkNodeRecursive(iSrf,nodeToCheck,nIdListVisited,idListVisited,
        nIdListHasNode,idListHasNode,edgeList,edgeEndPoint,anyActiveEdge);

    // each element that shares the node contributes two edges
    nEdgeList = 2*nIdListHasNode;

    // get max ID
    for(int i = 0; i < nIdListHasNode; i++)
        maxId = MathExtraLiggghts::max(maxId,idListHasNode[i]);

    // check if any 2 edges coplanar
    //NP cross-check with edge end-point to avoid false positives
    hasTwoColinearEdges = false;
    for(int i = 0; i < nEdgeList; i++)
    {
        for(int j = i+1; j < nEdgeList; j++)
        {
            /*NL*///printVec3D(this->screen,"vec1",edgeList[i]);
            /*NL*///printVec3D(this->screen,"vec2",edgeList[j]);
            if(edgeVecsColinear(edgeList[i],edgeList[j]) && !this->nodesAreEqual(edgeEndPoint[i],edgeEndPoint[j]))
                hasTwoColinearEdges = true;
        }
    }

    /*NL*/ //fprintf(this->screen,"result mesh id %d, node %d, nNodeNeighs %d anyActiveEdge %s hasTwoColinearEdges %s is highest id sharing this node %s final %s \n",
    /*NL*/ //      TrackingMesh<NUM_NODES>::id(iSrf),iNode,nIdListHasNode,anyActiveEdge?"y":"n",hasTwoColinearEdges?"y":"n",
    /*NL*/ //       TrackingMesh<NUM_NODES>::id(iSrf) == maxId?"yes":"no",
    /*NL*/ //       ((hasTwoColinearEdges || !anyActiveEdge)||(TrackingMesh<NUM_NODES>::id(iSrf) != maxId))?"deact":"act");

    //NP deactiveate all "incarnations" of a corner if
    //NP (a) any 2 colinear edges OR
    //NP (b) all edges deactivated

    //NP deactivate all
    if(hasTwoColinearEdges || !anyActiveEdge)
        cornerActive(iSrf)[iNode] = false;
    //NP deactivate all but one, let the highest ID live
    else if(TrackingMesh<NUM_NODES>::id(iSrf) == maxId)
        cornerActive(iSrf)[iNode] = true;
    else
        cornerActive(iSrf)[iNode] = false;

    return nIdListHasNode;
}

/* ---------------------------------------------------------------------- */

//NP recursive algorithm returns anyActiveEdge
//NP anyActiveEdge = true if any of the edges ending at the node is active

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::checkNodeRecursive(int iSrf,double *nodeToCheck,
        int &nIdListVisited,int *idListVisited,int &nIdListHasNode,int *idListHasNode,
        double **edgeList,double **edgeEndPoint,bool &anyActiveEdge)
{
    int idNeigh, iNeigh, nEdgeList = 2*nIdListHasNode, nEdgeEndPoint = 2*nIdListHasNode;

    /*NL*/ //fprintf(this->screen,"called with id %d, looking for %f %f %f, nIdListVisited %d\n",TrackingMesh<NUM_NODES>::id(iSrf),nodeToCheck[0],nodeToCheck[1],nodeToCheck[2],nIdListVisited);
    /*NL*/ //for(int i = 0; i < nIdListVisited; i++) fprintf(this->screen,"was already at %d\n",idListVisited[i]);

    // check if I have been here already
    for(int i = 0; i < nIdListVisited; i++)
        if(idListVisited[i] == TrackingMesh<NUM_NODES>::id(iSrf)) return;

    // add to visited list
    idListVisited[nIdListVisited++] = TrackingMesh<NUM_NODES>::id(iSrf);

    // if contains node, add to list and call neighbors
    int iNode = this->containsNode(iSrf, nodeToCheck);
    if(iNode >= 0)
    {
        /*NL*/ //fprintf(this->screen," found node at id %d\n",TrackingMesh<NUM_NODES>::id(iSrf));

        idListHasNode[nIdListHasNode++] = TrackingMesh<NUM_NODES>::id(iSrf);
        // node iNode is associated with edge iNode and iNode-1
        vectorCopy3D(edgeVec(iSrf)[iNode],edgeList[nEdgeList++]);
        vectorCopy3D(edgeVec(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeList[nEdgeList++]);
        vectorCopy3D(this->node_(iSrf)[(iNode+1)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        vectorCopy3D(this->node_(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES],edgeEndPoint[nEdgeEndPoint++]);
        if(edgeActive(iSrf)[iNode]) anyActiveEdge = true;
        else if(edgeActive(iSrf)[(iNode-1+NUM_NODES)%NUM_NODES]) anyActiveEdge = true;

        // only call recursive if have neighbor and if I have data of neigh element (own or ghost)

        for(int iN = 0; iN < nNeighs_(iSrf); iN++)
        {
            idNeigh = neighFaces_(iSrf)[iN];
            if(idNeigh < 0) return;
            iNeigh = this->map(idNeigh);
            if(iNeigh >= 0)
                checkNodeRecursive(iNeigh,nodeToCheck,nIdListVisited,idListVisited,nIdListHasNode,
                                    idListHasNode,edgeList,edgeEndPoint,anyActiveEdge);
        }
    }
    /*NL*/ //else fprintf(this->screen," did not find node at id %d\n",TrackingMesh<NUM_NODES>::id(iSrf));
    /*NL*/ //fprintf(this->screen,"RET 1\n");
}


/* ----------------------------------------------------------------------
   move mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::move(double *vecTotal, double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecTotal,vecIncremental);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::move(double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecIncremental);
}

/* ----------------------------------------------------------------------
   scale mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::scale(double factor)
{
    TrackingMesh<NUM_NODES>::scale(factor);

    /*NP
    dont have to do this here
    areaMesh_(0) = 0.;
    for(int i=0;i<this->size();i++)
    {
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      area(i) = calcArea(i);
      areaMesh_(0) += area(i);
      areaAcc(i) = area(i);
      if(i > 0) areaAcc(i) += areaAcc(i-1);
    }
    */
}

/* ----------------------------------------------------------------------
   rotate mesh
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::rotate(double *totalQ, double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(totalQ,dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    /*NP
    for(int i=0;i<this->center_.size();i++)
    {
      printVec3D(this->screen,"edgeLen from autorotate",edgeLen(i));
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[0]);
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[1]);
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[2]);
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      printVec3D(this->screen,"edgeLen from recalc",edgeLen(i));
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[0]);
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[1]);
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[2]);
    }
    */
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
void SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::rotate(double *dQ,double *origin)
{
    TrackingMesh<NUM_NODES>::rotate(dQ,origin);

    // find out if rotating every property is cheaper than
    // re-calculating them from the new nodes
    /*NP
    dont have to do this here
    for(int i=0;i<this->center_.size();i++)
    {
      printVec3D(this->screen,"edgeLen from autorotate",edgeLen(i));
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[0]);
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[1]);
      printVec3D(this->screen,"edgeVec from autorotate",edgeVec(i)[2]);
      calcEdgeVecLen(i, edgeLen(i), edgeVec(i));
      calcSurfaceNorm(i, surfaceNorm(i));
      calcEdgeNormals(i, edgeNorm(i));
      printVec3D(this->screen,"edgeLen from recalc",edgeLen(i));
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[0]);
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[1]);
      printVec3D(this->screen,"edgeVec from recalc",edgeVec(i)[2]);
    }
    this->error->all(FLERR,"end");
    */
}

/* ----------------------------------------------------------------------
   check if faces is planar
   used to check if a face can be used for particle insertion
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::isPlanar()
{
    int id_j;
    int flag = 0;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
        if(flag) break;

        for(int ineigh = 0; ineigh < nNeighs_(i); ineigh++)
        {
            id_j = neighFaces_(i)[ineigh];
            if(!areCoplanarNeighs(TrackingMesh<NUM_NODES>::id(i),id_j))
                flag = 1;
        }
    }

    MPI_Max_Scalar(flag,this->world);

    if(flag) return false;
    return true;
}

/* ----------------------------------------------------------------------
   check if point on surface - only valid if pos is in my subbox
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
bool SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::isOnSurface(double *pos)
{
    bool on_surf = false;

    int nall = this->sizeLocal()+this->sizeGhost();

    // brute force
    // loop over ghosts as well as they might overlap my subbox
    for(int i = 0; i < nall; i++)
    {
        on_surf = on_surf || isInElement(pos,i);
        if(on_surf) break;
    }

    return on_surf;
}

/* ----------------------------------------------------------------------
   return number of active edges and corners for debugging
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::n_active_edges(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(edgeActive(i)[0]) n++;
    if(edgeActive(i)[1]) n++;
    if(edgeActive(i)[2]) n++;
    return n;
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
int SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::n_active_corners(int i)
{
    int n = 0;
    if(i > this->size()) return n;

    if(cornerActive(i)[0]) n++;
    if(cornerActive(i)[1]) n++;
    if(cornerActive(i)[2]) n++;
    return n;
}

/* ----------------------------------------------------------------------
   edge-edge, edge-node, edge-point distance
------------------------------------------------------------------------- */

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeEdgeDist(int iSrf, int iEdge, int jSrf, int jEdge)
{
    double d1,d2,d3,d4;
    d1 = edgeNodeDist(iSrf,iEdge,jSrf,jEdge);
    d2 = edgeNodeDist(iSrf,iEdge,jSrf,(jEdge+1)%NUM_NODES);
    d3 = edgeNodeDist(jSrf,jEdge,iSrf,(iEdge+1)%NUM_NODES);
    d4 = edgeNodeDist(jSrf,jEdge,iSrf,(iEdge+1)%NUM_NODES);
    return MathExtraLiggghts::min(d1,d2,d3,d4);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgeNodeDist(int iSrf, int iEdge, int jSrf, int jNode)
{
    return edgePointDist(iSrf, iEdge, MultiNodeMesh<NUM_NODES>::node_(jSrf)[jNode]);
}

template<int NUM_NODES, int NUM_NEIGH_MAX>
double SurfaceMesh<NUM_NODES,NUM_NEIGH_MAX>::edgePointDist(int iSrf, int iEdge, double *point)
{
        double nodeToP[3], dot;

        vectorSubtract3D(point,MultiNodeMesh<NUM_NODES>::node_(iSrf)[iEdge],nodeToP);
        dot = vectorDot3D(edgeVec(iSrf)[iEdge],nodeToP);

        //NP node iEdge is closest
        if(dot < 0)
            return vectorMag3D(nodeToP);
        //NP node iEdge+1 is closest
        else if(dot > edgeLen(iSrf)[iEdge])
        {
            vectorSubtract3D(point,MultiNodeMesh<NUM_NODES>::node_(iSrf)[(iEdge+1)%NUM_NODES],nodeToP);
            return vectorMag3D(nodeToP);
        }
        //NP between both nodes
        else
            return MathExtraLiggghts::abs(vectorDot3D(edgeNorm(iSrf)[iEdge],nodeToP));
}

#endif
