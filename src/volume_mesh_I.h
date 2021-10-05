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

#ifndef LMP_VOLUME_MESH_I_H
#define LMP_VOLUME_MESH_I_H

#define NTRY_MC_VOLUME_MESH_I_H 30000
#define NITER_MC_VOLUME_MESH_I_H 5
#define TOLERANCE_MC_VOLUME_MESH_I_H 0.05

/* ----------------------------------------------------------------------
   constructor(s), destructor
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::VolumeMesh(LAMMPS *lmp)
:   TrackingMesh<NUM_NODES>(lmp),

    // TODO should keep volMeshSubdomain up-to-date more often for insertion faces
    volMesh_(*this->prop().template addGlobalProperty < ScalarContainer<double> >  ("volMesh", "comm_none","frame_trans_rot_invariant","restart_no",3)),

    //NP neigh topology and volume is communicated at exchange and borders
    //NP  (neigh topology is created once and never changed)

    vol_            (*this->prop().template addElementProperty< ScalarContainer<double> >        ("vol",   "comm_exchange_borders","frame_trans_rot_invariant", "restart_no",3)),
    volAcc_         (*this->prop().template addElementProperty< ScalarContainer<double> >        ("volAcc","comm_exchange_borders","frame_trans_rot_invariant", "restart_no",3)),

    faceNodes_      (*this->prop().template addElementProperty< MultiVectorContainer<int,NUM_FACES,NUM_NODES_PER_FACE> > ("faceNodes ","comm_exchange_borders","frame_invariant", "restart_no")),
    faceNormals_    (*this->prop().template addElementProperty< MultiVectorContainer<double,NUM_FACES,3> >               ("faceNormals ","comm_none","frame_scale_trans_invariant", "restart_no")),
    isBoundaryFace_ (*this->prop().template addElementProperty< VectorContainer<bool,NUM_FACES> >                        ("isBoundaryFace","comm_exchange_borders","frame_invariant", "restart_no")),

    nNeighs_        (*this->prop().template addElementProperty< ScalarContainer<int> >           ("nNeighs",        "comm_exchange_borders","frame_invariant", "restart_no")),
    neighElems_     (*this->prop().template addElementProperty< VectorContainer<int,NUM_FACES> > ("neighElems",     "comm_exchange_borders","frame_invariant", "restart_no"))
{
    //NP allocate 4 scalar spaces
    volMesh_.add(0.);
    volMesh_.add(0.);
    volMesh_.add(0.);
    volMesh_.add(0.);

    this->error->all(FLERR,"have to add neigh topology with comm_exchange_borders");
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::~VolumeMesh()
{}

/* ----------------------------------------------------------------------
   add / delete element
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
bool VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::addElement(double **nodeToAdd)
{
    TrackingMesh<NUM_NODES>::addElement(nodeToAdd,-1);
    this->error->one(FLERR,"TODO line #");
    this->error->one(FLERR,"TODO auto remove dupl");

    //NP need to do this because some classes may access data before
    //NP setup() is called
    calcVolPropertiesOfNewElement();

    return true;
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::deleteElement(int n)
{
    TrackingMesh<NUM_NODES>::deleteElement(n);
}

/* ----------------------------------------------------------------------
   calculate properties when adding new element
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::calcVolPropertiesOfNewElement()
{
    //NP IMPORTANT: do not use add() functions here
    //NP rather use set
    //NP this is b/c elements have been added already
    //NP in TrackingMesh::addElement()
    //NP via customValues_.grow()

    int f0[3],f1[3],f2[3],f3[3];

    int n = MultiNodeMesh<NUM_NODES>::node_.size()-1;

    // flip node 3 if on the wrong side
    checkOrientation(n);

    // assign nodes to faces
    vectorConstruct3D(f0,0,1,2);
    vectorConstruct3D(f1,0,3,1);
    vectorConstruct3D(f2,0,2,3);
    vectorConstruct3D(f3,1,3,2);
    faceNodes_.set(n,0,f0);
    faceNodes_.set(n,1,f1);
    faceNodes_.set(n,2,f2);
    faceNodes_.set(n,3,f3);

    // calc face normals
    calcFaceNormals(n);

    // calc volume
    double vol_elem = calcVol(n);
    volMesh_(0) += vol_elem;
    vol_(n) = vol_elem;
    volAcc_(n) = vol_elem;
    if(n > 0) volAcc_(n) += volAcc_(n-1);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
inline void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::checkOrientation(int n)
{
    double v01[3],v02[3],v03[3],tmp[3];
    double **node = this->node_.begin()[n];

    vectorSubtract3D(node[1],node[0],v01);
    vectorSubtract3D(node[2],node[0],v02);
    vectorSubtract3D(node[3],node[0],v03);

    vectorCross3D(v01,v02,tmp);

    // if wrong orientation, switch node 0 and 1
    if(vectorDot3D(tmp,v03) > 0.)
    {
        vectorCopy3D(node[0],tmp);
        vectorCopy3D(node[1],node[0]);
        vectorCopy3D(tmp,node[1]);
    }
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::calcFaceNormals(int n)
{
    double v01[3],v02[3],fnormal[3];
    double **node = this->node_.begin()[n];
    int **facenodes = this->faceNodes_(n);

    for(int iFace = 0; iFace < NUM_FACES; iFace++)
    {
        vectorSubtract3D(node[facenodes[iFace][1]],node[facenodes[iFace][0]],v01);
        vectorSubtract3D(node[facenodes[iFace][2]],node[facenodes[iFace][0]],v02);
        vectorCross3D(v01,v02,fnormal);
        vectorNormalize3D(fnormal);

        faceNormals_.set(n,iFace,fnormal);
    }
}

/* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::refreshOwned(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshOwned(setupFlag);
    // (re)calculate all properties for owned elements
    //NP calculates properties for newly arrived elements
    //NP also removes round-off isues for moving mesh

    recalcLocalVolProperties();
}
template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::refreshGhosts(int setupFlag)
{
    TrackingMesh<NUM_NODES>::refreshGhosts(setupFlag);

    recalcGhostVolProperties();
}

/* ----------------------------------------------------------------------
   recalculate properties of local elements
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::recalcLocalVolProperties()
{
    //NP could use this function instead of rotating
    //NP all properties
    //NP execute this after re-neighboring via
    //NP refresh() so no round-off issues occur


    // volMeshGlobal [volMesh_(0)] and volMeshOwned [volMesh_(1)]
    // calculated here

    volMesh_(0) = 0.;
    volMesh_(1) = 0.;

    int nlocal = this->sizeLocal();

    for(int i = 0; i < nlocal; i++)
    {
      calcFaceNormals(i);

      vol(i) = calcVol(i);
      volAcc(i) = vol(i);
      if(i > 0) volAcc(i) += volAcc(i-1);

      // add to local volume
      volMesh_(1) += vol(i);
      /*NL*///if (this->screen) fprintf(this->screen,"triangle %d: area %f, areaacc %f, mesharea %f\n",i,area_(i),areaAcc_(i),areaMesh_);
    }

    // mesh vol must be summed up
    MPI_Sum_Scalar(volMesh_(1),volMesh_(0),this->world);


    /*NL*/// if (this->screen) fprintf(this->screen,"proc %d, areaMeshGlobal() %f,areaMeshOwned() %f,areaMeshGhost() %f\n",
    /*NL*///         this->comm->me,areaMeshGlobal(),areaMeshOwned(),areaMeshGhost());
    /*NL*/// this->error->all(FLERR,"check this");
}

/* ----------------------------------------------------------------------
   recalculate properties of ghost elements
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::recalcGhostVolProperties()
{
    double pos[3];
    int n_succ, n_iter;
    int nlocal = this->sizeLocal();
    int nall = this->sizeLocal()+this->sizeGhost();

    // volMeshGhost [volMesh_(2)] and volMeshSubdomain [volMesh_(3)]
    // calculated here

    // accumulated vol includes owned and ghosts
    volMesh_(2) = 0.;
    for(int i = nlocal; i < nall; i++)
    {
      calcFaceNormals(i);

      vol(i) = calcVol(i);
      volAcc(i) = vol(i);
      if(i > 0) volAcc(i) += volAcc(i-1);

      // add to ghost area
      volMesh_(2) += vol(i);
    }

    /*NL*/// if (this->screen) fprintf(this->screen,"proc %d, areaMeshGlobal() %f,areaMeshOwned() %f,areaMeshGhost() %f\n",
    /*NL*///         this->comm->me,areaMeshGlobal(),areaMeshOwned(),areaMeshGhost());
    /*NL*/// this->error->all(FLERR,"check this");

    /*NL*///if (this->screen) fprintf(this->screen,"proc %d: areaMeshOwned+Ghost %f areaAcc(lastGhost) %f SHOULD BE EQUAL\n",
    /*NL*///        this->comm->me,areaMeshOwned()+areaMeshGhost(),areaAcc(nall-1));
    /*NL*///this->error->all(FLERR,"CHECK this");

    // calc area of owned and ghost elements in my subdomain
    //NP use monte carlo
    volMesh_(3) = 0.;
    double volCheck = 0.;

    if(this->isInsertionMesh())
    {
        n_succ = 0;
        n_iter = 0;

        // iterate long enough so MC has the desired tolerance
        while( (n_iter < NITER_MC_VOLUME_MESH_I_H) &&
               (fabs((volCheck-volMeshGlobal()))/volMeshGlobal() > TOLERANCE_MC_VOLUME_MESH_I_H) )
        {
            // only generate random positions if I have any mesh elements
            if(nall)
            {
                for(int i = 0; i < NTRY_MC_VOLUME_MESH_I_H; i++)
                {
                    // pick a random position on owned or ghost element
                    if((generateRandomOwnedGhost(pos) >= 0) && (this->domain->is_in_subdomain(pos)))
                        n_succ++;
                }
            }
            n_iter++;
            volMesh_(3) = static_cast<double>(n_succ)/static_cast<double>(NTRY_MC_VOLUME_MESH_I_H*n_iter) * (volMeshOwned()+volMeshGhost());

            MPI_Sum_Scalar(volMesh_(3),volCheck,this->world);
        }

        if(fabs((volCheck-volMeshGlobal()))/volMeshGlobal() > TOLERANCE_MC_VOLUME_MESH_I_H)
            this->error->all(FLERR,"Local mesh volume calculation failed, try boosting NITER_MC_VOLUME_MESH_I_H");

        // correct so sum of all owned vols is equal to global area
        volMesh_(3) *= volMeshGlobal()/volCheck;
    }

    /*NL*///if (this->screen) fprintf(this->screen,"proc %d: sizeGlobal() %d, sizeLocal() %d, sizeGhost() %d,  areaMeshGlobal %f areaMeshOwned %f  areaMeshGhost %f areaMeshLocal %f n_iter %d tolerance %f\n",
    /*NL*///        this->comm->me,this->sizeGlobal(),this->sizeLocal(),this->sizeGhost(),areaMeshGlobal(),areaMeshOwned(),areaMeshGhost(),areaMeshLocal(),
    /*NL*///        n_iter,((areaCheck-areaMeshGlobal()))/areaMeshGlobal());
    /*NL*///this->error->all(FLERR,"CHECK this");
}

/* ----------------------------------------------------------------------
   generate a random Element by volAcc
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
int VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::randomOwnedGhostElement()
{
    //NP disallow to use this unless this is an insertion mesh
    if(!this->isInsertionMesh()) this->error->one(FLERR,"Illegal call for non-insertion mesh");
    double r = this->random_->uniform() * (volMeshOwned()+volMeshGhost());
    int nall = this->sizeLocal()+this->sizeGhost()-1;
    return searchElementByVolAcc(r,0,nall);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
int VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::searchElementByVolAcc(double vol,int lo, int hi)
{
    if( (lo < 1 || vol > volAcc(lo-1)) && (vol <= volAcc(lo)) )
        return lo;
    if( (hi < 1 || vol > volAcc(hi-1)) && (vol <= volAcc(hi)) )
        return hi;

    int mid = static_cast<int>((lo+hi)/2);
    if(vol > volAcc(mid))
        return searchElementByVolAcc(vol,mid,hi);
    else
        return searchElementByVolAcc(vol,lo,mid);
}

/* ----------------------------------------------------------------------
   build neighlist, generate mesh topology
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::buildNeighbours()
{
    int iFace,jFace;

    // iterate over all elems, over ghosts as well
    int nall = this->sizeLocal()+this->sizeGhost();

    // inititalize neigh topology - reset to default, ~n

    int neighs[NUM_FACES];
    bool isb[NUM_FACES];
    for(int i=0;i<NUM_FACES;i++)
    {
        neighs[i] = -1;
        isb[i] = true;
    }

    for(int i = 0; i < nall; i++)
    {
        nNeighs_.set(i,0);
        neighElems_.set(i,neighs);
        isBoundaryFace_.set(i,isb);
    }

    // build neigh topology, ~n*n/2
    for(int i = 0; i < nall; i++)
    {
      for(int j = i+1; j < nall; j++)
      {
        //NP continue of do not share any node
        if(0 == this->nSharedNodes(i,j)) continue;

        if(shareFace(i,j,iFace,jFace))
        {
            neighElems_(i)[nNeighs_(i)] = this->id(i);
            neighElems_(j)[nNeighs_(j)] = this->id(j);
            nNeighs_(i)++;
            nNeighs_(j)++;
            isBoundaryFace_(i)[iFace] = false;
            isBoundaryFace_(j)[jFace] = false;
        }
      }
    }
}

/* ----------------------------------------------------------------------
   isInside etc
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
bool VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::isInside(const double *p)
{
    // check subdomain
    if(!this->domain->is_in_subdomain(p)) return false;

    // check bbox
    if(!this->bbox_.isInside(p)) return false;

    int nall = this->size();

    // brute force
    for(int i = 0; i < nall; i++)
        if(isInside(i,p)) return true;

    return false;
}

/* ----------------------------------------------------------------------
   move, rotate, scale mesh
------------------------------------------------------------------------- */

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::move(const double *vecTotal, const double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecTotal,vecIncremental);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::move(const double *vecIncremental)
{
    TrackingMesh<NUM_NODES>::move(vecIncremental);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::scale(double factor)
{
    TrackingMesh<NUM_NODES>::scale(factor);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::rotate(const double *totalQ, const double *dQ, const double *totalDispl, const double *dDisp)
{
    TrackingMesh<NUM_NODES>::rotate(totalQ,dQ,totalDispl,dDisp);
}

template<int NUM_NODES,int NUM_FACES,int NUM_NODES_PER_FACE>
void VolumeMesh<NUM_NODES,NUM_FACES,NUM_NODES_PER_FACE>::rotate(const double *dQ, const double *dDispl)
{
    TrackingMesh<NUM_NODES>::rotate(dQ,dDispl);
}

#endif
