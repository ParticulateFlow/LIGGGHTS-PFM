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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_I_H

#define BIG_MNMP 1.0e20
#define BUFFACTOR_MNMP 1.5
#define BUFMIN_MNMP 2000
#define BUFEXTRA_MNMP 2000

/*NL*/ #define LMP_MULTI_NODE_MESH_PARALLEL_I_H_DEBUG false
/*NL*/ #define LMP_MULTI_NODE_MESH_PARALLEL_I_H_ID 88

  /* ----------------------------------------------------------------------
   consturctors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::MultiNodeMeshParallel(LAMMPS *lmp)
  : MultiNodeMesh<NUM_NODES>(lmp),
    nLocal_(0), nGhost_(0), nGlobal_(0), nGlobalOrig_(0),
    maxsend_(0), maxrecv_(0),
    buf_send_(0), buf_recv_(0),
    half_atom_cut_(0.),
    size_exchange_(0),
    size_forward_(0),
    size_border_(0),
    maxforward_(0),maxreverse_(0),
    nswap_(0),
    maxswap_(0),
    sendnum_(0),recvnum_(0),
    firstrecv_(0),
    sendproc_(0),recvproc_(0),
    size_forward_recv_(0),
    size_reverse_recv_(0),
    slablo_(0),slabhi_(0),
    sendlist_(0),
    maxsendlist_(0),
    pbc_flag_(0),
    pbc_(0)
  {
      // initialize comm buffers & exchange memory
      //NP as in Comm constructor

      maxsend_ = BUFMIN_MNMP;
      this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      maxrecv_ = BUFMIN_MNMP;
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");

      maxswap_ = 6;
      allocate_swap(maxswap_);

      sendlist_ = (int **) this->memory->smalloc(maxswap_*sizeof(int *),"MultiNodeMeshParallel:sendlist");
      this->memory->create(maxsendlist_,maxswap_,"MultiNodeMeshParallel:maxsendlist");
      for (int i = 0; i < maxswap_; i++) {
        maxsendlist_[i] = BUFMIN_MNMP;
        this->memory->create(sendlist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendlist[i]");
      }
  }

  /* ----------------------------------------------------------------------
   destructor
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  MultiNodeMeshParallel<NUM_NODES>::~MultiNodeMeshParallel()
  {
      free_swap();

      if (sendlist_)
        for (int i = 0; i < maxswap_; i++)
            this->memory->destroy(sendlist_[i]);

      this->memory->sfree(sendlist_);
      this->memory->destroy(maxsendlist_);

      this->memory->destroy(buf_send_);
      this->memory->destroy(buf_recv_);
  }

  /* ----------------------------------------------------------------------
   add and delete elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::addElement(double **nodeToAdd)
  {
    MultiNodeMesh<NUM_NODES>::addElement(nodeToAdd);
    nLocal_++;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteElement(int n)
  {
    //NP make sure the call is valid
    //NP may not delete local element while ghosts are present because
    //NP ghost would then be mixed into the list of local elements
    if(n < nLocal_ && nGhost_ != 0)
        this->error->one(FLERR,"Illegal call to MultiNodeMeshParallel<NUM_NODES>::deleteElement");

    MultiNodeMesh<NUM_NODES>::deleteElement(n);

    if(n >= nLocal_)
        nGhost_--;
    else
        nLocal_--;
  }

  /* ----------------------------------------------------------------------
   recalculate properties on setup (on start and during simulation)
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshOwned(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshOwned(setupFlag);
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::refreshGhosts(int setupFlag)
  {
    MultiNodeMesh<NUM_NODES>::refreshGhosts(setupFlag);
  }

  /* ----------------------------------------------------------------------
   completely clear ghosts - called in borders()
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhosts()
  {
      // delete ghost data from container classes

      /*NL*///MPI_Barrier(this->world);
      /*NL*/// int ndel = 0, ni = nGhost_;
      /*NL*/// fprintf(this->screen,"**mesh %s proc %d having %d ghosts\n",this->mesh_id_,this->comm->me,nGhost_);

      while(nGhost_ > 0)
      {
          /*NL*/// if(ndel > ni-10)fprintf(this->screen,"proc %d deleting ghost # %d of %d (element %d)\n",this->comm->me,ndel,ni,nLocal_);
          /*NL*/// ndel++;
          deleteElement(nLocal_);
      }
  }

  /* ----------------------------------------------------------------------
   clear ghost data that is communicated via forward comm - called in forw comm
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearGhostForward(bool scale,bool translate,bool rotate)
  {
      // delete ghost data from container classes
      // delete only data that is communicated afterwards

      for(int i = this->sizeLocal()+this->sizeGhost()-1; i >= this->sizeLocal(); i--)
      {
          // clear ghost data that belongs to this class
          // must match push/pop implementation for forward comm in this class
          if(translate || rotate || scale)
          {
            this->node_.del(i);
            this->center_.del(i);
          }
          if(scale)
            this->rBound_.del(i);
      }
  }

  /* ----------------------------------------------------------------------
   check if all elements are in domain
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  bool MultiNodeMeshParallel<NUM_NODES>::allNodesInsideSimulationBox()
  {
    int flag = 0;
    for(int i=0;i<sizeLocal();i++)
      for(int j=0;j<NUM_NODES;j++)
      {
        /*NL*/ //printVec3D(screen,"node ",node_(i)[j]);
        if(!this->domain->is_in_domain(this->node_(i)[j]))
        {
            flag = 1;
            break;
        }
      }

    MPI_Max_Scalar(flag,this->world);
    if(flag) return false;
    else return true;
  }

  /* ----------------------------------------------------------------------
   setup of communication
  ------------------------------------------------------------------------- */

   template<int NUM_NODES>
   void MultiNodeMeshParallel<NUM_NODES>::setup()
   {
       double sublo[3],subhi[3], cut, extent_acc;
       double rBound_max, cut_ghost;
       double **sublo_all, **subhi_all;

       int nprocs = this->comm->nprocs;
       int me = this->comm->me;
       int myloc[3], loc_dim, nextproc, need_this;

       // get required size of communication per element
       bool scale = this->isScaling();
       bool translate = this->isTranslating();
       bool rotate = this->isRotating();

       size_exchange_ = elemBufSize(OPERATION_COMM_EXCHANGE,scale,translate,rotate) + 1;
       size_border_ = elemBufSize(OPERATION_COMM_BORDERS,scale,translate,rotate);
       size_forward_ = elemBufSize(OPERATION_COMM_FORWARD,scale,translate,rotate);
       size_reverse_ = elemBufSize(OPERATION_COMM_REVERSE,scale,translate,rotate);

       // maxforward = # of datums in largest forward communication
       // maxreverse = # of datums in largest reverse communication

       //NP is directly equal since communication is per mesh
       maxforward_ = MathExtraLiggghts::max(size_exchange_,size_border_,size_forward_);
       maxreverse_ = size_reverse_;

       /*NP
       fprintf(this->screen,"size_forward_ %d\n",size_forward_);
       fprintf(this->screen,"size_border_ %d\n",size_border_);
       this->error->all(FLERR,"check if this is reasonable");
       */

       // copy comm and domain data
       vectorCopy3D(this->comm->myloc,myloc);
       vectorCopy3D(this->domain->sublo,sublo);
       vectorCopy3D(this->domain->subhi,subhi);

       this->memory->create(sublo_all,nprocs,3,"MultiNodeMeshParallel::setup() sublo_all");
       this->memory->create(subhi_all,nprocs,3,"MultiNodeMeshParallel::setup() subhi_all");

       // ghost elements are for computing interaction with owned particles
       // so need to aquire ghost elements that overlap my subbox extened by
       // half neigh cutoff
       //NP do not need to account for cutghostuser or extend_cut_ghost() by fixes
       //NP neighbor->cutneighmax is set in neighbor->init()
       half_atom_cut_ = this->neighbor->cutneighmax / 2.;

       //NP add half skin in case of moving mesh as mesh itself might move half skin
       //NP before rebuild
       if(this->isMoving())
         half_atom_cut_+= this->neighbor->skin / 2.;

       // calculate maximum bounding radius of elements across all procs
       rBound_max = 0.;
       for(int i = 0; i < sizeLocal(); i++)
           rBound_max = MathExtraLiggghts::max(this->rBound_(i),rBound_max);
      MPI_Max_Scalar(rBound_max,this->world);

       // mesh element ghost cutoff is element bounding radius plus half atom neigh cut
       cut_ghost = rBound_max + half_atom_cut_;

       // set up maxneed_, sendneed_
       // account for non-uniform boundaries due to load-balancing
       // so aquire sub-box bounds from all processors

       MPI_Allgather(sublo,3,MPI_DOUBLE,&(sublo_all[0][0]),3,MPI_DOUBLE,this->world);
       MPI_Allgather(subhi,3,MPI_DOUBLE,&(subhi_all[0][0]),3,MPI_DOUBLE,this->world);

       /*NL*///printVec3D(this->screen,"sublo",sublo);
       /*NL*///printVec3D(this->screen,"sublo_all",sublo_all[me]);
       /*NL*///printVec3D(this->screen,"subhi",subhi);
       /*NL*///printVec3D(this->screen,"subhi_all",subhi_all[me]);
       /*NL*///this->error->all(FLERR,"check if this is the same");


       // set up maxneed_ and sendneed_
       // assume element with max bound radius is in my subbox
       //NP do not know exacly if this is the case since exchange()
       //NP is called after setup()
       //NP TODO possible improvement: test which rbound my neighs are going
       //NP      to send to me
       //NP calculate how far away I could have to communicate
       //NP account for non-uniform subboxes

       for(int dim = 0; dim < 3; dim++)
       {
           bool is_x = dim == 0 ? true : false;
           bool is_y = dim == 1 ? true : false;
           bool is_z = dim == 2 ? true : false;

           // go each direction (N-S-E-W-UP-DN)
           //NP does similar thing than Comm::updown()

           maxneed_[dim] = 0;
           for(int way = -1; way <= 1; way += 2)
           {
               // start from location of myself
               // reset accumulated extent
               loc_dim = myloc[dim];
               extent_acc = 0.;
               need_this = 0;
               sendneed_[dim][way == -1 ? 0 : 1] = 0;

               while(extent_acc < cut_ghost)
               {
                   // increase or decrease location
                   loc_dim += way;

                   // break if at dead end and non-pbc
                   if( (loc_dim < 0 && !this->domain->periodicity[dim]) ||
                       (loc_dim > this->comm->procgrid[dim]-1 && !this->domain->periodicity[dim]) )
                           break;

                   // wrap around PBCs
                   if(loc_dim < 0 && this->domain->periodicity[dim])
                      loc_dim = this->comm->procgrid[dim]-1;

                   if(loc_dim > this->comm->procgrid[dim]-1)
                      loc_dim = 0;

                   // increase counters
                   need_this++;
                   sendneed_[dim][way == -1 ? 0 : 1]++;

                   // go to next proc in proc grid and add its extent
                   nextproc = this->comm->grid2proc[is_x ? loc_dim : myloc[0]]
                                                   [is_y ? loc_dim : myloc[1]]
                                                   [is_z ? loc_dim : myloc[2]];
                   extent_acc += subhi_all[nextproc][dim] - sublo_all[nextproc][dim];
               }

               maxneed_[dim] = MathExtraLiggghts::max(maxneed_[dim],need_this);
           }

           // limit maxneed for non-pbc
           //NP not sure if this is needed here, but it cannot harm...
           if(maxneed_[dim] > this->comm->procgrid[dim]-1 && !this->domain->periodicity[dim])
               maxneed_[dim] = this->comm->procgrid[dim]-1;
       }

       // maxneed_ summed accross all processors
      MPI_Max_Vector(maxneed_,3,this->world);

       destroy(sublo_all);
       destroy(subhi_all);

       /*NL*/// fprintf(this->screen,"proc %d: sendneed_[0]: %d %d  sendneed_[1]: %d %d  sendneed_[2]: %d %d   maxneed %d %d %d \n",
       /*NL*///            this->comm->me,sendneed_[0][0],sendneed_[0][1],sendneed_[1][0],sendneed_[1][1],sendneed_[2][0],sendneed_[2][1],
       /*NL*///            maxneed_[0],maxneed_[1],maxneed_[2]);
       /*NL*/// this->error->all(FLERR,"check if this is reasonable");

       /*NP
       fprintf(this->screen,"cut_ghost %f\n",cut_ghost);
       printVec3D(this->screen,"maxneed_",maxneed_);
       this->error->all(FLERR,"check if this is reasonable");
       */

       // allocate comm memory
       //NP 2* for going left and right

       nswap_ = 2 * (maxneed_[0]+maxneed_[1]+maxneed_[2]);
       if (nswap_ > maxswap_) grow_swap(nswap_);

       // setup parameters for each exchange:
       //   slablo_/slabhi_ = boundaries for slab of elements to send at each swap
       //   use -BIG/midpt/BIG to insure all elements included even if round-off occurs
       //   if round-off, atoms elements across PBC can be < or > than subbox boundary
       //   note that borders() only loops over subset of elements during each swap

       // treat all as PBC here, non-PBC is handled in borders() via r/s need[][]
       // pbc_flag_: 0 = nothing across a boundary, 1 = something across a boundary
       // pbc_ = -1/0/1 for PBC factor in each of 3/6 orthogonal/triclinic dirs
       // 1st part of if statement is sending to the west/south/down
       // 2nd part of if statement is sending to the east/north/up

       int dim,ineed;

       int iswap = 0;
       for (dim = 0; dim < 3; dim++)
       {
         for (ineed = 0; ineed < 2*maxneed_[dim]; ineed++)
         {
           pbc_flag_[iswap] = 0;
           vectorZeroizeN(pbc_[iswap],6);

           // send left, receive right
           if (ineed % 2 == 0)
           {
               sendproc_[iswap] = this->comm->procneigh[dim][0];
               recvproc_[iswap] = this->comm->procneigh[dim][1];

               if (ineed < 2) slablo_[iswap] = -BIG_MNMP;
               else slablo_[iswap] = 0.5 * (this->domain->sublo[dim] + this->domain->subhi[dim]);

               // use half cut here, since rBound is used (added) in checkBorderElement()
               slabhi_[iswap] = this->domain->sublo[dim] + half_atom_cut_;

               if (myloc[dim] == 0)
               {
                   pbc_flag_[iswap] = 1;
                   pbc_[iswap][dim] = 1;
               }
           }
           // send right, receive left
           else
           {
               sendproc_[iswap] = this->comm->procneigh[dim][1];
               recvproc_[iswap] = this->comm->procneigh[dim][0];

               // use half cut here, since rBound is used (added) in checkBorderElement()
               slablo_[iswap] = this->domain->subhi[dim] -  half_atom_cut_;
               if (ineed < 2) slabhi_[iswap] = BIG_MNMP;
               else slabhi_[iswap] = 0.5 * (this->domain->sublo[dim] + this->domain->subhi[dim]);

               if (myloc[dim] == this->comm->procgrid[dim]-1)
               {
                   pbc_flag_[iswap] = 1;
                   pbc_[iswap][dim] = -1;
               }
           }
           iswap++;
         }
       }
   }

/* ----------------------------------------------------------------------
   realloc the buffers needed for communication and swaps
------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_swap(int n)
  {
      free_swap();
      allocate_swap(n);

      sendlist_ = (int **)
        this->memory->srealloc(sendlist_,n*sizeof(int *),"MultiNodeMeshParallel:sendlist_");
      this->memory->grow(maxsendlist_,n,"MultiNodeMeshParallel:maxsendlist_");
      for (int i = maxswap_; i < n; i++)
      {
        maxsendlist_[i] = BUFMIN_MNMP;
        this->memory->create(sendlist_[i],BUFMIN_MNMP,"MultiNodeMeshParallel:sendlist_[i]");
      }
      maxswap_ = n;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::allocate_swap(int n)
  {
      this->memory->create(sendnum_,n,"MultiNodeMeshParallel:sendnum_");
      this->memory->create(recvnum_,n,"MultiNodeMeshParallel:recvnum_");
      this->memory->create(sendproc_,n,"MultiNodeMeshParallel:sendproc_");
      this->memory->create(recvproc_,n,"MultiNodeMeshParallel:recvproc_");
      this->memory->create(size_forward_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(size_reverse_recv_,n,"MultiNodeMeshParallel:size");
      this->memory->create(slablo_,n,"MultiNodeMeshParallel:slablo_");
      this->memory->create(slabhi_,n,"MultiNodeMeshParallel:slabhi_");
      this->memory->create(firstrecv_,n,"MultiNodeMeshParallel:firstrecv");
      this->memory->create(pbc_flag_,n,"MultiNodeMeshParallel:pbc_flag_");
      this->memory->create(pbc_,n,6,"MultiNodeMeshParallel:pbc_");
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::free_swap()
  {
      this->memory->destroy(sendnum_);
      this->memory->destroy(recvnum_);
      this->memory->destroy(sendproc_);
      this->memory->destroy(recvproc_);
      this->memory->destroy(size_forward_recv_);
      this->memory->destroy(size_reverse_recv_);
      this->memory->destroy(slablo_);
      this->memory->destroy(slabhi_);
      this->memory->destroy(firstrecv_);
      this->memory->destroy(pbc_flag_);
      this->memory->destroy(pbc_);
  }

  /* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_send(int n, int flag)
  {
      maxsend_ = static_cast<int> (BUFFACTOR_MNMP * n);
      /*NL*///fprintf(this->screen,"grow_send at proc %d with flag %d, maxsend_ %d\n",this->comm->me,flag,maxsend_);
      if (flag)
        this->memory->grow(buf_send_,(maxsend_+BUFEXTRA_MNMP),"MultiNodeMeshParallel:buf_send");
      else {
        this->memory->destroy(buf_send_);
        this->memory->create(buf_send_,maxsend_+BUFEXTRA_MNMP,"MultiNodeMeshParallel:buf_send");
      }
  }

  /* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_recv(int n)
  {
      maxrecv_ = static_cast<int> (BUFFACTOR_MNMP * n);
      this->memory->destroy(buf_recv_);
      this->memory->create(buf_recv_,maxrecv_,"MultiNodeMeshParallel:buf_recv");
  }

  /* ----------------------------------------------------------------------
   realloc the size of the iswap sendlist as needed with BUFFACTOR
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::grow_list(int iswap, int n)
  {
    maxsendlist_[iswap] = static_cast<int> (BUFFACTOR_MNMP * n)+1;
    this->memory->grow(sendlist_[iswap],maxsendlist_[iswap],"MultiNodeMeshParallel:sendlist[iswap]");
  }

  /* ----------------------------------------------------------------------
   parallelization -
   initially, all processes have read the whole data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::initalSetup()
  {
      nGlobalOrig_ = sizeLocal();

      // delete all elements that do not belong to this processor
      deleteUnowned();

      if(sizeGlobal() != sizeGlobalOrig())
      {
        /*NL*/ fprintf(this->screen,"orig %d now %d\n", sizeGlobalOrig(),sizeGlobal());
        this->error->all(FLERR,"Mesh elements have been lost");
      }

      // set-up mesh parallelism
      setup();

      // re-calculate properties for owned particles
      //NP i.e. re-calculate from nodes, removing round-off issues for moving mesh
      refreshOwned(1);

      // identify elements that are near borders
      // forward communicate them
      //NP forward-communicates all properties which have been calculated before
      borders();

      // re-calculate properties for ghost particles
      refreshGhosts(1);

      // build mesh topology and neigh list
      //NP operations performed in parallel at this point
      //NP already have IDs at this point which is important
      buildNeighbours();

      /*NL*/// fprintf(this->screen,"INITIALSETUP: proc %d, mesh %s - nLocal %d, nGhost %d\n",
      /*NL*///                      this->comm->me,this->mesh_id_,nLocal_,nGhost_);

      //NP TODO: calc properties of all new elements, including ghosts
      //NP not sure if have to do this since borders communicates everything
      /*NL*///if(this->comm->nprocs > 1) this->error->one(FLERR,"TODO");
  }

  /* ----------------------------------------------------------------------
   parallelization - aggregates pbc, exchange and borders
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbcExchangeBorders(int setupFlag)
  {
      // need not do this during simulation for non-moving mesh and non-changing simulation box
      //NP in this case borders() and refresh() is done once
      //NP via parallelize()

      //NP have to reset to nodeOrig twice for steps where consecutive runs overlap
      //NP this is for consecutive runs
      if(setupFlag) this->reset_stepLastReset();

      //NP check for setupFlag since some mesh properties might have changed
      //NP so need to refresh everything once
      if(!setupFlag && !this->isMoving() && !this->domain->box_change) return;

      // set-up mesh parallelism
      setup();

      // enforce pbc
      pbc();

      // communicate particles
      exchange();

      if(sizeGlobal() != sizeGlobalOrig())
      {
        /*NL*/ fprintf(this->screen,"orig %d now %d\n", sizeGlobalOrig(),sizeGlobal());
        this->error->all(FLERR,"Mesh elements have been lost");
      }

      // re-calculate properties for owned particles
      //NP i.e. re-calculate from nodes, removing round-off issues for moving mesh
      refreshOwned(setupFlag);

      // identify elements that are near borders
      // forward communicate them
      //NP forward-communicates all properties which have been calculated before
      borders();

      // re-calculate properties for ghosts
      refreshGhosts(setupFlag);

      /*NL*/// fprintf(this->screen,"PBCEXCHBRDRS: proc %d, mesh %s - nLocal %d, nGhost %d\n",
      /*NL*///                      this->comm->me,this->mesh_id_,nLocal_,nGhost_);

      //NP build mesh topology and neigh list
      //NP operations performed in parallel at this point
      //NP need not do this because neighs are based on IDs and are
      //NP forward communicated
      /*NL*///buildNeighbours();
  }

  /* ----------------------------------------------------------------------
   parallelization - clear data of reverse comm properties
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::clearReverse()
  {
      // nothing to do here
  }

  /* ----------------------------------------------------------------------
   delete all particles which are not owned on this proc
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::deleteUnowned()
  {
      /*NL*///fprintf(this->screen,"proc %d has %d owned elements (is equal to nglobal here)\n",this->comm->me,nLocal_);

      int i = 0;

      while(i < nLocal_)
      {
          if(!this->domain->is_in_subdomain(this->center_(i)))
              this->deleteElement(i);
          else i++;
      }

      // calculate nGlobal for the first time
     MPI_Sum_Scalar(nLocal_,nGlobal_,this->world);

      /*NL*/// fprintf(this->screen,"proc %d has %d owned elements, nglobal is %d\n",this->comm->me,nLocal_,nGlobal_);
      /*NL*/// this->error->all(FLERR,"check this\n");
  }

  /* ----------------------------------------------------------------------
   enforce periodic boundary conditions
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::pbc()
  {
      double centerNew[3], delta[3];

      for(int i = 0; i < this->sizeLocal(); i++)
      {
          vectorCopy3D(this->center_(i),centerNew);
          this->domain->remap(centerNew);
          vectorSubtract3D(centerNew,this->center_(i),delta);

          // move element i incremental
          if(vectorMag3DSquared(delta) > 1e-9)
            this->moveElement(i,delta);
      }
  }

  /* ----------------------------------------------------------------------
   exchange elements with nearby processors
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::exchange()
  {
      int nrecv, nsend = 0;
      int nrecv1,nrecv2;
      double *buf;
      MPI_Request request;
      MPI_Status status;
      MPI_Comm world = this->world;

      int nprocs = this->comm->nprocs;
      int *procgrid = this->comm->procgrid;
      int procneigh[3][2];

      // clear global->local map for owned and ghost atoms
      //NP b/c elements migrate to new procs in exchange() and
      //NP new ghosts are created in borders()
      //NP generateMap() is done at end of borders()
      clearMap();

      // clear old ghosts
      //NP will be re-created in borders()
      clearGhosts();

      // copy procneigh
      for (int i = 0; i < 3; i++)
        for( int j = 0; j < 2; j++)
            procneigh[i][j] = this->comm->procneigh[i][j];

      for (int dim = 0; dim < 3; dim++)
      {
          // push data to buffer
          //NP fill buffer with elements leaving my box, using < and >=
          //NP when elements is deleted, fill it in with last atom

          nsend = pushExchange(dim);

          /*NL*/ //fprintf(this->screen,"proc %d pushing %d doubles to buf, dim %d on step %d\n",this->comm->me,nsend,dim,this->update->ntimestep);

          // send/recv in both directions
          // if 1 proc in dimension, no send/recv, set recv buf to send buf
          // if 2 procs in dimension, single send/recv
          // if more than 2 procs in dimension, send/recv to both neighbors

          if (procgrid[dim] == 1)
          {
            nrecv = nsend;
            buf = buf_send_;
          }
          else
          {
            MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,&nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
            nrecv = nrecv1;

            /*NL*/ //fprintf(this->screen,"nrecv1 %d on proc %d, dim %d on step %d\n",nrecv1,this->comm->me,dim,this->update->ntimestep);

            if (this->comm->procgrid[dim] > 2)
            {
                MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,&nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
                nrecv += nrecv2;
            }

            if (nrecv > maxrecv_) grow_recv(nrecv);

            MPI_Irecv(buf_recv_,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,world,&request);
            MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
            MPI_Wait(&request,&status);

            if (procgrid[dim] > 2)
            {
                MPI_Irecv(&buf_recv_[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,world,&request);
                MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
                MPI_Wait(&request,&status);
            }

            buf = buf_recv_;
          }

          // check incoming elements to see if they are in my box
          // if so, add on this proc

          /*NL*/// fprintf(this->screen,"proc %d has %d owned elems, dim %d on step %d\n",this->comm->me,nLocal_,dim,this->update->ntimestep);
          /*NL*/// fprintf(this->screen,"received %d doubles on proc %d, dim %d on step %d\n",nrecv,this->comm->me,dim,this->update->ntimestep);
          popExchange(nrecv,dim, buf);
          /*NL*/// fprintf(this->screen,"proc %d has %d owned elems, dim %d on step %d\n",this->comm->me,nLocal_,dim,this->update->ntimestep);
      }

      // re-calculate nGlobal as some element might have been lost
     MPI_Sum_Scalar(nLocal_,nGlobal_,world);
  }

  /* ----------------------------------------------------------------------
   generate ghost elements, refresh global map
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::borders()
  {
      int iswap, twoneed, nfirst, nlast, n, nsend, nrecv, smax, rmax;
      bool sendflag, dummy = false;
      double *buf;
      double lo,hi;
      MPI_Request request;
      MPI_Status status;

      iswap = 0;
      smax = rmax = 0;

      for (int dim = 0; dim < 3; dim++)
      {
          nlast = 0;

          // need to go left and right in each dim
          twoneed = 2*maxneed_[dim];
          for (int ineed = 0; ineed < twoneed; ineed++)
          {
              lo = slablo_[iswap];
              hi = slabhi_[iswap];

              // find elements within slab boundaries lo/hi using <= and >=
              //NP   for first swaps in a dim, check owned and ghost
              //NP   for later swaps in a dim, only check newly arrived ghosts
              //NP store sent atom indices in list for use in future timesteps

              if (ineed % 2 == 0)
              {
                  nfirst = nlast;
                  nlast = sizeLocal() + sizeGhost();
              }

              nsend = 0;

              // sendflag = 0 if I do not send on this swap
              //NP deviating from Comm, it is only zero if would send
              //NP through periodic or non-periodic BC

              //NP TODO: differentiate between non-PBC (do not send)
              //NP       PBC (send). if send through PBC, add box offset
              //NP       to ghost position

              /*NP OLD - from Comm
              if (ineed/2 >= sendneed_[dim][ineed % 2])
                  sendflag = false;
              else sendflag = true;
              */

              sendflag = true;
              //NP comm left and on left end of proc grid
              if(ineed % 2 == 0 && this->comm->myloc[dim] == 0)
                sendflag = false;

              //NP comm right and on right end of proc grid
              if(ineed % 2 == 1 && this->comm->myloc[dim] == this->comm->procgrid[dim]-1)
                sendflag = false;

              /*NL*/// if(this->comm->me == 1 && iswap == 2)
              /*NL*///    fprintf(this->screen,"   mesh %s: proc %d swap # %d: sendflag %s\n",this->mesh_id_,this->comm->me,iswap,sendflag?"true":"false");

              // find send elements
              if(sendflag)
              {
                  /*NL*/// if(true/*this->comm->me == 1 && iswap == 2*/)
                  /*NL*///    fprintf(this->screen,"   mesh %s: proc %d checks elements from %d to %d for send on swap # %d\n",this->mesh_id_,this->comm->me,nfirst,nlast,iswap);

                  for (int i = nfirst; i < nlast; i++)
                  {
                      if( ((ineed % 2 == 0) && checkBorderElementLeft(i,dim,lo,hi))  ||
                          ((ineed % 2 != 0) && checkBorderElementRight(i,dim,lo,hi))  )
                      {
                          if (nsend >= maxsendlist_[iswap])
                              grow_list(iswap,nsend);
                          sendlist_[iswap][nsend++] = i;

                          /*NL*/ if(LMP_MULTI_NODE_MESH_PARALLEL_I_H_DEBUG && LMP_MULTI_NODE_MESH_PARALLEL_I_H_ID == id(i))
                          /*NL*/    fprintf(this->screen,"     mesh %s: proc %d prepares element id %d for send on swap # %d\n",this->mesh_id_,this->comm->me,id(i),iswap);
                      }
                  }
              }

              /*NL*/// fprintf(this->screen,"mesh %s: proc %d sends %d elements on swap # %d\n",this->mesh_id_,this->comm->me,nsend,iswap);

              // pack up list of border elements

              if(nsend*size_border_ > maxsend_)
                grow_send(nsend*size_border_,0);

              n = pushElemListToBuffer(nsend, sendlist_[iswap], buf_send_, OPERATION_COMM_BORDERS,dummy,dummy,dummy);

              // swap atoms with other proc
              // no MPI calls except SendRecv if nsend/nrecv = 0
              // put incoming ghosts at end of my atom arrays
              // if swapping with self, simply copy, no messages

              if (sendproc_[iswap] != this->comm->me)
              {
                  MPI_Sendrecv(&nsend,1,MPI_INT,sendproc_[iswap],0,&nrecv,1,MPI_INT,recvproc_[iswap],0,this->world,&status);
                  if (nrecv*size_border_ > maxrecv_)
                      grow_recv(nrecv*size_border_);
                  if (nrecv)
                      MPI_Irecv(buf_recv_,nrecv*size_border_,MPI_DOUBLE,recvproc_[iswap],0,this->world,&request);

                  if (n)
                      MPI_Send(buf_send_,n,MPI_DOUBLE,sendproc_[iswap],0,this->world);

                  if (nrecv)
                      MPI_Wait(&request,&status);

                  buf = buf_recv_;
              }
              else
              {
                  nrecv = nsend;
                  buf = buf_send_;
              }

              // unpack buffer

              /*NL*/// fprintf(this->screen,"mesh %s: proc %d received %d elements on swap # %d\n",this->mesh_id_,this->comm->me,nrecv,iswap);

              n = popElemListFromBuffer(nLocal_+nGhost_,nrecv,buf_recv_,OPERATION_COMM_BORDERS,dummy,dummy,dummy);

              /*NL*/// if(LMP_MULTI_NODE_MESH_PARALLEL_I_H_DEBUG && LMP_MULTI_NODE_MESH_PARALLEL_I_H_ID == id(nLocal_+nGhost_-1))
              /*NL*///    fprintf(this->screen,"     mesh %s: proc %d recvs element id %d on swap # %d\n",this->mesh_id_,this->comm->me,id(nLocal_+nGhost_-1),iswap);

              // set pointers & counters

              smax = MAX(smax,nsend);
              rmax = MAX(rmax,nrecv);
              sendnum_[iswap] = nsend;
              recvnum_[iswap] = nrecv;
              size_forward_recv_[iswap] = nrecv*size_forward_;
              size_reverse_recv_[iswap] = nsend*size_reverse_;
              firstrecv_[iswap] = nLocal_+nGhost_;
              nGhost_ += nrecv;
              iswap++;
          }
      }

      // insure send/recv buffers are long enough for all forward & reverse comm
      int max = MAX(maxforward_*smax,maxreverse_*rmax);
      if (max > maxsend_) grow_send(max,0);
      max = MAX(maxforward_*rmax,maxreverse_*smax);
      if (max > maxrecv_) grow_recv(max);

      // build global-local map
      this->generateMap();
  }

  /* ----------------------------------------------------------------------
   check if element qualifies as ghost
  ------------------------------------------------------------------------- */

  //NP IMPORTANT need bounding radius  as criterion so have all ghosts to correctly
  //NP calculate edge activation/deactivation in case edge of owned element is
  //NP completely off the processor

  template<int NUM_NODES>
  inline bool MultiNodeMeshParallel<NUM_NODES>::checkBorderElementLeft(int i,int dim,double lo, double hi)
  {
      /*NL*/// fprintf(this->screen,"this->center_(i)[dim] %f, rbound is %f, lo %f hi %f\n",this->center_(i)[dim],this->rBound_(i),lo,hi);
      if (this->center_(i)[dim] >= lo                     &&
          this->center_(i)[dim] <= hi + this->rBound_(i)     )
        return true;

      return false;

      /*NP  does not work - dont know why...
      bool all_left = true, all_right = true;

      for(int j = 0; j < NUM_NODES; j++)
         if (this->node_(i)[j][dim] >= lo)
            all_left = false;

      for(int j = 0; j < NUM_NODES; j++)
         if (this->node_(i)[j][dim] <= hi)
            all_right = false;

      if(all_left || all_right) return false;
      return true;*/
  }

  template<int NUM_NODES>
  inline bool MultiNodeMeshParallel<NUM_NODES>::checkBorderElementRight(int i,int dim,double lo, double hi)
  {
      if (this->center_(i)[dim] >= lo - this->rBound_(i) &&
          this->center_(i)[dim] <= hi                       )
        return true;

      return false;

      /*NP  does not work - dont know why...
      bool all_left = true, all_right = true;

      for(int j = 0; j < NUM_NODES; j++)
         if (this->node_(i)[j][dim] >= lo)
            all_left = false;

      for(int j = 0; j < NUM_NODES; j++)
         if (this->node_(i)[j][dim] <= hi)
            all_right = false;

      if(all_left || all_right) return false;
      return true;*/
  }

  /* ----------------------------------------------------------------------
   communicate properties to ghost elements
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::forwardComm()
  {
      int n,iElem;
      MPI_Request request;
      MPI_Status status;
      int me = this->comm->me;

      //NP get info about mesh movement to know what to communicate
      bool scale = this->isScaling();
      bool translate = this->isTranslating();
      bool rotate = this->isRotating();

      /*NL*/// fprintf(this->screen,"mesh %s: size_forward_ %d\n",this->mesh_id_,size_forward_);

      // exit here if no forward communication at all
      //NP done for efficiency reasons
      if(size_forward_ == 0)
        return;

      //NP THIS IS OLD, now using firstrecv_ as in Comm
      //NP clear ghost data, communication will recreate up-to-date data
      //NP requires ghosts to be on end of the container list
      //NP need to know which properties to delete, therefore need flags
      /*NL*/ //clearGhostForward(scale,translate,rotate);

      // exchange data with another proc
      // if other proc is self, just copy

      for (int iswap = 0; iswap < nswap_; iswap++)
      {
          if (sendproc_[iswap] != me)
          {
                if (size_forward_recv_[iswap])
                    MPI_Irecv(buf_recv_,size_forward_recv_[iswap],MPI_DOUBLE,recvproc_[iswap],0,this->world,&request);

                n = pushElemListToBuffer(sendnum_[iswap],sendlist_[iswap],buf_send_,OPERATION_COMM_FORWARD,scale,translate,rotate);
                /*NL*/// fprintf(this->screen,"proc %d: sendnum_[iswap] %d, recvnum_[iswap] %d,size_forward_recv_[iswap] %d\n",this->comm->me,sendnum_[iswap],recvnum_[iswap],size_forward_recv_[iswap]);
                /*NL*/// fprintf(this->screen,"proc %d: pushing n %d doubles to fwd buff, sendnum_[iswap] %d\n",this->comm->me,n,sendnum_[iswap]);

                if (n)
                    MPI_Send(buf_send_,n,MPI_DOUBLE,sendproc_[iswap],0,this->world);

                if (size_forward_recv_[iswap])
                    MPI_Wait(&request,&status);

                n = popElemListFromBuffer(firstrecv_[iswap],recvnum_[iswap],buf_recv_,OPERATION_COMM_FORWARD,scale,translate,rotate);
                /*NL*/// fprintf(this->screen,"proc %d: popping n %d doubles from fwd buff, recvnum_[iswap] %d\n",this->comm->me,n,recvnum_[iswap]);
          }
          else
          {
              n = pushElemListToBuffer(sendnum_[iswap],sendlist_[iswap],buf_send_,OPERATION_COMM_FORWARD,scale,translate,rotate);

              n = popElemListFromBuffer(firstrecv_[iswap],recvnum_[iswap],buf_recv_,OPERATION_COMM_FORWARD,scale,translate,rotate);
          }
      }
  }

  /* ----------------------------------------------------------------------
   reverse communication of properties on atoms every timestep
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::reverseComm()
  {
      int n;
      MPI_Request request;
      MPI_Status status;
      double *buf;
      int me = this->comm->me;

      //NP get info about mesh movement to know what to communicate
      bool scale = this->isScaling();
      bool translate = this->isTranslating();
      bool rotate = this->isRotating();

      // exchange data with another proc
      // if other proc is self, just copy

      for (int iswap = nswap_-1; iswap >= 0; iswap--)
      {
          if (sendproc_[iswap] != me)
          {
              if (size_reverse_recv_[iswap])
                  MPI_Irecv(buf_recv_,size_reverse_recv_[iswap],MPI_DOUBLE,sendproc_[iswap],0,this->world,&request);

              //NP n = avec->pack_reverse(recvnum[iswap],firstrecv[iswap],buf_send);
              n = pushElemListToBufferReverse(firstrecv_[iswap],recvnum_[iswap],buf_send_,OPERATION_COMM_REVERSE,scale,translate,rotate);

              if (n) MPI_Send(buf_send_,n,MPI_DOUBLE,recvproc_[iswap],0,this->world);
              if (size_reverse_recv_[iswap]) MPI_Wait(&request,&status);

              //NP avec->unpack_reverse(sendnum[iswap],sendlist[iswap],buf_recv);
              n = popElemListFromBufferReverse(sendnum_[iswap],sendlist_[iswap],buf_recv_,OPERATION_COMM_REVERSE,scale,translate,rotate);
          }
          else
          {
              n = pushElemListToBufferReverse(firstrecv_[iswap],recvnum_[iswap],buf_send_,OPERATION_COMM_REVERSE,scale,translate,rotate);
              n = popElemListFromBufferReverse(sendnum_[iswap],sendlist_[iswap],buf_recv_,OPERATION_COMM_REVERSE,scale,translate,rotate);
          }
      }
  }

#endif
