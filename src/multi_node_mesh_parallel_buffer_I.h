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

#ifndef LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H
#define LMP_MULTI_NODE_MESH_PARALLEL_BUFFER_I_H

  /* ----------------------------------------------------------------------
   push / pop for exchange
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushExchange(int dim)
  {
      //NP do NOT make a local copy of buf_send_ as this fct calls re-allocation!!

      // scale translate rotate not needed here
      bool dummy = false;
      double checklo,checkhi;

      //NP consistent with Domain::is_in_subdomain()
      checklo = this->domain->sublo[dim];
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

      int nsend = 0, nsend_this = 0;
      int i = 0;
      while(i < nLocal_)
      {
          if(!(this->center_(i)[dim] >= checklo && this->center_(i)[dim] < checkhi))
          {
              nsend_this = pushElemToBuffer(i,&(buf_send_[nsend+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
              buf_send_[nsend] = static_cast<double>(nsend_this+1);
              nsend += (nsend_this+1);
              /*NL*/// fprintf(this->screen,"step %d: proc %d pushes element id %d with center %f %f %f\n",
              /*NL*///         this->update->ntimestep,this->comm->me,this->id_slow(i),this->center_(i)[0],this->center_(i)[1],this->center_(i)[2]);
              if (nsend > maxsend_)
                  grow_send(nsend,1);
              this->deleteElement(i); //NP deleteElement() decreases nLocal
          }
          else i++;
      }
      return nsend;
  }

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::popExchange(int nrecv,int dim,double *buf)
  {
      double center_elem[3];
      double checklo,checkhi;
      int m = 0, nrecv_this;

      // scale translate rotate not needed here
      bool dummy = false;

      //NP consistent with Domain::is_in_subdomain()
      checklo = this->domain->sublo[dim];
      if(this->domain->subhi[dim] == this->domain->boxhi[dim])
        checkhi = this->domain->boxhi[dim] + SMALL_DMBRDR;
      else
        checkhi = this->domain->subhi[dim];

      while (m < nrecv)
      {
          // number of values is first in buffer
          nrecv_this = static_cast<int>(buf[m]);

          // center is next in buffer, test it
          vectorCopy3D(&(buf[m+1]),center_elem);

          /*NL*/// fprintf(this->screen,"proc %d: nrecv_this %d center %f %f %f xlo %f xhi %f\n",
          /*NL*///         this->comm->me,nrecv_this,center_elem[0],center_elem[1],center_elem[2],this->domain->sublo[0],this->domain->subhi[0]);

          //NP do not ask if the center is completely in, just ask for center_elem[dim]
          //NP this makes elements go around the corner
          if(center_elem[dim] >= checklo && center_elem[dim] < checkhi)
          {
            popElemFromBuffer(&(buf[m+1]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);
            nLocal_++;
            /*NL*/// fprintf(this->screen,"proc %d pops element id %d with center %f %f %f\n",
            /*NL*///                      this->comm->me,this->id_slow(nLocal_-1),center_elem[0],center_elem[1],center_elem[2]);
          }
          /*NL*/// else
          /*NL*///   fprintf(this->screen,"proc %d DID NOT pop one with center %f %f %f\n",
          /*NL*///                      this->comm->me,center_elem[0],center_elem[1],center_elem[2]);

          /*NL*/// this->error->one(FLERR,"end");

          m += nrecv_this;
      }
  }

  /* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::writeRestart(FILE *fp)
  {
      int size_this;

      // # elements
      int nlocal = this->sizeLocal();
      int nglobal = sizeGlobal();

      // buffer sizes
      int sizeMesh, sizeElements, sizeElements_all;

      sizeMesh = sizeRestartMesh();
      sizeElements = nlocal * (sizeRestartElement() + 1); //NP need 1 extra space for pre-element buffer length

      double *bufMesh = NULL, *sendbufElems = NULL, *recvbufElems = NULL;
      bool dummy = false;

      /*NL*/ //fprintf(this->screen,"buffersize mesh %d buffersize elems local %d\n",sizeMesh,sizeElements);

      // pack global data into buffer
      // do this only on proc 0
      if(this->comm->me == 0)
      {
          this->memory->create(bufMesh,sizeMesh,"MultiNodeMeshParallel::writeRestart:bufMesh");
          pushMeshPropsToBuffer(bufMesh, OPERATION_RESTART,dummy,dummy,dummy);
      }

      // allocate send buffer and pack element data
      // all local elements are in list
      this->memory->create(sendbufElems,sizeElements,"MultiNodeMeshParallel::writeRestart:sendbuf");
      sizeElements = 0;
      for(int i = 0; i < nlocal; i++)
      {
          size_this = pushElemToBuffer(i,&(sendbufElems[sizeElements+1]),OPERATION_RESTART,dummy,dummy,dummy);
          sendbufElems[sizeElements] = static_cast<double>(size_this+1);
          sizeElements += (size_this+1);
      }

      // gather the pre-element data
      //NP send from all to proc 0
      sizeElements_all = MPI_Gather0_Vector(sendbufElems,sizeElements,recvbufElems,this->world);

      /*NL*/ //fprintf(this->screen,"buffersize mesh %d buffersize elems local %d, buffersize elems global %d\n",sizeMesh,sizeElements, sizeElements_all);

      // actually write data to restart file
      // do this only on proc 0
      if(this->comm->me == 0)
      {
        double nG = static_cast<double>(nglobal);

        // for error check
        double sE = static_cast<double>(sizeRestartElement());
        double sM = static_cast<double>(sizeRestartMesh());

        /*NL*///fprintf(this->screen,"nglobal %f, sE %f, sM %f\n",nG,sE,sM);

        // size with 3 extra values
        int size = (sizeMesh+sizeElements_all+3) * sizeof(double);

        // write size
        fwrite(&size,sizeof(int),1,fp);

        // write 3 extra values
        fwrite(&nG,sizeof(double),1,fp);
        fwrite(&sE,sizeof(double),1,fp);
        fwrite(&sM,sizeof(double),1,fp);

        // write per-element and mesh data
        fwrite(recvbufElems,sizeof(double),sizeElements_all,fp);
        fwrite(bufMesh,sizeof(double),sizeMesh,fp);
      }

      // free mem

      if(bufMesh)
        this->memory->destroy(bufMesh);

      this->memory->destroy(sendbufElems);

      //NP need to use simple delete [] b/c MPI_Gather0_Vector uses new
      if(recvbufElems)
        delete []recvbufElems;
  }

  /* ----------------------------------------------------------------------
   restart functionality - write all required data into restart buffer
   executed on all processes, but only proc 0 writes into writebuf
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  void MultiNodeMeshParallel<NUM_NODES>::restart(double *list)
  {
      int m, nglobal, nrecv_this, nO, sE, sM;
      bool dummy = false;

      m = 0;

      nglobal = static_cast<int> (list[m++]);
      sE = static_cast<int> (list[m++]);
      sM = static_cast<int> (list[m++]);

      /*NL*///fprintf(this->screen,"nglobal %d, sE %d, sM %d, sizeRestartElement() %d, sizeRestartMesh() %d\n",nglobal,sE,sM,sizeRestartElement(),sizeRestartMesh());

      if(sE != sizeRestartElement() || sM != sizeRestartMesh())
          this->error->all(FLERR,"Incompatible mesh restart file - mesh has different properties in restarted simulation");

      for(int i = 0; i < nglobal; i++)
      {
          nrecv_this = static_cast<int>(list[m]);
          /*NL*///fprintf(this->screen,"nrecv_this %d\n",nrecv_this);
          popElemFromBuffer(&(list[m+1]),OPERATION_RESTART,dummy,dummy,dummy);
          m += nrecv_this;
      }

      popMeshPropsFromBuffer(&list[m],OPERATION_RESTART,dummy,dummy,dummy);
  }

  /* ----------------------------------------------------------------------
   size of restart data
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartMesh()
  {
      //NP get info about mesh movement to know what to communicate
      bool dummy = false;
      return meshPropsBufSize(OPERATION_RESTART,dummy,dummy,dummy);
  }

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::sizeRestartElement()
  {
      //NP get info about mesh movement to know what to communicate
      bool dummy = false;
      return elemBufSize(OPERATION_RESTART,dummy,dummy,dummy);
  }

  /*NP ----------------------------------------------------------------------
   manual communication for this class, so lots of if / else statements here
   to evaluate operation flag

   need not do that for derived classes that implement their properties
   via TrackingMesh, for these cases the operation flag is evaluated
   in the container class

   on OPERATION_BORDERS, communicate center and rbound so dont have to
   refresh them
  ------------------------------------------------------------------------- */


  /* ----------------------------------------------------------------------
   return required buffer size for a list of elements for borders(),forwardComm()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
      return n*elemBufSize(operation,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   push a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemListToBuffer(int n, int *list, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      /*NL*/// fprintf(this->screen,"pushListToBuffer n %d\n",n);
      /*NL*/// for(int ii = 0; ii < n; ii++)
      /*NL*///     fprintf(this->screen,"  list[%d] %d\n",ii,list[ii]);

      int nsend = 0;

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          //NP push center first to test against
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemListToBuffer(n,list,&(buf[nsend]),operation);
          return nsend;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          /*NL*///fprintf(this->screen,"comm forward, translate is %s\n",translate?"yes":"no");
          //NP node_orig cannot change during a run
          //NP currently proc updates owned and ghost elements for moving mesh
          //NP if(translate || rotate || scale)
          //NP  nsend += MultiNodeMesh<NUM_NODES>::node_.pushListToBuffer(n,list,&(buf[nsend]),operation);
          return nsend;
      }

      //NP OPERATION_RESTART not implemented, is per-element operation, not a list operation
      //NP OPERATION_COMM_REVERSE is implemented in pushElemListToBufferReverse
      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for borders(), forwardComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemListFromBuffer(int first, int n, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nrecv = 0;

      //NP OPERATION_RESTART not implemented, is per-element operation, not a list operation

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemListFromBuffer(first,n,&(buf[nrecv]),operation);
          return nrecv;
      }

      if(operation == OPERATION_COMM_FORWARD)
      {
          //NP node_orig cannot change during a run
          //NP currently proc updates owned and ghost elements for moving mesh
          //NP if(translate || rotate || scale)
          //    nrecv += MultiNodeMesh<NUM_NODES>::node_.popListFromBuffer(first,n,&(buf[nrecv]),operation);
          return nrecv;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemListToBufferReverse(int first, int n, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nsend = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        //NP no reverse comm here
        return nsend;
      }

      //NP other stuff implemented in pushElemListToBuffer
      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop a list of elements for reverseComm()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemListFromBufferReverse(int n, int *list, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nrecv = 0;

      if(operation == OPERATION_COMM_REVERSE)
      {
        //NP no reverse comm here
        return nrecv;
      }

      //NP other stuff implemented in pushElemListToBuffer
      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   return required buffer size for one element for exchange()
   must match push / pop implementation
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
      int size_buf = 0;

      //NP need to implement for per-element and per-list calls here
      //NP since per-list = n * per_element

      if(operation == OPERATION_RESTART)
      {
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          return size_buf;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          size_buf += MultiNodeMesh<NUM_NODES>::center_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          size_buf += MultiNodeMesh<NUM_NODES>::rBound_.elemBufSize();
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            size_buf += MultiNodeMesh<NUM_NODES>::node_orig_->elemBufSize();
          return size_buf;
      }

      //NP OPERATION_COMM_FORWARD, OPERATION_COMM_REVSERSE are list operations, not per-element operations
      //NP need to implement this here since elemListBufSize() refers to here
      if(operation == OPERATION_COMM_FORWARD)
      {
          //NP node_orig cannot change during a run
          //NP if(translate || rotate || scale)
          //NP  size_buf += MultiNodeMesh<NUM_NODES>::node_.elemBufSize();
          return size_buf;
      }

      if(operation == OPERATION_COMM_REVERSE)
      {
        //NP nothing todo here
        //NP need to implement this
        return size_buf;
      }

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::elemBufSize");
      return 0;
  }

  /* ----------------------------------------------------------------------
   push one element for exchange()
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer(int i, double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nsend = 0;

      //NP OPERATION_COMM_REVERSE is a list operation, not per-element operation

      if(operation == OPERATION_RESTART)
      {
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);

          return nsend;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          //NP push center first to test against
          nsend += MultiNodeMesh<NUM_NODES>::center_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::node_.pushElemToBuffer(i,&(buf[nsend]),operation);
          nsend += MultiNodeMesh<NUM_NODES>::rBound_.pushElemToBuffer(i,&(buf[nsend]),operation);
          if(this->node_orig_)
              nsend += this->node_orig_->pushElemToBuffer(i,&(buf[nsend]),operation);
          return nsend;
      }

      //NP OPERATION_COMM_FORWARD, OPERATION_COMM_REVSERSE are list operations, not per-element operations

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::pushElemToBuffer");
      return 0;
  }

  /* ----------------------------------------------------------------------
   pop one element for exchange, restart or restart
   depending on operation and if mesh scales, translates or rotates,
   different properties are communicated
  ------------------------------------------------------------------------- */

  template<int NUM_NODES>
  int MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer(double *buf,int operation,bool scale,bool translate,bool rotate)
  {
      int nrecv = 0;

      //NP OPERATION_COMM_REVERSE is a list operation, not per-element operation

      if(operation == OPERATION_RESTART)
      {
          MultiVectorContainer<double,NUM_NODES,3> nodeTmp;
          //NP pop node info from buffer and call addElement so everything is initialized correctly
          nrecv += nodeTmp.popElemFromBuffer(&(buf[nrecv]),operation);
          this->addElement(nodeTmp.begin()[0],-1);

          //NP this->addElement() creates properties, calculates everything
          //NP also adds custom values for which restart data is available
          //NP these are now cleared and re-created with correct restart data
          //NP in TrackingMesh::popElemFromBuffer
          this->prop().deleteRestartElement(nLocal_-1,false,false,false);

          return nrecv;
      }

      if(operation == OPERATION_COMM_EXCHANGE || operation == OPERATION_COMM_BORDERS)
      {
          nrecv += MultiNodeMesh<NUM_NODES>::center_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::node_.popElemFromBuffer(&(buf[nrecv]),operation);
          nrecv += MultiNodeMesh<NUM_NODES>::rBound_.popElemFromBuffer(&(buf[nrecv]),operation);
          if(MultiNodeMesh<NUM_NODES>::node_orig_)
            nrecv += MultiNodeMesh<NUM_NODES>::node_orig_->popElemFromBuffer(&(buf[nrecv]),operation);
          return nrecv;
      }

      //NP OPERATION_COMM_FORWARD, OPERATION_COMM_REVSERSE are list operations, not per-element operations

      this->error->one(FLERR,"Illegal operation in MultiNodeMeshParallel<NUM_NODES>::popElemFromBuffer");
      return 0;
  }

#endif
