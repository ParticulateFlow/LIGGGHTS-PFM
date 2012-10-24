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

#define DELTA 10000

#include "multisphere_parallel.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "memory.h"

//NP same as in comm
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 1000

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

MultisphereParallel::MultisphereParallel(LAMMPS *lmp) :
  Multisphere(lmp),

  // initialize comm buffers & exchange memory
  maxsend_(BUFMIN),
  maxrecv_(BUFMIN),
  buf_send_((double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"frm:buf_send_")),
  buf_recv_((double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"frm:buf_send_"))
{

}

MultisphereParallel::~MultisphereParallel()
{
    memory->sfree(buf_send_);
    memory->sfree(buf_recv_);
}

/* ----------------------------------------------------------------------
   realloc the size of the send buffer as needed with BUFFACTOR & BUFEXTRA
   if flag = 1, realloc
   if flag = 0, don't need to realloc with copy, just free/malloc
------------------------------------------------------------------------- */

//NP from Comm::grow_send
void MultisphereParallel::grow_send(int n, int flag)
{
  maxsend_ = static_cast<int> (BUFFACTOR * n);
  if (flag)
    buf_send_ = (double *) memory->srealloc(buf_send_,(maxsend_+BUFEXTRA)*sizeof(double),"comm:buf_send_");
  else {
    memory->sfree(buf_send_);
    buf_send_ = (double *) memory->smalloc((maxsend_+BUFEXTRA)*sizeof(double),"comm:buf_send_");
  }
}

/* ----------------------------------------------------------------------
   free/malloc the size of the recv buffer as needed with BUFFACTOR
------------------------------------------------------------------------- */
//NP from Comm::grow_recv
void MultisphereParallel::grow_recv(int n)
{
  maxrecv_ = static_cast<int> (BUFFACTOR * n);
  memory->sfree(buf_recv_);
  buf_recv_ = (double *) memory->smalloc(maxrecv_*sizeof(double),        "comm:buf_recv_");
}

/* ----------------------------------------------------------------------
   exchange bodies with neighbor procs
------------------------------------------------------------------------- */

//NP from void Comm::exchange()
//NP assign bodies to procs based on x_bound, not xcm
//NP do this because this will lead to smaller bounding sphere radius
//NP and thus to a smaller corona extension

void MultisphereParallel::exchange()
{
  int i,m,nsend,nrecv,nrecv1,nrecv2;
  double lo,hi,value;
  double x[3];
  double *sublo,*subhi,*buf;
  MPI_Request request;
  MPI_Status status;

  // subbox bounds for orthogonal
  // triclinic not implemented

  sublo = domain->sublo;
  subhi = domain->subhi;

  // loop over dimensions

  /*NL*/// fprintf(screen,"step %d: proc %d has %d bodies\n",update->ntimestep,comm->me,nbody);
  /*NL*/// if(nbody) fprintf(screen,"step %d: proc %d body 0: xcv %f vcm %f\n",update->ntimestep,comm->me,xcm[0][2],vcm[0][2]);

  for (int dim = 0; dim < 3; dim++) {

    // fill buffer with atoms leaving my box, using < and >=
    // when atom is deleted, fill it in with last atom

    lo = sublo[dim];
    hi = subhi[dim];
    i = nsend = 0;

    while (i < nbody_) {

          MathExtraLiggghts::local_coosys_to_cartesian(x,xcm_to_xbound_(i),ex_space_(i),ey_space_(i),ez_space_(i));
          vectorAdd3D(xcm_(i),x,x);

          /*NL*///fprintf(screen,"step %d proc %d x %f %f %f, lo %f, hi %f\n",update->ntimestep,comm->me,x[0],x[1],x[2],lo,hi);

          if (x[dim] < lo || x[dim] >= hi)
          {
            if (nsend > maxsend_)
                grow_send(nsend,1);
            nsend += pack_exchange_rigid(i,&buf_send_[nsend]);
            /*NL*/// fprintf(screen,"lengths: xcm %d, vcm %d fcm %d id %d nbody_ %d\n",
            /*NL*///                  xcm_.size(),vcm_.size(),fcm_.size(),id_.size(),nbody_);
            remove_body(i);
            /*NL*/// fprintf(screen,"lengths: xcm %d, vcm %d fcm %d id %d nbody_ %d\n",
            /*NL*///                  xcm_.size(),vcm_.size(),fcm_.size(),id_.size(),nbody_);
          }
          else i++;
    }

    // send/recv atoms in both directions
    // if 1 proc in dimension, no send/recv, set recv buf to send buf
    // if 2 procs in dimension, single send/recv
    // if more than 2 procs in dimension, send/recv to both neighbors

    int procneigh[3][2];
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 2; j++)
            procneigh[i][j] = comm->procneigh[i][j];

    int *procgrid = comm->procgrid;

    if (procgrid[dim] == 1) {
      nrecv = nsend;
      buf = buf_send_;

    }
    else
    {
          MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][0],0,&nrecv1,1,MPI_INT,procneigh[dim][1],0,world,&status);
          nrecv = nrecv1;
          if (procgrid[dim] > 2) {
             MPI_Sendrecv(&nsend,1,MPI_INT,procneigh[dim][1],0,&nrecv2,1,MPI_INT,procneigh[dim][0],0,world,&status);
             nrecv += nrecv2;
          }

          if (nrecv > maxrecv_) grow_recv(nrecv);

          MPI_Irecv(buf_recv_,nrecv1,MPI_DOUBLE,procneigh[dim][1],0,world,&request);
          MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][0],0,world);
          MPI_Wait(&request,&status);

          if (procgrid[dim] > 2) {
            MPI_Irecv(&buf_recv_[nrecv1],nrecv2,MPI_DOUBLE,procneigh[dim][0],0,world,&request);
            MPI_Send(buf_send_,nsend,MPI_DOUBLE,procneigh[dim][1],0,world);
            MPI_Wait(&request,&status);
          }

          buf = buf_recv_;
    }

    // check incoming atoms to see if they are in my box
    // if so, add to my list

    m = 0;

    while (m < nrecv) {
      value = buf[m+dim+1];
      /*NL*///fprintf(screen,"value %f, lo %f hi %f\n",value,lo,hi);
      //NP unpack data if body in my subbox,  also adds the body
      if (value >= lo && value < hi) m += unpack_exchange_rigid(&buf[m]);
      else m += static_cast<int> (buf[m]);
    }
  }
}
