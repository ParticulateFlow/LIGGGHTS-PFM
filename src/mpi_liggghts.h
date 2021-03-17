/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-2016 DCS Computing GmbH, Linz
   Copyright 2017-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_MPI_LIGGGHTS_H
#define LMP_MPI_LIGGGHTS_H

#include <mpi.h>
#include <stdio.h>
#include "lmptype.h"

/* ---------------------------------------------------------------------- */
// a poor man's inline MPI wrappers for LIGGGHTS
/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS
{
void mpi_create_custom_operations();
void mpi_free_custom_operations();

inline int mpi_tag_upper_bound(MPI_Comm comm)
{
    int flag, *mpi_tag_ub;
    MPI_Comm_get_attr(comm, MPI_TAG_UB, &mpi_tag_ub, &flag);
    return *mpi_tag_ub;
}

/* ----------------------------------------------------------------------
   Helper function to be able to templatize wrappers
------------------------------------------------------------------------- */

template<typename T>
inline MPI_Datatype mpi_type()
{
  printf("\n\n\n**************LIGGGHTS MPI: ILLEGAL CALL TO mpi_type()*************\n\n\n");
  return 0;
}

template<>
inline MPI_Datatype mpi_type<double>()
{
  return MPI_DOUBLE;
}

template<>
inline MPI_Datatype mpi_type<int>()
{
  return MPI_INT;
}

template<>
inline MPI_Datatype mpi_type<int64_t>()
{
  return MPI_LONG_LONG_INT;
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Vector(T* vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Sum_Scalar(T& scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_SUM, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Scalar(T scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Scalar(T& scalar, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, &scalar, 1, mpi_type<T>(), MPI_MAX, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Scalar(T scalar, T& scalar_all, MPI_Comm comm)
{
  MPI_Allreduce(&scalar, &scalar_all, 1, mpi_type<T>(), MPI_MAX, comm);
}

/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Max_Vector(T *vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_MAX, comm);
}


/* ---------------------------------------------------------------------- */

template<typename T>
inline void MPI_Min_Vector(T* vector, int len, MPI_Comm comm)
{
  MPI_Allreduce(MPI_IN_PLACE, vector, len, mpi_type<T>(), MPI_MIN, comm);
}

/* ---------------------------------------------------------------------- */

inline void MPI_Allgather_Sum_Scalar(int scalar,int &scalar_acc,MPI_Comm comm)
{
    int rank,size, *allg;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    allg = new int[size];

    MPI_Allgather(&scalar,1,MPI_INT,allg,1,MPI_INT,comm);

    scalar_acc = 0;
    for (int iproc = 1; iproc < rank; iproc++)
       scalar_acc = scalar_acc + allg[iproc-1];

    delete []allg;
}

/* ----------------------------------------------------------------------
   Gather vector data from all processors at proc 0
   returns allocated and populated array vector0 to caller
------------------------------------------------------------------------- */

template<typename T>
inline int MPI_Gather0_Vector(T *vector, int size ,T *&vector_0,MPI_Comm comm)
{
    int me,nprocs, *recvcnts, *displs;
    int size_0;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    recvcnts = new int[nprocs];
    displs = new int[nprocs];

    //NP gather recvcount for each processor
    //NP MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm)
    MPI_Allgather(&size,1,MPI_INT,recvcnts,1,MPI_INT,comm);

    size_0 = 0;
    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
    {
        size_0 += recvcnts[iproc-1];
        displs[iproc] = displs[iproc-1] + recvcnts[iproc-1];
    }
    size_0 += recvcnts[nprocs-1];


    //NP allocate
    if(me == 0)
        vector_0 = new T[size_0];
    else
        vector_0 = 0;

    //NP use MPI_Gatherv to gather vector data at proc 0
    //NP MPI_Gatherv(sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, root, comm)

    MPI_Gatherv(vector,size,mpi_type<T>(),vector_0, recvcnts, displs, mpi_type<T>(),0, comm);

    delete []recvcnts;
    delete []displs;

    return size_0;
}

/* ----------------------------------------------------------------------
   Allgather vector data from all processors
   returns allocated and populated array vector_all to caller
------------------------------------------------------------------------- */

template<typename T>
inline int MPI_Allgather_Vector(T *vector, int size ,T *&vector_all,MPI_Comm comm)
{
    int me,nprocs, *recvcnts, *displs;
    int size_all;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &me);

    recvcnts = new int[nprocs];
    displs = new int[nprocs];

    //NP gather recvcount for each processor
    MPI_Allgather(&size,1,MPI_INT,recvcnts,1,MPI_INT,comm);

    size_all = 0;
    displs[0] = 0;
    for (int iproc = 1; iproc < nprocs; iproc++)
    {
        size_all += recvcnts[iproc-1];
        displs[iproc] = displs[iproc-1] + recvcnts[iproc-1];
    }
    size_all += recvcnts[nprocs-1];

    //NP allocate
    vector_all = new T[size_all];

    //NP use MPI_Allgatherv to gather vector data at each proc
    MPI_Allgatherv(vector,size,mpi_type<T>(),vector_all, recvcnts, displs, mpi_type<T>(), comm);

    delete []recvcnts;
    delete []displs;

    return size_all;
}

} // end namespace LAMMPS_NS


#endif
