/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2018-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "mpi_liggghts.h"

/* ---------------------------------------------------------------------- */
// register custom MPI operations
/* ---------------------------------------------------------------------- */

MPI_Op MPI_ABSMIN_OP;
MPI_Op MPI_ABSMAX_OP;
namespace LAMMPS_NS
{

template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
void abs_min(void *invec, void *inoutvec, int len)
{
  T* in = (T*)invec;
  T* io = (T*)inoutvec;
  for (int i = 0; i < len; ++i) {
    const T x = std::abs(in[i]);
    const T y = std::abs(io[i]);
    io[i] = std::min(x,y);
  }
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
void abs_min(void *invec, void *inoutvec, int len)
{
  T* in = (T*)invec;
  T* io = (T*)inoutvec;
  for (int i = 0; i < len; ++i) {
    const T x = std::fabs(in[i]);
    const T y = std::fabs(io[i]);
    io[i] = std::min(x,y);
  }
}

/* ----------------------------------------------------------------------
   MPI reduce operator that computes the minimum absolute value.
------------------------------------------------------------------------- */

void mpi_abs_min_op(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
  const int    count = *len;
  MPI_Datatype dt    = *datatype;
  if (dt == MPI_INT) {
    abs_min<int>(invec, inoutvec, count);
  } else if (dt == MPI_LONG) {
    abs_min<long>(invec, inoutvec, count);
  } else if (dt == MPI_LONG_LONG) {
    abs_min<long long>(invec, inoutvec, count);
  } else if (dt == MPI_FLOAT) {
    abs_min<float>(invec, inoutvec, count);
  } else if (dt == MPI_DOUBLE) {
    abs_min<double>(invec, inoutvec, count);
  } else {
    printf("unsupported data type\n");
  }
}


template<typename T, typename std::enable_if<std::is_integral<T>::value, int>::type = 0>
void abs_max(void *invec, void *inoutvec, int len)
{
  T* in = (T*)invec;
  T* io = (T*)inoutvec;
  for (int i = 0; i < len; ++i) {
    const T x = std::abs(in[i]);
    const T y = std::abs(io[i]);
    io[i] = std::max(x,y);
  }
}

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
void abs_max(void *invec, void *inoutvec, int len)
{
  T* in = (T*)invec;
  T* io = (T*)inoutvec;
  for (int i = 0; i < len; ++i) {
    const T x = std::fabs(in[i]);
    const T y = std::fabs(io[i]);
    io[i] = std::max(x,y);
  }
}

/* ----------------------------------------------------------------------
   MPI reduce operator that computes the maximum absolute value.
------------------------------------------------------------------------- */

void mpi_abs_max_op(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
{
  const int count = *len;
  MPI_Datatype dt = *datatype;

  if (dt == MPI_INT) {
    abs_max<int>(invec, inoutvec, count);
  } else if (dt == MPI_LONG) {
    abs_max<long>(invec, inoutvec, count);
  } else if (dt == MPI_LONG_LONG) {
    abs_max<long long>(invec, inoutvec, count);
  } else if (dt == MPI_FLOAT) {
    abs_max<float>(invec, inoutvec, count);
  } else if (dt == MPI_DOUBLE) {
    abs_max<double>(invec, inoutvec, count);
  } else {
    printf("unsupported data type\n");
  }
}

void mpi_create_custom_operations()
{
  MPI_Op_create(mpi_abs_min_op, 1, &MPI_ABSMIN_OP);
  MPI_Op_create(mpi_abs_max_op, 1, &MPI_ABSMAX_OP);
}

void mpi_free_custom_operations()
{
  MPI_Op_free(&MPI_ABSMIN_OP);
  MPI_Op_free(&MPI_ABSMAX_OP);
}

} // end namespace LAMMPS_NS

