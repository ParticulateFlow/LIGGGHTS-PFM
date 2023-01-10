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

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(mpi,CfdDatacouplingMPI)

#else

#ifndef LMP_CFD_DATACOUPLING_MPI_H
#define LMP_CFD_DATACOUPLING_MPI_H

#define MULTI_PARTITION_CFD

#include "cfd_datacoupling.h"
#include "multisphere_parallel.h"
#include "error.h"
#include "properties.h"
#include "universe.h"
#include <mpi.h>

namespace LAMMPS_NS {

class CfdDatacouplingMPI : public CfdDatacoupling {
 public:
  CfdDatacouplingMPI(class LAMMPS *, int,int, char **,class FixCfdCoupling*);
  ~CfdDatacouplingMPI();

  void exchange();

  virtual void pull(const char *name, const char *type, void *&ptr, const char *datatype,int iworld=0);
  virtual void push(const char *name, const char *type, void *&ptr, const char *datatype,int iworld=0);

  template <typename T> void pull_mpi(const char *,const char *,void *&,int iworld=0);
  template <typename T> void push_mpi(const char *,const char *,void *&,int iworld=0);

  virtual bool error_push()
  { return false;}

  virtual void allocate_external(int    **&data, int len2,int len1,     int initvalue);
  virtual void allocate_external(int    **&data, int len2,const char *keyword,int initvalue);
  virtual void allocate_external(double **&data, int len2,int len1,     double initvalue);
  virtual void allocate_external(double **&data, int len2,const char *keyword,double initvalue);

 private:
  template <typename T> T* check_grow(int len);
  template <typename T> MPI_Datatype mpi_type_dc();

  // 1D helper array needed to allreduce the quantities
  int len_allred_double;
  double *allred_double;

  int len_allred_int;
  int *allred_int;
};

/* ---------------------------------------------------------------------- */
//NP OF to LIGGGHTS, global to local
//NP ptr to is local, ptr from is global

template <typename T>
void CfdDatacouplingMPI::pull_mpi(const char *name,const char *type,void *&from,int iworld)
{
    int len1 = -1, len2 = -1, m;
    // len1 = atom->tag_max(); except for scalar-global, vector-global, matrix-global

    // get reference where to write the data
    void * to = find_pull_property(name,type,len1,len2);

#ifdef MULTI_PARTITION_CFD
    if (universe->existflag == 0 || iworld == universe->iworld)
#endif
    {
    if (atom->nlocal && (!to || len1 < 0 || len2 < 0))
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->one(FLERR,"This is fatal");
    }
    }

    int total_len = len1*len2;
#ifdef MULTI_PARTITION_CFD
    if (universe->existflag == 1) // enforce size of requested world
        MPI_Bcast(&total_len, 1, MPI_INT, universe->root_proc[iworld], universe->uworld);
#endif

    // return if no data to transmit
    if(total_len < 1) return;

    // check memory allocation
    T* allred = check_grow<T>(total_len);

    // perform allreduce on incoming data
    T **from_t = (T**)from;
#ifndef MULTI_PARTITION_CFD
    MPI_Allreduce(&(from_t[0][0]), &(allred[0]), total_len, mpi_type_dc<T>(), MPI_SUM, world);
#else
    //MPI_Reduce + MPI_Bcast
    // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    // int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    //reduce from all procs in universe onto root proc of iworld
    MPI_Reduce   (&(from_t[0][0]), &(allred[0]), total_len, mpi_type_dc<T>(), MPI_SUM, universe->root_proc[iworld], universe->uworld);
    //bcast to iworld
    if (iworld == universe->iworld) {
        MPI_Bcast    (&(allred[0]), total_len, mpi_type_dc<T>(), 0, world);
    }
#endif

#ifdef MULTI_PARTITION_CFD
    if(iworld == universe->iworld) // only copy to requested world
#endif
    {
    // copy data - loops over max # global atoms, bodies
    if(strcmp(type,"scalar-atom") == 0)
    {
        T *to_t = (T*) to;
        for (int i = 0; i < len1; i++)
            if ((m = atom->map(i+1)) >= 0)
                to_t[m] = allred[i];
    }
    else if(strcmp(type,"vector-atom") == 0)
    {
        T **to_t = (T**) to;
        for (int i = 0; i < len1; i++)
            if ((m = atom->map(i+1)) >= 0)
                for (int j = 0; j < len2; j++)
                    to_t[m][j] = allred[i*len2 + j];
    }
    else if(strcmp(type,"scalar-multisphere") == 0)
    {
        T *to_t = (T*) to;
        MultisphereParallel *ms_data = properties_->ms_data();
        if(!ms_data)
            error->one(FLERR,"Transferring a multisphere property from/to LIGGGHTS requires a fix multisphere");
        for (int i = 0; i < len1; i++)
            if ((m = ms_data->map(i+1)) >= 0)
                to_t[m] = allred[i];
    }
    else if(strcmp(type,"vector-multisphere") == 0)
    {
        T **to_t = (T**) to;
        MultisphereParallel *ms_data = properties_->ms_data();
        if(!ms_data)
            error->one(FLERR,"Transferring a multisphere property from/to LIGGGHTS requires a fix multisphere");
        for (int i = 0; i < len1; i++)
            if ((m = ms_data->map(i+1)) >= 0)
                for (int j = 0; j < len2; j++)
                    to_t[m][j] = allred[i*len2 + j];
    }
    else if(strcmp(type,"scalar-global") == 0 || strcmp(type,"vector-global") == 0)
    {
        T *to_t = (T*) to;
        for (int i = 0; i < len1; i++)
            to_t[i] = allred[i];
    }
    else if(strcmp(type,"matrix-global") == 0)
    {
        T **to_t = (T**) to;
        for (int i = 0; i < len1; i++)
            for (int j = 0; j < len2; j++)
                to_t[i][j] = allred[i*len2 + j];
    }
    else error->one(FLERR,"Illegal data type in CfdDatacouplingMPI::pull");
    }
}

/* ---------------------------------------------------------------------- */
//NP LIGGGHTS to OF, local to global
//NP ptr from is local, ptr to is global
//NP len1 is global # of datums (max tag)

template <typename T>
void CfdDatacouplingMPI::push_mpi(const char *name,const char *type,void *&to,int iworld)
{
    int len1 = -1, len2 = -1, id;

    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int nbodies = 0;

    MultisphereParallel *ms_data = properties_->ms_data();
    if(ms_data) nbodies = ms_data->n_body();

    // get reference where to write the data
    void * from = find_push_property(name,type,len1,len2);

#ifdef MULTI_PARTITION_CFD
    if (universe->existflag == 0 || iworld == universe->iworld)
#endif
    if (atom->nlocal && (!from || len1 < 0 || len2 < 0))
    {
        /*NL*/ //if(screen) fprintf(screen,"nlocal %d, len1 %d lens2 %d\n",atom->nlocal,len1,len2);
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->one(FLERR,"This is fatal");
    }

    int total_len = len1*len2;
#ifdef MULTI_PARTITION_CFD
    if (universe->existflag == 1) // enforce size of requested world
        MPI_Bcast(&total_len, 1, MPI_INT, universe->root_proc[iworld], universe->uworld);
#endif
    // return if no data to transmit
    if(total_len < 1) return;

    // check memory allocation
    T * allred = check_grow<T>(total_len);

    // zeroize before using allreduce
    vectorZeroizeN(allred,total_len);

#ifdef MULTI_PARTITION_CFD
    if(iworld == universe->iworld) // only copy from requested world
#endif
    {
    // copy data - loop local # atoms, bodies
    if(strcmp(type,"scalar-atom") == 0)
    {
        T *from_t = (T*) from;
        for (int i = 0; i < nlocal; i++)
        {
            id = tag[i];
            allred[id-1] = from_t[i];
        }
    }
    else if(strcmp(type,"vector-atom") == 0)
    {
        T **from_t = (T**) from;
        for (int i = 0; i < nlocal; i++)
        {
            id = tag[i];
            for (int j = 0; j < len2; j++)
                allred[(id-1)*len2 + j] = from_t[i][j];
        }
    }
    else if(strcmp(type,"scalar-multisphere") == 0)
    {
        T *from_t = (T*) from;
        if(!ms_data)
            error->one(FLERR,"Transferring a multisphere property from/to LIGGGHTS requires a fix multisphere");
        for (int i = 0; i < nbodies; i++) // loops over # local bodies
        {
            id = ms_data->tag(i);
            allred[id-1] = from_t[i];
        }
    }
    else if(strcmp(type,"vector-multisphere") == 0)
    {
        T **from_t = (T**) from;
        if(!ms_data)
            error->one(FLERR,"Transferring a multisphere property from/to LIGGGHTS requires a fix multisphere");
        for (int i = 0; i < nbodies; i++) // loops over # local bodies
        {
            id = ms_data->tag(i);
            for (int j = 0; j < len2; j++)
            {
                allred[(id-1)*len2 + j] = from_t[i][j];
                /*NL*/ //if(screen) fprintf(screen,"id  %d from %d %d: %f\n",id ,i,j,from_t[i][j]);
            }
        }
        /*NL*/ //if (screen) printVecN(screen,"allred",allred,total_len);
    }
    else if(strcmp(type,"scalar-global") == 0 || strcmp(type,"vector-global") == 0 || strcmp(type,"matrix-global") == 0)
    {
        T **from_t = (T**) from;
        for (int i = 0; i < len1; i++)
            for (int j = 0; j < len2; j++)
                allred[i*len2 + j] = from_t[i][j];
    }
    else error->one(FLERR,"Illegal data type in CfdDatacouplingMPI::pull");
    }

    // perform allreduce on outgoing data
    T **to_t = (T**)to;
#ifndef MULTI_PARTITION_CFD
    MPI_Allreduce(&(allred[0]),&(to_t[0][0]),total_len,mpi_type_dc<T>(),MPI_SUM,world);
#else
    //MPI_Reduce + MPI_Bcast
    // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
    // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    // reduce from request world onto root of requested world
    if (iworld == universe->iworld) {
        MPI_Reduce(&(allred[0]), &(to_t[0][0]), total_len, mpi_type_dc<T>(), MPI_SUM, 0, world);
    }
    // bcast from root of requested world to all procs in universe
    MPI_Bcast(&(to_t[0][0]), total_len, mpi_type_dc<T>(), universe->root_proc[iworld], universe->uworld);
#endif

}

/* ---------------------------------------------------------------------- */

template<typename T>
inline MPI_Datatype CfdDatacouplingMPI::mpi_type_dc()
{
    error->all(FLERR,"Illegal call to mpi_type_dc(), valid types are int and double");
    return 0;
}

template<>
inline MPI_Datatype CfdDatacouplingMPI::mpi_type_dc<double>()
{
    return MPI_DOUBLE;
}

template<>
inline MPI_Datatype CfdDatacouplingMPI::mpi_type_dc<int>()
{
    return MPI_INT;
}

/* ---------------------------------------------------------------------- */

template <typename T>
T* CfdDatacouplingMPI::check_grow(int len)
{
    error->all(FLERR,"Illegal call to template <typename T> T* check_grow(int len)");
    return NULL;
}

template <>
inline double* CfdDatacouplingMPI::check_grow<double>(int len)
{
    while(len > len_allred_double)
        len_allred_double += 10000;

    allred_double = (double*) memory->srealloc(allred_double,len_allred_double*sizeof(double),"fix_cfd_coupling:allred_double");
    vectorZeroizeN(allred_double,len_allred_double);
    return allred_double;
}

template <>
inline int* CfdDatacouplingMPI::check_grow<int>(int len)
{
    while(len > len_allred_int)
        len_allred_int += 10000;

    allred_int = (int*) memory->srealloc(allred_int,len_allred_int*sizeof(int),"fix_cfd_coupling:allred_int");
    vectorZeroizeN(allred_int,len_allred_int);
    return allred_int;
}

}

#endif
#endif
