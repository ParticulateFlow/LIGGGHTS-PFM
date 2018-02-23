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

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include <mpi.h>

double** o2o_liggghts_get_boundingbox(void *ptr);

int o2o_liggghts_get_maxtag(void *ptr);
/*
int o2o_liggghts_get_maxtag_ms(void *ptr);
int o2o_liggghts_get_ntypes_ms(void *ptr);
double* o2o_liggghts_get_vclump_ms(void *ptr);
*/
void* o2o_locate_coupling_fix(void *ptr);
void o2o_data_of_to_liggghts
(
    const char *name,
    const char *type,
    void *ptr,
    void *data,
    const char *datatype,
    const int* ids,
    const int ncollected
);
void o2o_update_rm(void *ptr);
void o2o_check_datatransfer(void *ptr);

void o2o_allocate_external_int(int    **&data, int len2, int len1, int initvalue, void *ptr);
void o2o_allocate_external_int(int    **&data, int len2, const char *,   int initvalue, void *ptr);


void o2o_allocate_external_double(double **&data, int len2, int len1, double initvalue, void *ptr);
void o2o_allocate_external_double(double **&data, int len2, const char *,   double initvalue, void *ptr);