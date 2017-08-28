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

#ifndef LMP_MULTISPHERE_PARALLEL_I_H
#define LMP_MULTISPHERE_PARALLEL_I_H

/* ---------------------------------------------------------------------- */

inline int MultisphereParallel::pack_exchange_rigid(int i, double *buf)
{
  /*NL*///if (screen) fprintf(screen,"packing tag %d\n",local2global[i]);

  int m = 1; //NP # of elements in buffer
  double xbound[3];
  bool dummy = false;

  // calculate xbound in global coo sys
  MathExtraLiggghts::local_coosys_to_cartesian
  (
    xbound,xcm_to_xbound_(i),
    ex_space_(i),ey_space_(i),ez_space_(i)
  );

  vectorAdd3D(xcm_(i),xbound,xbound);

  // have to pack xbound first because exchange() tests against first 3 values in buffer
  vectorToBuf3D(xbound,buf,m);
  /*NL*/// if (screen) fprintf(screen,"proc %d pack xbound %f %f %f id %d at step %d\n",comm->me,xbound[0],xbound[1],xbound[2],id_(i),update->ntimestep);

  m += customValues_.pushElemToBuffer(i,&(buf[m]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);

  buf[0] = m;
  /*NL*/ //if (screen) for(int i = 0; i < m; i++) fprintf(screen,"pack buf value %d: %f\n",i,buf[i]);
  return m;
}

/* ---------------------------------------------------------------------- */

inline int MultisphereParallel::unpack_exchange_rigid(double *buf)
{
  double xbound[3];
  bool dummy = false;
  int m = 0;

  int nvalues = buf[m++];

  bufToVector3D(xbound,buf,m);
  m += customValues_.popElemFromBuffer(&(buf[m]),OPERATION_COMM_EXCHANGE,dummy,dummy,dummy);

  nbody_++;

  /*NL*/ //if (screen) {
  /*NL*/ // fprintf(screen,"proc %d unpack xbound %f %f %f id %d at step %d\n",comm->me,xbound[0],xbound[1],xbound[2],id_(nbody_-1),update->ntimestep);

  /*NL*/ // for(int i = 0; i < nvalues; i++) fprintf(screen,"unpack buf value %d: %f\n",i,buf[i]);
  /*NL*/ // fprintf(screen,"nvalues %d m %d lengths: xcm %d, vcm %d fcm %d id %d nbody_ %d\n",
  /*NL*/ //                 nvalues,m,xcm_.size(),vcm_.size(),fcm_.size(),id_.size(),nbody_);
  /*NL*/ // for (int i = 0; i < nbody_; i++) {
  /*NL*/ //  fprintf(screen,"proc %d: id_(%d) %d\n",comm->me,i,id_(i));
  /*NL*/ //  printVec3D(screen,"xcm_",xcm_(i));
  /*NL*/ //  printVec3D(screen,"vcm_",vcm_(i));
  /*NL*/ // }
  /*NL*/ //}

  return nvalues;
}

#endif
