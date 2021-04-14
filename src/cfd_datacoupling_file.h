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

   CfdDataCouplingStyle(file,CfdDatacouplingFile)

#else

#ifndef LMP_CFD_DATACOUPLING_FILE_H
#define LMP_CFD_DATACOUPLING_FILE_H

#include "cfd_datacoupling.h"

namespace LAMMPS_NS {

class CfdDatacouplingFile : public CfdDatacoupling {
 public:
  CfdDatacouplingFile(class LAMMPS *, int, int, char **,class FixCfdCoupling* fc);
  ~CfdDatacouplingFile();
  friend class FixTempFromFile;

  void pull(const char *, const char *, void *&, const char *, int iworld = 0);
  void push(const char *, const char *, void *&, const char *, int iworld = 0);
  virtual void post_create();

  void exchange();


  private:
   char* filepath;
   //NP if set, expect 0 or 1 to be appended to file path
   //NP setting to 0 does also deactivate op_complete()
   int append;

   bool firstexec;
   int t0;

   char * getFilePath(const char *name, bool flag);
   void op_complete(const char *name);
   void writeVectorData(const char *name,  double ** field);
   void writeScalarData(const char *name,  double * field);
   void writeGlobalVectorData(const char *name,  double * field,int);
   void writeGlobalArrayData(const char *name,  double ** field,int,int);

   void readVectorData(const char *name,  double ** field);
   void readScalarData(const char *name,  double * field);
   void readGlobalVectorData(const char* name, double *field, int &len);
   void readGlobalArrayData(const char *name, double ** field, int &len1, int &len2);

   //NP would need clump-type things here
};

}

#endif
#endif
