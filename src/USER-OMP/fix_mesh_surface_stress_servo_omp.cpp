/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
   Copyright 2013-     JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "fix_mesh_surface_stress_servo_omp.h"
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "group.h"
#include "force.h"
#include "bounding_box.h"
#include "input_mesh_tri.h"
#include "fix_contact_history.h"
#include "fix_neighlist_mesh_omp.h"
#include "multi_node_mesh.h"
#include "modify.h"
#include "comm.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixMeshSurfaceStressServoOMP::FixMeshSurfaceStressServoOMP(LAMMPS *lmp, int narg, char **arg)
: FixMeshSurfaceStressServo(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

FixMeshSurfaceStressServoOMP::~FixMeshSurfaceStressServoOMP()
{
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   called from fix wall/gran out of post_create()
------------------------------------------------------------------------- */

void FixMeshSurfaceStressServoOMP::createWallNeighList(int igrp)
{
  if(fix_mesh_neighlist_) return;
  char *neighlist_name = new char[strlen(id)+1+20];
  sprintf(neighlist_name,"wall_neighlist_%s",id);

  char **fixarg = new char*[4];
  fixarg[0]= neighlist_name;
  fixarg[1]= (char *) "all";
  fixarg[2]= (char *) "neighlist/mesh/omp";
  fixarg[3]= id;
  modify->add_fix(4,fixarg);

  fix_mesh_neighlist_ =
      static_cast<FixNeighlistMesh*>(modify->find_fix_id(neighlist_name));

  // fix added with "all", correct this now
  fix_mesh_neighlist_->igroup = igrp;
  fix_mesh_neighlist_->groupbit = group->bitmask[igrp];

  delete []fixarg;
  delete []neighlist_name;
}

/* ----------------------------------------------------------------------
   called from fix messflow/mesh out of post_create()
------------------------------------------------------------------------- */

class FixNeighlistMesh* FixMeshSurfaceStressServoOMP::createOtherNeighList(int igrp,const char *nId)
{
  FixNeighlistMesh* neighlist;

  char *neighlist_name = new char[strlen(id)+1+20+strlen(nId)+1];
  sprintf(neighlist_name,"neighlist_%s_%s",nId,id);

  char **fixarg = new char*[4];
  fixarg[0]= neighlist_name;
  fixarg[1]= (char *) "all";
  fixarg[2]= (char *) "neighlist/mesh/omp";
  fixarg[3]= id;
  modify->add_fix(4,fixarg);

  neighlist =
      static_cast<FixNeighlistMesh*>(modify->find_fix_id(neighlist_name));

  // fix added with "all", correct this now
  neighlist->igroup = igrp;
  neighlist->groupbit = group->bitmask[igrp];

  delete []fixarg;
  delete []neighlist_name;

  return neighlist;
}
