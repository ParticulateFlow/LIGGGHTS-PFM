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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */
#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

#include <string>

namespace LIGGGHTS {
namespace ContactModels {

// data available in noCollision() and collision()

struct ContactData {
  double radi;
  double radj;
  double radsum;
  double rsq; // squared distance between centers of collision partners
  double delta[3]; // line of centers

  double area_ratio;

  int * touch;
  double * contact_history;

  int i; // local particle index
  int j;

  bool is_wall;
  bool has_force_update;

  bool is_non_spherical; // superquadric

#ifdef NONSPHERICAL_ACTIVE_FLAG
  double contact_point[3];
#endif

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  double reff;
#endif

  ContactData() :
    radi(0.0),
    radj(0.0),
    radsum(0.0),
    rsq(0.0),
    area_ratio(1.0),
    touch(NULL),
    contact_history(NULL),
    i(-1),
    j(-1),
    is_wall(false),
    has_force_update(false),
    is_non_spherical(false)
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    ,reff(0.0)
#endif
  {}
};

// data available in collision() only

struct CollisionData: ContactData {
  double r; // distance between centers of collision partners
  double rinv; // one over r
  double en[3]; // direction of r
  double * v_i;
  double * v_j;
  double * omega_i;
  double * omega_j;

  double kt; // tangential stiffness
  double kn; // normal stiffness
  double gammat; // tangential damping
  double gamman; // normal damping

  double Fn; // normal force
  double Ft; // tangential force

  double vn; // relative normal velocity
  double deltan; // normal overlap
  double cri; // contact radius i
  double crj; // contact radius j
  double wr1;
  double wr2;
  double wr3;

  double vtr1; // relative tangential velocity (including rotational velocity)
  double vtr2;
  double vtr3;

  double mi; // particle mass
  double mj;
  double meff; // effective mass

  int computeflag;
  int shearupdate;
  int itype; // atom/material type
  int jtype;

  CollisionData() : Fn(0.0), Ft(0.0) {}
};

struct ForceData {
  double delta_F[3];       // total force acting on particle
  double delta_torque[3];  // torque acting on a particle

  ForceData()
  {
    reset();
  }

  inline void reset() {
    delta_F[0] = 0.0;
    delta_F[1] = 0.0;
    delta_F[2] = 0.0;
    delta_torque[0] = 0.0;
    delta_torque[1] = 0.0;
    delta_torque[2] = 0.0;
  }
};

struct ForceUpdate {
  double * f;
  double * torque;
  ForceData forces;

  ForceUpdate(double * f, double * torque, ForceData & forces) : f(f), torque(torque), forces(forces)
  {
  }
};

}

class IContactHistorySetup {
public:
  virtual int add_history_value(const std::string & name, const std::string & newtonflag) = 0;
};

}


#endif /* CONTACT_INTERFACE_H_ */
