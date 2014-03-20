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
#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_HAMAKER,hamaker,3)
#else
#ifndef COHESION_MODEL_HAMAKER_H_
#define COHESION_MODEL_HAMAKER_H_
#include "contact_models.h"
#include "domain.h"

namespace ContactModels {
  using namespace std;
  using namespace LAMMPS_NS;
  
  template<typename Style>
  class CohesionModel<COHESION_HAMAKER, Style> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), aHamakerEff(NULL), hCutEff(NULL), hMaxEff(NULL)
    {
      if(domain->dimension!=3)
          error->all(FLERR,"Cohesion model valid for 3d simulations only");
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("aHamakerEff", &MODEL_PARAMS::createHamakerConstant);
      registry.registerProperty("hCutEff", &MODEL_PARAMS::createHamakerMinimumParticleDistance);
      registry.registerProperty("hMaxEff", &MODEL_PARAMS::createHamakerMaxEff);

      registry.connect("aHamakerEff", aHamakerEff);
      registry.connect("hCutEff", hCutEff);
      registry.connect("hMaxEff", hMaxEff);
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) //NP modified C.K.
    {
      // since the particles are in contact:
      //     * the maximum cohesive force acts between the particles
      //     * the calculation of the displacement (hIJ) is not required (no cdata.r)
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double reff = (ri*rj)/(ri+rj);
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      // since the particles are in contact, the maximum cohesive force acts
      // using the equation directly without function call may be faster, since you can skip one if-statement
      const double Fn_coh = calcCohesiveForce(0.0, reff, itype, jtype);

      cdata.Fn += Fn_coh;

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn_coh * cdata.area_ratio;
        i_forces.delta_F[0] += Fn_ * cdata.en[0];
        i_forces.delta_F[1] += Fn_ * cdata.en[1];
        i_forces.delta_F[2] += Fn_ * cdata.en[2];
      } else {
        const double fx = Fn_coh * cdata.en[0];
        const double fy = Fn_coh * cdata.en[1];
        const double fz = Fn_coh * cdata.en[2];

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
      }
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

    void noCollision(ContactData& cdata, ForceData& i_forces, ForceData& j_forces)
    {
      const double r = sqrt(cdata.rsq);
      const double rinv = 1.0 / r;

      // unit normal vector
      const double enx = cdata.delta[0] * rinv;
      const double eny = cdata.delta[1] * rinv;
      const double enz = cdata.delta[2] * rinv;
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const int itype = atom->type[cdata.i];
      const int jtype = atom->type[cdata.j];

      const double reff = (ri*rj)/(ri+rj);
      const double hIJ = (r-ri-rj);

      const double Fn_coh = calcCohesiveForce(hIJ, reff, itype, jtype);

      // apply normal force
      if(cdata.is_wall) {
        const double Fn_ = Fn_coh * cdata.area_ratio;
        i_forces.delta_F[0] += Fn_ * enx;
        i_forces.delta_F[1] += Fn_ * eny;
        i_forces.delta_F[2] += Fn_ * enz;
      } else {
        const double fx = Fn_coh * enx;
        const double fy = Fn_coh * eny;
        const double fz = Fn_coh * enz;

        i_forces.delta_F[0] += fx;
        i_forces.delta_F[1] += fy;
        i_forces.delta_F[2] += fz;

        j_forces.delta_F[0] -= fx;
        j_forces.delta_F[1] -= fy;
        j_forces.delta_F[2] -= fz;
      }

      cdata.has_force_update = true;
    }

  private:
    inline double calcCohesiveForce(const double hIJ, const double reff, const int itype, const int jtype)
    {
        if (hIJ > hCutEff[itype][jtype]) {
            return -aHamakerEff[itype][jtype]/((6*hIJ*hIJ)*reff);
        }
        return -aHamakerEff[itype][jtype]/(6*hCutEff[itype][jtype]*hCutEff[itype][jtype]*reff);
    }

    double **aHamakerEff;
    double **hCutEff;
    double **hMaxEff;
  };
}
#endif // COHESION_MODEL_HAMAKER_H_
#endif