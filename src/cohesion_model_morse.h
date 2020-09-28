/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department of Particulate Flow Modelling
   Copyright 2020- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Achuth Balachandran Nair (JKU Linz)
------------------------------------------------------------------------- */
#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_MORSE,morse,5)
#else
#ifndef COHESION_MODEL_MORSE_H_
#define COHESION_MODEL_MORSE_H_
#include "contact_models.h"
#include "domain.h"


namespace MODEL_PARAMS
{
  static MatrixProperty* createMorseConstant(PropertyRegistry & registry, const char * caller, bool)
  {
    return createPerTypePairProperty(registry, "morseConstant", caller);
  }

  static MatrixProperty* createBetaConstant(PropertyRegistry & registry, const char * caller, bool)
  {
    return createPerTypePairProperty(registry, "betaConstant", caller);
  }

  static MatrixProperty* createCutOffDistance(PropertyRegistry & registry, const char * caller, bool)
  {
    return createPerTypePairProperty(registry, "cutoffDist", caller);
  }
}

namespace LIGGGHTS {
namespace ContactModels {
  using namespace LAMMPS_NS;

  template<>
  class CohesionModel<COHESION_MORSE> : protected Pointers {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    CohesionModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp), morse(NULL), beta(NULL), cutoff(NULL)
    {
      if(domain->dimension!=3)
          error->all(FLERR,"Cohesion model valid for 3d simulations only");
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "COHESION/MORSE loaded\n");
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("morseConstant", &MODEL_PARAMS::createMorseConstant);
      registry.registerProperty("betaConstant", &MODEL_PARAMS::createBetaConstant);
      registry.registerProperty("cutoffDist", &MODEL_PARAMS::createCutOffDistance);

      registry.connect("morseConstant", morse,"cohesion_model morse");
      registry.connect("betaConstant", beta,"cohesion_model morse");
      registry.connect("cutoffDist", cutoff,"cohesion_model morse");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model morse");
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) //NP modified C.K.
    {
      // since the particles are in contact:
      //     * the maximum cohesive force acts between the particles
      //     * the calculation of the displacement (hIJ) is not required (no cdata.r)
      const double r = sqrt(cdata.rsq);
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double dist = cdata.is_wall ? (r-ri) : (r-ri-rj);
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      // since the particles are in contact, the maximum cohesive force acts
      // using the equation directly without function call may be faster, since you can skip one if-statement
      const double Fn_coh = calcCohesiveForce(cdata, dist, r, itype, jtype);

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

      const double dist = cdata.is_wall ? (r-ri) : (r-ri-rj);

      const double Fn_coh = calcCohesiveForce(cdata, dist, r, itype, jtype);

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
    inline double calcCohesiveForce(ContactData& cdata,const double dist, const double r, const int itype, const int jtype)
    {
        double dr = dist - cutoff[itype][jtype];
        double dexp = exp(-beta[itype][jtype]*dr);

        double f_coh = morse[itype][jtype] * (dexp*dexp - dexp)/r;

        return f_coh;
    }

    double **morse;
    double **beta;
    double **cutoff;
  };
}
}
#endif // COHESION_MODEL_MORSE_H_
#endif
