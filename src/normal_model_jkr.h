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
   Andreas Aigner (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(JKR,jkr,5)
#else
#ifndef NORMAL_MODEL_JKR_H_
#define NORMAL_MODEL_JKR_H_
#include "global_properties.h"
#include <math.h>
#include <limits>

namespace MODEL_PARAMS
{
  static MatrixProperty * createCoefContactRadius(PropertyRegistry & registry, const char * caller, bool)
  {
    const int max_type = registry.max_type();

    registry.registerProperty("Yeff", &createYeff);
    registry.registerProperty("cohEnergyDens", &createCohesionEnergyDensity);

    MatrixProperty * matrix = new MatrixProperty(max_type+1, max_type+1);
    MatrixProperty * YeffProp = registry.getMatrixProperty("Yeff", caller);
    MatrixProperty * cohEnergyDensProp = registry.getMatrixProperty("cohEnergyDens", caller);
    double * const * const Yeff = YeffProp->data;
    double * const * const cohEnergyDens = cohEnergyDensProp->data;

    for(int i=1;i< max_type+1; i++)
      {
        for(int j=1;j<max_type+1;j++)
          {
            matrix->data[i][j] = sqrt(M_PI*cohEnergyDens[i][j]/Yeff[i][j]);
          }
      }

    return matrix;
  }
}

namespace LIGGGHTS {

namespace ContactModels
{
  template<>
  class NormalModel<JKR> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup*) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      betaeff(NULL),
      tangential_damping(false),
      cohEnergyDens(NULL),
      coefContactRadius(NULL)
    {
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "JKR loaded\n");
    }

    void registerSettings(Settings & settings)
    {
      settings.registerOnOff("tangential_damping", tangential_damping, true);
    }

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff);
      registry.registerProperty("betaeff", &MODEL_PARAMS::createBetaEff);
      registry.registerProperty("cohEnergyDens", &MODEL_PARAMS::createCohesionEnergyDensity);
      registry.registerProperty("coefContactRadius", &MODEL_PARAMS::createCoefContactRadius);

      registry.connect("Yeff", Yeff, "normal_model jkr");
      registry.connect("Geff", Geff, "normal_model jkr");
      registry.connect("betaeff", betaeff, "normal_model jkr");
      registry.connect("cohEnergyDens", cohEnergyDens, "normal_model jkr");
      registry.connect("coefContactRadius", coefContactRadius, "normal_model jkr");

      //NP modified C.K.
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model jkr");
    }

    // effective exponent for stress-strain relationship
    //NP used for area correction of heat transfer
    inline double stressStrainExponent()
    {
      return 1.5;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff=cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      double meff=cdata.meff;

      double deltan=cdata.deltan;

      double sqrtval = sqrt(reff*cdata.deltan);

      double Sn=2.*Yeff[itype][jtype]*sqrtval;
      double St=8.*Geff[itype][jtype]*sqrtval;

      // original Hertz model (required for other models)
      double kn=4./3.*Yeff[itype][jtype]*sqrtval;
      double kt=St;
      const double sqrtFiveOverSix = 0.91287092917527685576161630466800355658790782499663875;
      double gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
      double gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
      if (!tangential_damping) gammat = 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      // calculate contact radius
      const double cRad = calcContactRadius(itype,jtype,reff,deltan);

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact_hertz = kn*cdata.deltan;
      const double cRad3 = cRad*cRad*cRad;
      const double Fn_contact = 4./3.*Yeff[itype][jtype]*cRad3/reff - sqrt(16.*M_PI*cohEnergyDens[itype][jtype]*Yeff[itype][jtype]*cRad3);
      const double Fn = Fn_contact + Fn_damping;
      const double Fn_hertz = Fn_contact_hertz + Fn_damping;

      //NP current implementation: As simplification we use Hertz force for the calculation of the tangential force
      //NP  Cohesive forces (negative forces) may corrupt Coloumb's law!
      cdata.Fn = Fn_hertz; // save original normal Hertz force for tangential forces and rolling friction
      //NP save stiffness and damping coefficient of the original Hertz model for rolling resistance model, ...
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;

      // apply normal force (use JKR Force)
      if(cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
      } else {
        i_forces.delta_F[0] = Fn * cdata.en[0];
        i_forces.delta_F[1] = Fn * cdata.en[1];
        i_forces.delta_F[2] = Fn * cdata.en[2];

        j_forces.delta_F[0] = -i_forces.delta_F[0];
        j_forces.delta_F[1] = -i_forces.delta_F[1];
        j_forces.delta_F[2] = -i_forces.delta_F[2];
      }
    }

    void noCollision(ContactData& cdata, ForceData& i_forces, ForceData& j_forces){
      //NP missing non contact force!!
      //NP There is no way for non-contact, particle-wall forces!!
      const int i = cdata.i;
      const int j = cdata.j;
      const int itype = atom->type[i];
      const int jtype = atom->type[j];
      const double ri = cdata.radi;
      const double rj = cdata.radj;
      const double reff=cdata.is_wall ? cdata.radi : (ri*rj/(ri+rj));
      const double radsum = cdata.radsum;

      // minimum deltan
      const double minDeltaN = -1.5*cbrt(reff*cohEnergyDens[itype][jtype]*cohEnergyDens[itype][jtype]*M_PI*M_PI/(2*Yeff[itype][jtype]*Yeff[itype][jtype]));
      const double r = sqrt(cdata.rsq);
      const double deltan = radsum - r;
      /*NL*/ //if (screen) fprintf(screen,"Non-contact: minDeltaN = %f, deltan = %f\n",minDeltaN,deltan);

      if (minDeltaN < deltan) { // non-contact force
        const double cRad = calcContactRadius(itype,jtype,reff,deltan);

        //NP current state: no damping!
        const double cRad3 = cRad*cRad*cRad;
        const double Fn = 4./3.*Yeff[itype][jtype]*cRad3/reff - sqrt(16.*M_PI*cohEnergyDens[itype][jtype]*Yeff[itype][jtype]*cRad3);

        double **x = atom->x;
        const double rinv = 1.0 / r;
        const double enx = (x[i][0] - x[j][0]) * rinv;
        const double eny = (x[i][1] - x[j][1]) * rinv;
        const double enz = (x[i][2] - x[j][2]) * rinv;


        // apply normal force
        if(cdata.is_wall) {
          const double Fn_ = Fn * cdata.area_ratio;
          i_forces.delta_F[0] = Fn_ * enx;
          i_forces.delta_F[1] = Fn_ * eny;
          i_forces.delta_F[2] = Fn_ * enz;
        } else {
          i_forces.delta_F[0] = Fn * enx;
          i_forces.delta_F[1] = Fn * eny;
          i_forces.delta_F[2] = Fn * enz;

          j_forces.delta_F[0] = -i_forces.delta_F[0];
          j_forces.delta_F[1] = -i_forces.delta_F[1];
          j_forces.delta_F[2] = -i_forces.delta_F[2];
        }

        /*NL*/ //if (screen) fprintf(screen,"Non-contact force: cRad = %f, Fn = %f\n",cRad,Fn);

      } // no interaction

    }
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    inline double calcContactRadius(const int itype, const int jtype, const double reff, const double deltan){
      // contact radius calculation (Newton faster than analytic solution)
      //NP analytic solution of a quartic equation
      double cRad = reff;
      for (int k=0; k<10; k++) {
        const double cRadOld = cRad;
        const double cRadOld2 = cRadOld*cRadOld;
        const double coef = coefContactRadius[itype][jtype]*reff*sqrt(cRadOld);
        cRad = (cRadOld*cRadOld2 + cRadOld*coef + deltan*cRadOld*reff)/(2*cRadOld2 - coef); //NP modified sign of deltan
        /*NL*/ //if (screen) fprintf(screen,"Loop: deltan = %f, coefContactRadius = %f, cRad = %f, cRadOld = %f\n",deltan,coefContactRadius[itype][jtype],cRad,cRadOld);
        if (fabs(cRad-cRadOld) < std::numeric_limits<double>::epsilon()) break;
      }
      /*NL*/ //if (isnan(cRad)) lmp->error->all(FLERR,"cRad is NAN. Cohesive surface energy to high?");
      /*NL*/ //if (screen) fprintf(screen,"Final: coefContactRadius = %f, cRad = %f\n",coefContactRadius[itype][jtype],cRad);

      return cRad;
    }

    double ** Yeff;
    double ** Geff;
    double ** betaeff;

    bool tangential_damping;

    //NP pre-calculated coefficients for the jkr model
    double ** cohEnergyDens;
    double ** coefContactRadius;
  };
}

}
#endif
#endif
