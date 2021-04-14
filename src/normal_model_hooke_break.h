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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */
#ifdef NORMAL_MODEL
NORMAL_MODEL(HOOKE_BREAK,hooke/break,8)
#else
#ifndef NORMAL_MODEL_HOOKE_BREAK_H_
#define NORMAL_MODEL_HOOKE_BREAK_H_
#include "contact_models.h"
#include <math.h>
#include "atom.h"
#include "force.h"
#include "update.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class NormalModel<HOOKE_BREAK> : protected Pointers
  {
  public:
    static const int MASK = CM_REGISTER_SETTINGS | CM_CONNECT_TO_PROPERTIES | CM_COLLISION;

    NormalModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      Yeff(NULL),
      Geff(NULL),
      coeffRestMax(NULL),
      coeffRestLog(NULL),
      coeffMu(NULL),
      coeffStc(NULL),
      charVel(0.0),
      kappa(2./7.), // Landry et al., Phys. Rev. E 67 (041303), 1-9 (2003)
      viscous(false),
      tangential_damping(false),
      limitForce(false),
      ktToKn(false),
      displayedSettings(false)
    {
      history_offset = hsetup->add_history_value("deltaMax", "0");
      hsetup->add_history_value("sibling", "0");
      hsetup->add_history_value("siblingDeltaMax", "0");
      hsetup->add_history_value("collisionFactor", "0");
      hsetup->add_history_value("impactEnergy", "0");
      hsetup->add_history_value("forceMax", "0");
      hsetup->add_history_value("enx", "1");
      hsetup->add_history_value("eny", "1");
      hsetup->add_history_value("enz", "1");
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "HOOKE/BREAK loaded\n");
    }

    inline void registerSettings(Settings & settings)
    {
      settings.registerOnOff("viscous", viscous);
      settings.registerOnOff("tangential_damping", tangential_damping, true);
      settings.registerOnOff("limitForce", limitForce);
      settings.registerOnOff("ktToKnUser", ktToKn);
    }

    inline void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("Yeff", &MODEL_PARAMS::createYeff);
      registry.registerProperty("Geff", &MODEL_PARAMS::createGeff);
      registry.registerProperty("charVel", &MODEL_PARAMS::createCharacteristicVelocity);

      registry.connect("Yeff", Yeff,"model hooke/break");
      registry.connect("Geff", Geff,"model hooke/break");
      registry.connect("charVel", charVel,"model hooke/break");

      if(viscous) {
        registry.registerProperty("coeffMu", &MODEL_PARAMS::createCoeffMu);
        registry.registerProperty("coeffStc", &MODEL_PARAMS::createCoeffStc);
        registry.registerProperty("coeffRestMax", &MODEL_PARAMS::createCoeffRestMax);

        registry.connect("coeffMu", coeffMu,"model hooke/break viscous");
        registry.connect("coeffStc", coeffStc,"model hooke/break viscous");
        registry.connect("coeffRestMax", coeffRestMax,"model hooke/break viscous");
        //registry.connect("log(coeffRestMax)+coeffStc", logRestMaxPlusStc);
      } else {
        registry.registerProperty("coeffRestLog", &MODEL_PARAMS::createCoeffRestLog);

        registry.connect("coeffRestLog", coeffRestLog,"model hooke/break viscous");
      }

      if(ktToKn) {
        registry.registerProperty("stiffnessRatio", &MODEL_PARAMS::createStiffnessRatio);
        registry.connect("stiffnessRatio", kappa,"model hooke/break");
      }
      //NP modified C.K.
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"model hooke/break");
    }

    // effective exponent for stress-strain relationship
    //NP used for area correction of heat transfer
    inline double stressStrainExponent()
    {
      return 1.;
    }

    inline void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      const int itype = cdata.itype;
      const int jtype = cdata.jtype;
      double ri = cdata.radi;
      double rj = cdata.radj;
      double reff=cdata.is_wall ? ri : (ri*rj/(ri+rj));
      double meff=cdata.meff;
      double coeffRestLogChosen;

      double deltan = cdata.deltan;
      double * const history = &cdata.contact_history[history_offset];
      double * const deltaMax = &cdata.contact_history[history_offset];
      double sibling = history[1];
      double collisionFactor = history[3];
      double * const forceMax = &cdata.contact_history[history_offset+5];
      double * const en = &cdata.contact_history[history_offset+6];

      // limit forces between siblings
      //NP sibling contact flags must be set after insertion (pre_exchange)
      //NP and before this force calculation (pre_force)
      if (sibling) {
        double * const siblingDeltaMax = &cdata.contact_history[history_offset+2];
        if (deltan > *siblingDeltaMax) {
          if(*siblingDeltaMax == 0.0) {
            *siblingDeltaMax = deltan; // part of the overlap that needs to be scaled
            deltan *= collisionFactor;
          } else {
            deltan -= *siblingDeltaMax;
            deltan += *siblingDeltaMax*collisionFactor;
          }
        } else {
          *siblingDeltaMax = deltan;
          deltan *= collisionFactor;
        }

        cdata.deltan = deltan;
        cdata.shearupdate = false;
        cdata.vtr1 = cdata.vtr2 = cdata.vtr3 = 0.0;
      } else { // sibling contact is never an impact
        // detect new contact / maximum overlap
        if (*deltaMax == 0.0) { // new contact
          *deltaMax = deltan;
          history[4] = 0.5 * cdata.vn * cdata.vn; // mass specific impact energy
        } else if (*deltaMax > deltan || *deltaMax < 0.0) {
          *deltaMax = -deltan;
          history[4] = 0.0;
        } else { //if (*deltaMax <= deltan)
          *deltaMax = deltan;
          history[4] = 0.0;
        }
      }

      const double sqrtval = sqrt(reff);

      if (!displayedSettings) {
        displayedSettings = true;
        /*
        if (comm->me == 0 && screen) {
          if(ktToKn)
            fprintf(screen," NormalModel<HOOKE_BREAK>: will use user-modified ktToKn of %f.\n",kappa);
          if(tangential_damping)
            fprintf(screen," NormalModel<HOOKE_BREAK>: will apply tangential damping.\n");
          if(viscous)
            fprintf(screen," NormalModel<HOOKE_BREAK>: will apply damping based on Stokes number.\n");
          if(limitForce)
            fprintf(screen," NormalModel<HOOKE_BREAK>: will limit normal force.\n");
        }
        */
      }
      if (viscous) {
        // Stokes Number from MW Schmeeckle (2001)
        const double stokes=cdata.meff*cdata.vn/(6.0*M_PI*coeffMu[itype][jtype]*reff*reff);
        // Empirical from Legendre (2006)
        coeffRestLogChosen=log(coeffRestMax[itype][jtype])+coeffStc[itype][jtype]/stokes;
      } else {
        coeffRestLogChosen=coeffRestLog[itype][jtype];
      }

      const double k = (16./15.)*sqrtval*(Yeff[itype][jtype]);
      double kn = k*pow(meff*charVel*charVel/k,0.2);
      double kt = kn;
      if(ktToKn) kt *= kappa;
      const double coeffRestLogChosenSq = coeffRestLogChosen*coeffRestLogChosen;
      const double gamman = sqrt(4.*meff*kn*coeffRestLogChosenSq/(coeffRestLogChosenSq+M_PI*M_PI));
      const double gammat = tangential_damping ? gamman : 0.0;

      // convert Kn and Kt from pressure units to force/distance^2
      kn /= force->nktv2p;
      kt /= force->nktv2p;

      const double Fn_damping = -gamman*cdata.vn;
      const double Fn_contact = kn*deltan;
      double Fn = Fn_damping + Fn_contact;

      //limit force to avoid the artefact of negative repulsion force
      if (limitForce && (Fn<0.0) ) {
        Fn = 0.0;
      }

      cdata.Fn = Fn;
      cdata.kn = kn;
      cdata.kt = kt;
      cdata.gamman = gamman;
      cdata.gammat = gammat;

      // store for max force breakage criterion
      *forceMax = Fn;
      en[0] = cdata.en[0];
      en[1] = cdata.en[1];
      en[2] = cdata.en[2];

      // apply normal force
      if (cdata.is_wall) {
        const double Fn_ = Fn * cdata.area_ratio;
        i_forces.delta_F[0] = Fn_ * cdata.en[0];
        i_forces.delta_F[1] = Fn_ * cdata.en[1];
        i_forces.delta_F[2] = Fn_ * cdata.en[2];
      } else {
        const double fx = cdata.Fn * cdata.en[0];
        const double fy = cdata.Fn * cdata.en[1];
        const double fz = cdata.Fn * cdata.en[2];

        i_forces.delta_F[0] = fx;
        i_forces.delta_F[1] = fy;
        i_forces.delta_F[2] = fz;

        j_forces.delta_F[0] = -fx;
        j_forces.delta_F[1] = -fy;
        j_forces.delta_F[2] = -fz;
      }
    }

    inline void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
    double ** Yeff;
    double ** Geff;
    double ** coeffRestMax;
    double ** coeffRestLog;
    double ** coeffMu;
    double ** coeffStc;
    double charVel;
    double kappa;

    bool viscous;
    bool tangential_damping;
    bool limitForce;
    bool ktToKn;
    bool displayedSettings;
    int history_offset;
  };
}
}
#endif // NORMAL_MODEL_HOOKE_BREAK_H_
#endif
