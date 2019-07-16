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
#ifdef ROLLING_MODEL
ROLLING_MODEL(ROLLING_EPSD,epsd,2)
#else
#ifndef ROLLING_MODEL_EPSD_H_
#define ROLLING_MODEL_EPSD_H_
#include "contact_models.h"
#include <algorithm>
#include <math.h>
#include "domain.h"
#include "math_extra_liggghts.h"
#ifdef SUPERQUADRIC_ACTIVE_FLAG
#include "math_extra_liggghts_nonspherical.h"
#endif

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace LAMMPS_NS;

  template<>
  class RollingModel<ROLLING_EPSD> : protected Pointers
  {
  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    RollingModel(class LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp), coeffRollFrict(NULL), coeffRollVisc(NULL)
    {
      history_offset = hsetup->add_history_value("r_torquex_old", "1");
      hsetup->add_history_value("r_torquey_old", "1");
      hsetup->add_history_value("r_torquez_old", "1");
      /*NL*/ if(comm->me == 0 && screen) fprintf(screen, "EPSD loaded\n");
    }

    void registerSettings(Settings&) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("coeffRollFrict", &MODEL_PARAMS::createCoeffRollFrict);
      registry.registerProperty("coeffRollVisc", &MODEL_PARAMS::createCoeffRollVisc);
      registry.connect("coeffRollFrict", coeffRollFrict,"rolling_model epsd");
      registry.connect("coeffRollVisc", coeffRollVisc,"rolling_model epsd");

      //NP modified C.K.
      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"rolling model epsd");
    }

    void collision(CollisionData & cdata, ForceData & i_forces, ForceData & j_forces) //NP modified C.K.
    {
      double r_torque[3];
      vectorZeroize3D(r_torque);

      if(cdata.touch) *cdata.touch |= TOUCH_ROLLING_MODEL;

      const double radi = cdata.radi;
      const double radj = cdata.radj;
      double reff = cdata.is_wall ? radi : (radi*radj/(radi+radj));

#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(cdata.is_non_spherical && atom->superquadric_flag) {
        reff = cdata.reff;
      }
#endif
      if(cdata.is_wall) {
        const double wr1 = cdata.wr1;
        const double wr2 = cdata.wr2;
        const double wr3 = cdata.wr3;

        double r_inertia = 0.0;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(cdata.is_non_spherical) {
          const double rii = pointDistance(cdata.contact_point, atom->x[cdata.i]);
          const double omega_mag = sqrt(wr1*wr1 + wr2*wr2 + wr3*wr3); // TODO move sqrt after if
          if(omega_mag != 0.0) {
            const int i = cdata.i;
            double er[3];
            er[0] = wr1 / omega_mag;
            er[1] = wr2 / omega_mag;
            er[2] = wr3 / omega_mag;
            const double Ix = atom->inertia[i][0];
            const double Iy = atom->inertia[i][1];
            const double Iz = atom->inertia[i][2];
            double inertia_tensor[9];
            double inertia_tensor_local[9] = { Ix, 0.0, 0.0,
                                               0.0, Iy, 0.0,
                                               0.0, 0.0, Iz };
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local, atom->quaternion[cdata.i], inertia_tensor);
            double temp[3];
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor, er, temp);
            double Ii = MathExtra::dot3(temp, er);
            r_inertia = Ii + cdata.mi*rii*rii;
          }
        } else {
          if (domain->dimension == 2) r_inertia = 1.5*cdata.mi*reff*reff;
          else  r_inertia = 1.4*cdata.mi*reff*reff;
        }
#else
        if (domain->dimension == 2) r_inertia = 1.5*cdata.mi*reff*reff;
        else  r_inertia = 1.4*cdata.mi*reff*reff;
#endif

        calcRollTorque(r_torque,cdata,reff,wr1,wr2,wr3,r_inertia);


      } else {
        double  wr_roll[3];

        const int i = cdata.i;
        const int j = cdata.j;

        const double * const * const omega = atom->omega;
        // relative rotational velocity
        vectorSubtract3D(omega[i],omega[j],wr_roll);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        double r_inertia_red_i, r_inertia_red_j;
        double r_inertia = 0.0;
        if (cdata.is_non_spherical) {
          const double rii = pointDistance(cdata.contact_point, atom->x[i]);
          const double rjj = pointDistance(cdata.contact_point, atom->x[j]);
          const double omega_mag = vectorMag3D(wr_roll);
          if (omega_mag != 0.0) {
            double er[3];
            er[0] = wr_roll[0] / omega_mag;
            er[1] = wr_roll[1] / omega_mag;
            er[2] = wr_roll[2] / omega_mag;
            const double Ix_i = atom->inertia[i][0];
            const double Iy_i = atom->inertia[i][1];
            const double Iz_i = atom->inertia[i][2];

            const double Ix_j = atom->inertia[j][0];
            const double Iy_j = atom->inertia[j][1];
            const double Iz_j = atom->inertia[j][2];

            double inertia_tensor_i[9];
            double inertia_tensor_local_i[9] = { Ix_i, 0.0, 0.0,
                                                 0.0, Iy_i, 0.0,
                                                 0.0, 0.0, Iz_i };
            double inertia_tensor_j[9];
            double inertia_tensor_local_j[9] = { Ix_j, 0.0, 0.0,
                                                 0.0, Iy_j, 0.0,
                                                 0.0, 0.0, Iz_j };
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local_i, atom->quaternion[i], inertia_tensor_i);
            MathExtraLiggghtsNonspherical::tensor_quat_rotate(inertia_tensor_local_j, atom->quaternion[j], inertia_tensor_j);
            double temp[3];
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor_i, er, temp);
            double Ii = MathExtra::dot3(temp, er);
            MathExtraLiggghtsNonspherical::matvec(inertia_tensor_j, er, temp);
            double Ij = MathExtra::dot3(temp, er);
            r_inertia_red_i = Ii + cdata.mi*rii*rii;
            r_inertia_red_j = Ij + cdata.mj*rjj*rjj;
            r_inertia = r_inertia_red_i*r_inertia_red_j / (r_inertia_red_i + r_inertia_red_j);
          }
        } else {
          r_inertia_red_i = cdata.mi*radi*radi;
          r_inertia_red_j = cdata.mj*radj*radj;
          if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
          else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
        }
#else
        const double r_inertia_red_i = cdata.mi*radi*radi;
        const double r_inertia_red_j = cdata.mj*radj*radj;
        double r_inertia;
        if (domain->dimension == 2) r_inertia = 1.5 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
        else  r_inertia = 1.4 * r_inertia_red_i * r_inertia_red_j/(r_inertia_red_i + r_inertia_red_j);
#endif

        calcRollTorque(r_torque,cdata,reff,wr_roll[0],wr_roll[1],wr_roll[2],r_inertia);
      }

      i_forces.delta_torque[0] -= r_torque[0];
      i_forces.delta_torque[1] -= r_torque[1];
      i_forces.delta_torque[2] -= r_torque[2];
      j_forces.delta_torque[0] += r_torque[0];
      j_forces.delta_torque[1] += r_torque[1];
      j_forces.delta_torque[2] += r_torque[2];
    }

    void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      if(cdata.touch) *cdata.touch &= ~TOUCH_ROLLING_MODEL;
      double * const c_history = &cdata.contact_history[history_offset];
      c_history[0] = 0.0; // this is the r_torque_old
      c_history[1] = 0.0; // this is the r_torque_old
      c_history[2] = 0.0; // this is the r_torque_old
    }

    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  private:
    double ** coeffRollFrict;
    double ** coeffRollVisc;
    int history_offset;

    inline void calcRollTorque(double (&r_torque)[3],const CollisionData & cdata,double reff,double wr1,double wr2,double wr3,double r_inertia) {

      double wr_n[3],wr_t[3],dr_torque[3];

      const int itype = cdata.itype;
      const int jtype = cdata.jtype;

      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      const double dt = update->dt; //NP TODO: any transformation of the timestep needed for other unit systems?

      double * const c_history = &cdata.contact_history[history_offset]; // requires Style::TANGENTIAL == TANGENTIAL_HISTORY
      const double rmu= coeffRollFrict[itype][jtype];

      // remove normal (torsion) part of relative rotation
      // use only tangential parts for rolling torque
      const double wr_dot_delta = wr1*enx+ wr2*eny + wr3*enz;
      wr_n[0] = enx * wr_dot_delta;
      wr_n[1] = eny * wr_dot_delta;
      wr_n[2] = enz * wr_dot_delta;
      wr_t[0] = wr1 - wr_n[0];
      wr_t[1] = wr2 - wr_n[1];
      wr_t[2] = wr3 - wr_n[2];

      // spring
      const double kr = 2.25*cdata.kn*rmu*rmu*reff*reff; //NP modified A.A.; Not sure if kr is right for 3D;

      vectorScalarMult3D(wr_t,dt*kr,dr_torque);

      r_torque[0] = c_history[0] + dr_torque[0];
      r_torque[1] = c_history[1] + dr_torque[1];
      r_torque[2] = c_history[2] + dr_torque[2];

      // limit max. torque
      const double r_torque_mag = vectorMag3D(r_torque);
      const double r_torque_max = fabs(cdata.Fn)*reff*rmu;
      if(r_torque_mag > r_torque_max)
      {
        //printf("[%d] %e > %e\n", update->ntimestep, r_torque_mag, r_torque_max);
        const double factor = r_torque_max / r_torque_mag;

        r_torque[0] *= factor;
        r_torque[1] *= factor;
        r_torque[2] *= factor;

        // save rolling torque due to spring
        if (cdata.shearupdate && cdata.computeflag) {
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];
        }

        // no damping / no dashpot in case of full mobilisation rolling angle

      } else {
        // save rolling torque due to spring before adding damping torque
        if (cdata.shearupdate && cdata.computeflag) {
          c_history[0] = r_torque[0];
          c_history[1] = r_torque[1];
          c_history[2] = r_torque[2];
        }

        // dashpot
        /*NL*/ //if (screen) fprintf(screen,"Calc r_coef for types %i %i with coef= %e, r_inertia=%e and kr=%e\n",itype,jtype,coeffRollVisc[itype][jtype],r_inertia,kr);
        const double r_coef = coeffRollVisc[itype][jtype] * 2 * sqrt(r_inertia*kr);

        // add damping torque
        r_torque[0] += r_coef*wr_t[0];
        r_torque[1] += r_coef*wr_t[1];
        r_torque[2] += r_coef*wr_t[2];
      }
    }
  };
}
}
#endif // ROLLING_MODEL_EPSD_H_
#endif
