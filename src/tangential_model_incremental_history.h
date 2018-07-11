/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2018- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifdef TANGENTIAL_MODEL
TANGENTIAL_MODEL(TANGENTIAL_INCREMENTAL_HISTORY,incremental_history,2)
#else
#ifndef TANGENTIAL_MODEL_INCREMENTAL_HISTORY_H_
#define TANGENTIAL_MODEL_INCREMENTAL_HISTORY_H_
#include "contact_models.h"
#include <math.h>
#include "update.h"
#include "global_properties.h"
#include "atom.h"

namespace LIGGGHTS {
namespace ContactModels
{
#define TRACK_TANGENTIAL_OVERLAP
#define SCALE_OLD_TANGENTIAL_FORCE
  template<>
  class TangentialModel<TANGENTIAL_INCREMENTAL_HISTORY> : protected Pointers
  {
    double ** coeffFrict;
    int history_offset;

  public:
    static const int MASK = CM_CONNECT_TO_PROPERTIES | CM_COLLISION | CM_NO_COLLISION;

    TangentialModel(LAMMPS * lmp, IContactHistorySetup * hsetup) : Pointers(lmp),
      coeffFrict(NULL)
    {
#ifdef TRACK_TANGENTIAL_OVERLAP
      history_offset = hsetup->add_history_value("shearx", "1");
      hsetup->add_history_value("sheary", "1");
      hsetup->add_history_value("shearz", "1");
      hsetup->add_history_value("shearfx", "1");
      hsetup->add_history_value("shearfy", "1");
      hsetup->add_history_value("shearfz", "1");
#else
      history_offset = hsetup->add_history_value("shearfx", "1");
      hsetup->add_history_value("shearfy", "1");
      hsetup->add_history_value("shearfz", "1");
#endif
#ifdef SCALE_OLD_TANGENTIAL_FORCE
      hsetup->add_history_value("kt", "0");
      hsetup->add_history_value("Fn", "1");
#endif

      if (comm->me == 0 && screen) fprintf(screen, "TANGENTIAL/INCREMENTAL_HISTORY loaded\n");
    }

    inline void registerSettings(Settings&){}

    inline void connectToProperties(PropertyRegistry & registry)
    {
      registry.registerProperty("coeffFrict", &MODEL_PARAMS::createCoeffFrict);
      registry.connect("coeffFrict", coeffFrict,"tangential_model history");
    }

    inline void collision(const CollisionData & cdata, ForceData & i_forces, ForceData & j_forces)
    {
      // normal forces = Hookian contact + normal velocity damping
      const double enx = cdata.en[0];
      const double eny = cdata.en[1];
      const double enz = cdata.en[2];

      // shear history effects
      if (cdata.touch) *cdata.touch |= TOUCH_TANGENTIAL_MODEL;
#ifdef TRACK_TANGENTIAL_OVERLAP
      double * const shear = &cdata.contact_history[history_offset];
      double * const Ft_shear = &cdata.contact_history[history_offset+3];
#ifdef SCALE_OLD_TANGENTIAL_FORCE
      double* const kt_old = &cdata.contact_history[history_offset+6];
      double* const Fn_old = &cdata.contact_history[history_offset+7];
#endif
#else
      double * const Ft_shear = &cdata.contact_history[history_offset];
#ifdef SCALE_OLD_TANGENTIAL_FORCE
      double* const kt_old = &cdata.contact_history[history_offset+3];
      double* const Fn_old = &cdata.contact_history[history_offset+4];
#endif
#endif

      if (cdata.shearupdate && cdata.computeflag) {
        const double dt = update->dt;

        const double deltatx = cdata.vtr1 * dt;
        const double deltaty = cdata.vtr2 * dt;
        const double deltatz = cdata.vtr3 * dt;

#ifdef TRACK_TANGENTIAL_OVERLAP
        shear[0] += deltatx;
        shear[1] += deltaty;
        shear[2] += deltatz;

        // project vector into normal plane and
        // subtract normal component from shear displacement vector
        const double shear_dot_en = shear[0]*enx + shear[1]*eny + shear[2]*enz;
        shear[0] -= shear_dot_en * enx;
        shear[1] -= shear_dot_en * eny;
        shear[2] -= shear_dot_en * enz;
#endif
        const double Ft_shear_dot_en = Ft_shear[0]*enx + Ft_shear[1]*eny + Ft_shear[2]*enz;
        Ft_shear[0] -= Ft_shear_dot_en * enx;
        Ft_shear[1] -= Ft_shear_dot_en * eny;
        Ft_shear[2] -= Ft_shear_dot_en * enz;

        const double kt = cdata.kt;
        const double gammat = cdata.gammat;
        const double xmu = coeffFrict[cdata.itype][cdata.jtype];

#ifdef SCALE_OLD_TANGENTIAL_FORCE
        // rescale Ft_shear, see
        // Thornton et al., Powder Technol. 233 (2013) 30-46
        if((cdata.Fn - Fn_old[0]) < 0. && kt_old[0] > 0. && kt_old[0] != kt) {
          const double kt_ratio = kt/kt_old[0];
          Ft_shear[0] *= kt_ratio;
          Ft_shear[1] *= kt_ratio;
          Ft_shear[2] *= kt_ratio;
        }
        kt_old[0] = kt;
        Fn_old[0] = cdata.Fn;
#endif

        Ft_shear[0] -= kt * deltatx + gammat * cdata.vtr1;
        Ft_shear[1] -= kt * deltaty + gammat * cdata.vtr2;
        Ft_shear[2] -= kt * deltatz + gammat * cdata.vtr3;

        // rescale tangential forces if needed
        const double Ft_shear_mag = sqrt(Ft_shear[0]*Ft_shear[0] + Ft_shear[1]*Ft_shear[1] + Ft_shear[2]*Ft_shear[2]);
        const double Ft_friction = xmu * fabs(cdata.Fn);

        if (Ft_shear_mag > Ft_friction) {
          if (Ft_shear_mag > 0.0) {
            const double ratio = Ft_friction / Ft_shear_mag;
            Ft_shear[0] *= ratio;
            Ft_shear[1] *= ratio;
            Ft_shear[2] *= ratio;
          }
        }
      }

      const double Ft1 = Ft_shear[0];
      const double Ft2 = Ft_shear[1];
      const double Ft3 = Ft_shear[2];

      // forces & torques
      const double tor1 = eny * Ft3 - enz * Ft2;
      const double tor2 = enz * Ft1 - enx * Ft3;
      const double tor3 = enx * Ft2 - eny * Ft1;

      // return resulting forces
      if (cdata.is_wall) {
        const double area_ratio = cdata.area_ratio;
        i_forces.delta_F[0] += Ft1 * area_ratio;
        i_forces.delta_F[1] += Ft2 * area_ratio;
        i_forces.delta_F[2] += Ft3 * area_ratio;
        i_forces.delta_torque[0] = -cdata.cri * tor1 * area_ratio;
        i_forces.delta_torque[1] = -cdata.cri * tor2 * area_ratio;
        i_forces.delta_torque[2] = -cdata.cri * tor3 * area_ratio;
      } else {
        i_forces.delta_F[0] += Ft1;
        i_forces.delta_F[1] += Ft2;
        i_forces.delta_F[2] += Ft3;
        i_forces.delta_torque[0] = -cdata.cri * tor1;
        i_forces.delta_torque[1] = -cdata.cri * tor2;
        i_forces.delta_torque[2] = -cdata.cri * tor3;

        j_forces.delta_F[0] += -Ft1;
        j_forces.delta_F[1] += -Ft2;
        j_forces.delta_F[2] += -Ft3;
        j_forces.delta_torque[0] = -cdata.crj * tor1;
        j_forces.delta_torque[1] = -cdata.crj * tor2;
        j_forces.delta_torque[2] = -cdata.crj * tor3;
      }
    }

    inline void noCollision(ContactData & cdata, ForceData&, ForceData&)
    {
      // unset non-touching neighbors
      if (cdata.touch) *cdata.touch &= ~TOUCH_TANGENTIAL_MODEL;
      if (!cdata.contact_history) return;
      double * const shear = &cdata.contact_history[history_offset];
      int i = 0;
      shear[i++] = 0.0;
      shear[i++] = 0.0;
      shear[i++] = 0.0;
#ifdef TRACK_TANGENTIAL_OVERLAP
      shear[i++] = 0.0;
      shear[i++] = 0.0;
      shear[i++] = 0.0;
#endif
#ifdef SCALE_OLD_TANGENTIAL_FORCE
      shear[i++] = 0.0;
      shear[i++] = 0.0;
#endif
    }

    inline void beginPass(CollisionData&, ForceData&, ForceData&){}
    inline void endPass(CollisionData&, ForceData&, ForceData&){}
  };
}
}
#endif // TANGENTIAL_MODEL_INCREMENTAL_HISTORY_H_
#endif
