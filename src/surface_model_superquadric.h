/* ----------------------------------------------------------------------
    This file is part of the LIGGGHTS DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Alexander Podlozhnyuk (DCS Computing GmbH, Linz)

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_SUPERQUADRIC,superquadric,5)
#else
#ifndef SURFACE_MODEL_SUPERQUADRIC_H_
#define SURFACE_MODEL_SUPERQUADRIC_H_

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#include "contact_models.h"
#include <cmath>
#include <algorithm>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "math_extra_liggghts_nonspherical.h"
#include "math_extra_liggghts_superquadric.h"

namespace LIGGGHTS {
namespace ContactModels
{
  template<>
  class SurfaceModel<SURFACE_SUPERQUADRIC> : protected Pointers
  {
    int inequality_start_offset;
    int particles_were_in_contact_offset;
    int contact_point_offset;
    int alpha1_offset;
    int alpha2_offset;
    Superquadric particle_i;
    Superquadric particle_j;
    enum {SURFACES_FAR, SURFACES_CLOSE, SURFACES_INTERSECT};

  public:
    static const int MASK = CM_COLLISION;

    SurfaceModel(LAMMPS * lmp, IContactHistorySetup* hsetup) :
        Pointers(lmp)
    {
      if(!atom->superquadric_flag)
        error->one(FLERR,"Applying surface model superquadric to a non-superquadric particle!");

      inequality_start_offset = hsetup->add_history_value("inequality_obb","0");
      particles_were_in_contact_offset = hsetup->add_history_value("particles_were_in_contact","0");
      contact_point_offset = hsetup->add_history_value("cpx", "0");
      hsetup->add_history_value("cpy", "0");
      hsetup->add_history_value("cpz", "0");
      alpha1_offset = hsetup->add_history_value("a1", "0");
      alpha2_offset = hsetup->add_history_value("a2", "0");
    }

    inline void registerSettings(Settings& settings)
    {
      settings.registerDoubleSetting("curvatureLimitFactor",curvatureLimitFactor, 0.0);
      settings.registerYesNo("meanCurvature", meanCurvature, false);
      settings.registerYesNo("gaussianCurvature", gaussianCurvature, false);
      if(curvatureLimitFactor < 0.0)
        error->one(FLERR,"Curvature limiter cannot be negative!");
      if(!meanCurvature && !gaussianCurvature)
        curvatureLimitFactor = 0.0;
      if(meanCurvature && gaussianCurvature)
        error->one(FLERR,"meanCurvature and gaussianCurvature cannot be simultaneously yes !");
    }

    inline void connectToProperties(PropertyRegistry&) {}

    inline bool checkSurfaceIntersect(CollisionData & cdata)
    {
      cdata.is_non_spherical = true;
      bool particles_in_contact = false;
      double *const prev_step_point = &cdata.contact_history[contact_point_offset]; //contact points
      double *const inequality_start = &cdata.contact_history[inequality_start_offset];
      double *const particles_were_in_contact = &cdata.contact_history[particles_were_in_contact_offset];

      const int iPart = cdata.i;
      const int jPart = cdata.j;

      particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
      particle_j.set(atom->x[jPart], atom->quaternion[jPart], atom->shape[jPart], atom->blockiness[jPart]);

      unsigned int int_inequality_start = MathExtraLiggghtsNonspherical::round_int(*inequality_start);

      bool obb_intersect = false;
      if(*particles_were_in_contact == SURFACES_INTERSECT)
        obb_intersect = true; //particles had overlap on the previous time step, skipping OBB intersection check
      else
        obb_intersect = MathExtraLiggghtsNonspherical::obb_intersect(&particle_i, &particle_j, int_inequality_start);

      if(obb_intersect) {//OBB intersect particles in possible contact

        double fi, fj;
        const double ri = cbrt(particle_i.shape[0]*particle_i.shape[1]*particle_i.shape[2]);
        const double rj = cbrt(particle_j.shape[0]*particle_j.shape[1]*particle_j.shape[2]);
        double ratio = ri / (ri + rj);

        if(*particles_were_in_contact == SURFACES_FAR)
          MathExtraLiggghtsNonspherical::calc_contact_point_if_no_previous_point_avaialable(cdata, &particle_i, &particle_j, cdata.contact_point, fi, fj, this->error);
        else
          MathExtraLiggghtsNonspherical::calc_contact_point_using_prev_step(cdata, &particle_i, &particle_j, ratio, update->dt, prev_step_point, cdata.contact_point, fi, fj, this->error);
        vectorCopy3D(cdata.contact_point, prev_step_point); //store contact point in contact history for the next DEM time step

        particles_in_contact = std::max(fi, fj) < 0.0;

        if(particles_in_contact) {
          double contact_point_i[3], contact_point_j[3];
          vectorSubtract3D(particle_j.gradient, particle_i.gradient, cdata.en);
          vectorNormalize3D(cdata.en); //normalize

          double *const alpha_i = &cdata.contact_history[alpha1_offset];
          double *const alpha_j = &cdata.contact_history[alpha2_offset];

          cdata.deltan = MathExtraLiggghtsNonspherical::extended_overlap_algorithm(&particle_i, &particle_j, cdata.en, alpha_i, alpha_j,
                cdata.contact_point, contact_point_i, contact_point_j, cdata.delta);

          cdata.reff = cdata.radi*cdata.radj / (cdata.radi + cdata.radj);
          if(curvatureLimitFactor > 0.0) {
              int curvature_radius_method = meanCurvature ? 0 : 1;
              double koefi = particle_i.calc_curvature_coefficient(curvature_radius_method, contact_point_i); //mean curvature coefficient
              double koefj = particle_j.calc_curvature_coefficient(curvature_radius_method, contact_point_j); //mean curvature coefficient
              cdata.reff = MathExtraLiggghtsNonspherical::get_effective_radius(cdata, atom->blockiness[iPart], atom->blockiness[jPart], koefi, koefj, curvatureLimitFactor, this->error);
          }

          *particles_were_in_contact = SURFACES_INTERSECT;
        } else
          *particles_were_in_contact = SURFACES_CLOSE;
       } else
         *particles_were_in_contact = SURFACES_FAR;
      *inequality_start = static_cast<double>(int_inequality_start);
      return particles_in_contact;
    }

    inline void collision(CollisionData & cdata, ForceData&, ForceData&)
    {
      if(cdata.is_wall) {
        cdata.reff = cdata.radi;
        if(curvatureLimitFactor > 0.0) {
          const int iPart = cdata.i;
          particle_i.set(atom->x[iPart], atom->quaternion[iPart], atom->shape[iPart], atom->blockiness[iPart]);
          int curvature_radius_method = meanCurvature ? 0 : 1;
          double koefi = particle_i.calc_curvature_coefficient(curvature_radius_method, cdata.contact_point); //mean curvature coefficient
          cdata.reff = MathExtraLiggghtsNonspherical::get_effective_radius_wall(cdata, atom->blockiness[iPart], koefi, curvatureLimitFactor, this->error);
        }
      }
      MathExtraLiggghtsNonspherical::surfacesIntersectNonSpherical(cdata, atom->x);
    }

    inline void noCollision(ContactData&, ForceData&, ForceData&){}
    void beginPass(CollisionData&, ForceData&, ForceData&){}
    void endPass(CollisionData&, ForceData&, ForceData&){}

  protected:
     double curvatureLimitFactor;
     bool meanCurvature;
     bool gaussianCurvature;
  };
}
}

#else  // SUPERQUADRIC_ACTIVE_FLAG

// add dummy class
namespace LIGGGHTS
{
namespace ContactModels
{
template<>
class SurfaceModel<SURFACE_SUPERQUADRIC> : protected Pointers
{
public:
    static const int MASK = CM_COLLISION;

    SurfaceModel(LAMMPS * lmp, IContactHistorySetup * hsetup) :
      Pointers(lmp)
    {
        error->all(FLERR, "Surface model used without support for superquadric particles");
    }
    inline void beginPass(CollisionData& , ForceData& , ForceData& ) {}
    inline void endPass(CollisionData& , ForceData& , ForceData& ) {}
    inline void registerSettings(Settings&) {}
    inline void connectToProperties(PropertyRegistry&) {}
    inline bool checkSurfaceIntersect(CollisionData&) { return false; }
    inline void collision(CollisionData&, ForceData&, ForceData&) {}
    inline void noCollision(ContactData&, ForceData&, ForceData&) {}
  };
}
}
#endif // SUPERQUADRIC_ACTIVE_FLAG

#endif // SURFACE_MODEL_SUPERQUADRIC_H_
#endif // SURFACE_MODEL
