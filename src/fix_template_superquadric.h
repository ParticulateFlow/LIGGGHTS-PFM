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

    Copyright 2015-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef SUPERQUADRIC_ACTIVE_FLAG

#ifdef FIX_CLASS

FixStyle(particletemplate/superquadric,FixTemplateSuperquadric)

#else

#ifndef LMP_FIX_TEMPLATE_SUPERQUADRIC_H
#define LMP_FIX_TEMPLATE_SUPERQUADRIC_H

#include "fix.h"
#include "probability_distribution.h"
#include "fix_template_sphere.h"
#include "math_extra_liggghts_superquadric.h"

namespace PARTICLE_PACKING
{

class SQ: public Sphere
{
public:
    SQ():
        Sphere(),
        a(0.0),
        b(0.0),
        c(0.0),
        n1(0.0),
        n2(0.0)
    {}
    SQ(const double * const _x, const double * const _shape, const double * const _blockiness, const double _density, const int _id):
        Sphere(_x, 0.0, _density, _id),
        a(_shape[0]),
        b(_shape[1]),
        c(_shape[2]),
        n1(_blockiness[0]),
        n2(_blockiness[1])
    {
        MathExtraLiggghtsNonspherical::bounding_sphere_radius_superquadric(_shape, _blockiness, &radius);
    }

    double get_shape(int i) const
    {
        switch(i)
        {
        case 0:
            return a;
        case 1:
            return b;
        case 2:
            return c;
        }
        return -1;
    }

    double get_blockiness(int i) const
    {
        switch(i)
        {
        case 0:
            return n1;
        case 1:
            return n2;
        }
        return -1;
    }
  protected:
    double a, b, c, n1, n2;
};

}

namespace LAMMPS_NS {

class FixTemplateSuperquadric : public FixTemplateSphere {
 public:

  FixTemplateSuperquadric(class LAMMPS *, int, char **);
  ~FixTemplateSuperquadric();

  // access to protected properties
  virtual double min_rad() const;
  virtual double max_rad() const;
  virtual double max_r_bound() const;

  // many particle generation, used by fix insert commands
  virtual void init_ptilist(int n_random_max);
  virtual void randomize_ptilist(int,int);

 protected:

  bool relative;

  // properties of particle template
  class LMP_PROBABILITY_NS::PDF *pdf_shapex;
  class LMP_PROBABILITY_NS::PDF *pdf_shapey;
  class LMP_PROBABILITY_NS::PDF *pdf_shapez;
  class LMP_PROBABILITY_NS::PDF *pdf_size;
  class LMP_PROBABILITY_NS::PDF *pdf_blockiness1;
  class LMP_PROBABILITY_NS::PDF *pdf_blockiness2;
};

}

#endif // LMP_FIX_TEMPLATE_SUPERQUADRIC_H
#endif // FIX_CLASS

#endif // SUPERQUADRIC_ACTIVE_FLAG
