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


#ifdef FIX_CLASS

FixStyle(particletemplate/sphere,FixTemplateSphere)

#else

#ifndef LMP_FIX_TEMPLATE_SPHERE_H
#define LMP_FIX_TEMPLATE_SPHERE_H

#include <vector>
#include "fix.h"
#include "fix_property_atom.h"
#include "probability_distribution.h"
#include "math_const.h"

namespace PARTICLE_PACKING
{

class Sphere
{
public:
    Sphere() :
        pos_x(0.0),
        pos_y(0.0),
        pos_z(0.0),
        radius(0.0),
        density(0.0),
        is_local(false),
        id(0)
    {}

    Sphere(const double * const _x, const double _radius, const double _density, const int _id) :
        pos_x(_x[0]),
        pos_y(_x[1]),
        pos_z(_x[2]),
        radius(_radius),
        density(_density),
        is_local(false),
        id(_id)
    {}

    void move_particle(const double * const shift)
    {
        pos_x += shift[0];
        pos_y += shift[1];
        pos_z += shift[2];
    }

    inline void get_pos(double * const _x) const
    {
        _x[0] = pos_x;
        _x[1] = pos_y;
        _x[2] = pos_z;
    }

    void set_local()
    { is_local = true; }

    void unset_local()
    { is_local = false; }

    bool get_local() const
    { return is_local; }

    double get_radius() const
    { return radius; }

    double get_density() const
    { return density; }

    double get_volume() const
    { return radius*radius*radius*LAMMPS_NS::MathConst::MY_4PI3; }

    int get_id() const
    { return id; }

    size_t n_fix_properties() const
    { return fix_properties.size(); }

    size_t fix_property_nentries(const int i) const
    { return fix_property_values[i].size(); }

    LAMMPS_NS::FixPropertyAtom* get_fix_property(const int i) const
    { return fix_properties[i]; }

    double fix_property_value(const int i, const int j) const
    { return fix_property_values[i][j]; }

    void init_fix_properties(std::vector<std::pair<LAMMPS_NS::FixPropertyAtom*, int> > &fpa_list)
    {
        std::vector<std::pair<LAMMPS_NS::FixPropertyAtom*, int> >::iterator it = fpa_list.begin();
        size_t n = 0;
        for (; it != fpa_list.end(); it++)
        {
            if (it->first)
                n++;
        }
        fix_properties.resize(n);
        fix_property_values.resize(n);
        for (size_t i = 0, j = 0; i < fpa_list.size(); i++)
        {
            if (!fpa_list[i].first)
                continue;
            fix_properties[j] = fpa_list[i].first;
            size_t m = fpa_list[i].second;
            fix_property_values[j].resize(m);
            j++;
        }
    }

    void set_fix_property_values(const int i, const double * const data)
    {
        for (size_t j = 0; j < fix_property_values[i].size(); j++)
            fix_property_values[i][j] = data[j];
    }

protected:
    double pos_x, pos_y, pos_z, radius, density;
    bool is_local;
    int id;
    std::vector<LAMMPS_NS::FixPropertyAtom*> fix_properties;
    std::vector<std::vector<double> > fix_property_values;
};

}

namespace LAMMPS_NS {

class FixTemplateSphere : public Fix {
 public:

  FixTemplateSphere(class LAMMPS *, int, char **);
  ~FixTemplateSphere();

  // inherited from Fix
  virtual void post_create(){}
  virtual int setmask();
  void write_restart(FILE *);
  void restart(char *);

  // access to protected properties
  virtual double volexpect();           //NP expectancy value for particle volume
  virtual double massexpect();          //NP expectancy value for particle mass
  virtual double min_rad() const;
  virtual double max_rad() const;
  virtual double max_r_bound() const;
  virtual int number_spheres();
  virtual int maxtype();
  virtual int mintype();
  class Region *region();

  // many particle generation, used by fix insert commands
  virtual void init_ptilist(int);
  virtual void delete_ptilist();
  virtual void randomize_ptilist(int,int);
  int n_pti_max;
  class ParticleToInsert **pti_list;

  virtual void finalize_insertion() {}

  ParticleToInsert* get_single_random_pti(int distribution_groupbit);

protected:

  int iarg;

  class Region *reg;
  class FixRegionVariable *reg_var;

  // random generator
  class RanPark *random;
  int seed;

  // properties of particle template
  int atom_type;
  class LMP_PROBABILITY_NS::PDF *pdf_radius;   //NP radius of each sphere
  class LMP_PROBABILITY_NS::PDF *pdf_density;
  double volume_expect;
  double mass_expect;
  double vol_limit;

  void randomize_single_pti(ParticleToInsert* &pti, int distribution_groupbit);
};

}

#endif
#endif
