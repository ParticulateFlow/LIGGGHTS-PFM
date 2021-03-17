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

FixStyle(particletemplate/multiplespheres,FixTemplateMultiplespheres)

#else

#ifndef LMP_FIX_TEMPLATE_MULTIPLESPHERES_H
#define LMP_FIX_TEMPLATE_MULTIPLESPHERES_H

#include "fix.h"
#include "fix_template_sphere.h"
#include "fix_property_atom.h"

namespace LAMMPS_NS {

class FixTemplateMultiplespheres : public FixTemplateSphere {
 public:

  FixTemplateMultiplespheres(class LAMMPS *, int, char **);
  virtual ~FixTemplateMultiplespheres();

  virtual void post_create();
  double max_r_bound() const;
  double max_rad() const;
  double min_rad() const;
  int maxtype();
  int mintype();
  int number_spheres();
  bool is_bonded() const { return bonded; }

  // multi insertion
  virtual void init_ptilist(int);
  void randomize_ptilist(int ,int );

  virtual void finalize_insertion();

  inline double get_bond_id(const int i) const
  { return fix_bond_random_id ? fix_bond_random_id->vector_atom[i] : 0.0; }

 protected:

  // number of partners and partner array
  int *np;
  int **p;

  // template calculations
  virtual void calc_bounding_sphere();
  virtual void calc_center_of_mass();

  // sqr distance from x_sphere[j] to xtest
  double dist_sqr(int j,double *xtest);

  // generate random point in bbox
  void generate_xtry(double *xtry);

  // number of spheres in template
  int nspheres;

  // coords of each sphere with respect to center of mass
  double **x_sphere;

  // radius of each sphere
  double *r_sphere;

  // scale factor if read from a file
  double scale_fact;

  // atom type might be variable if read from file
  int *atom_type_sphere;

  // bounding box
  double x_min[3], x_max[3];

  // bounding sphere - radius and coordinates with respect to com
  double r_bound;
  double x_bound[3];

  // radius of sphere with equal volume
  double r_equiv;

  // number of tries for mc
  int ntry;

  bool bonded;
  int bond_type;

  class FixPropertyAtom *fix_bond_random_id;
};

}

#endif
#endif
