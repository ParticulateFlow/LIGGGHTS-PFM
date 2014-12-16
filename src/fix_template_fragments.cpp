/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2014- JKU Linz

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

#include "fix_template_fragments.h"
#include "atom.h"
#include "atom_vec.h"
#include "stdlib.h"
#include "vector_liggghts.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "particleToInsert_fragments.h"
#include "pair_gran.h"
#include "math_const.h"
#include "particleSizeDistribution.h"
#include "particleSpatialDistribution.h"
#include "primitive_wall.h"
#include <map>

using namespace MathConst;


using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;

/*NL*/#define LMP_DEBUGMODE_FRAGMENTS true

/* ---------------------------------------------------------------------- */

FixTemplateFragments::FixTemplateFragments(LAMMPS *lmp, int narg, char **arg) :
  FixTemplateSphere(lmp, narg, arg)
{
  if (pdf_density->rand_style() != RANDOM_CONSTANT)
    error->all(FLERR,"Fix particletemplate/fragments currently only supports constant density");

  if (pdf_radius)
    error->fix_error(FLERR,this,"currently does not support keyword 'radius'");

  bool hasargs = true;
  while (iarg < narg && hasargs) {
    hasargs = false;
    if (strcmp(arg[iarg],"breakage_index") == 0) {
      if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments for breakage_index");
      ++iarg;
      t10_max = atof(arg[iarg++]);
      rad_min_pct = atof(arg[iarg++]); // percentage of original size
      hasargs = true;
    } else if (strcmp(arg[iarg],"tn_family") == 0) {
      if (iarg+3 > narg) error->fix_error(FLERR,this,"not enough arguments for tn_family");
      int sizes = atoi(arg[iarg+1]);
      iarg += 2;
      if (iarg+sizes > narg) error->fix_error(FLERR,this,"not enough arguments for tn_family");
      while (sizes > 0) {
        radiiMassFractions[atoi(arg[iarg++])] = 0.0;
        --sizes;
      }
      hasargs = true;
    } else {
      std::string unknown("unknown keyword ");
      unknown += arg[iarg];
      error->fix_error(FLERR,this,unknown.c_str());
    }
  }

  // set default tn_family if none is specified
  if(radiiMassFractions.empty()) {
    radiiMassFractions[2] = 0.0;
    radiiMassFractions[3] = 0.0;
    radiiMassFractions[5] = 0.0;
    radiiMassFractions[10] = 0.0;
  }

  nspheres = 0;
  x_sphere = NULL;
  r_sphere = NULL;
}

/* ---------------------------------------------------------------------- */

FixTemplateFragments::~FixTemplateFragments()
{
  memory->destroy(x_sphere);
  delete [] r_sphere;
}

/* ---------------------------------------------------------------------- */

void FixTemplateFragments::post_create()
{
  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  max_type = pair_gran->mpg->max_type();
  dnum = pair_gran->dnum();

  // don't know these in advance ...
  volume_expect = 1.0;
  mass_expect = volume_expect*expectancy(pdf_density);
  r_equiv = 1.0;
  r_bound = 0.0;
}

/* ----------------------------------------------------------------------*/

double FixTemplateFragments::min_rad()
{
  return 0.0;
}

/* ----------------------------------------------------------------------*/

double FixTemplateFragments::max_rad()
{
  return 0.0;
}

/* ----------------------------------------------------------------------*/

double FixTemplateFragments::max_r_bound()
{
  return r_bound;
}

/* ----------------------------------------------------------------------*/
// used by pdd to check which template has the most spheres
int FixTemplateFragments::number_spheres()
{
  if (LMP_DEBUGMODE_FRAGMENTS) fprintf(screen,"FixTemplateFragments::number_spheres: nspheres = %d\n", nspheres);
  return nspheres;
}

/* ----------------------------------------------------------------------*/

void FixTemplateFragments::randomize_single()
{
  error->all(FLERR,"re-implement FixTemplateFragments::randomize_single()!");
}

/* ----------------------------------------------------------------------*/

void FixTemplateFragments::init_ptilist(int n_random_max)
{
  // n_random_max is local
  if (pti_list) error->all(FLERR,"invalid FixTemplateSphere::init_list()");
  n_pti_max = n_random_max;
  pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsertFragments*),"pti_list");
  memset(pti_list,0,n_pti_max*sizeof(ParticleToInsertFragments*));
}

/* ----------------------------------------------------------------------*/

void FixTemplateFragments::randomize_ptilist(int n_random,int distribution_groupbit)
{
  for (int i = 0; i < n_random; ++i) {

    nspheres = r_sphere_list[i].size();
    delete pti_list[i];
    pti_list[i] = new ParticleToInsertFragments(lmp,nspheres);

    ParticleToInsertFragments *pti = static_cast<ParticleToInsertFragments*>(pti_list[i]);

    pti->density_ins = expectancy(pdf_density);
    pti->volume_ins = 0.0;
    pti->volume_ins = 0.0;
    pti->collision_factor = collision_factor[i];
    pti->r_bound_ins = r_bound;
    pti->atom_type = atom_type;

    for (int j = 0; j < nspheres; ++j) {
      pti->radius_ins[j] = r_sphere_list[i][j];
      pti->x_ins[j][0] = x_sphere_list[i][j][0];
      pti->x_ins[j][1] = x_sphere_list[i][j][1];
      pti->x_ins[j][2] = x_sphere_list[i][j][2];
      pti->volume_ins += MY_4PI3 * r_sphere_list[i][j] * r_sphere_list[i][j] * r_sphere_list[i][j];
    }

    pti->mass_ins = pti->volume_ins * pti->density_ins;

    vectorZeroize3D(pti->v_ins);
    vectorZeroize3D(pti->omega_ins);

    pti->groupbit = groupbit | distribution_groupbit; //NP also contains insert_groupbit
  }
}


void FixTemplateFragments::pre_insert(
    int n_total,
    double **breakdata,
    const std::multimap<int,std::vector<double> > &contacting_atoms_mm,
    const std::multimap<int,PrimitiveWall*> &prim_walls_mm,
    const std::multimap<int,TriMeshContacts*> &meshes_mm)
{
  // generate particle size/spatial distributions
  r_sphere_list.clear();
  collision_factor.clear();
  r_bound = 0.0;
  nspheres = 0;

  if (n_total <= 0) {
    x_sphere_list.clear();
    return;
  }

  // breakdata[i][6] ... radius
  double density = expectancy(pdf_density);
  x_sphere_list.resize(n_total);

  for (int i = 0; i < n_total; ++i) {
    ParticleSizeDistribution psd(breakdata[i][7], density, breakdata[i][6], breakdata[i][6]*rad_min_pct, breakdata[i][6]-0.5*breakdata[i][10], t10_max);
    psd.range_mass_fractions(radiiMassFractions);
    std::vector<double> radii;
    psd.radii(radiiMassFractions, radii);
    r_sphere_list.push_back(radii);

    if (r_bound < breakdata[i][6]) {
      r_bound = breakdata[i][6];
    }
    if (nspheres < static_cast<int>(radii.size())) {
      nspheres = static_cast<int>(radii.size());
    }

    // overlapping atoms
    std::vector<std::vector<double> > ext_atoms;
    {
      std::pair<std::multimap<int,std::vector<double> >::const_iterator, std::multimap<int,std::vector<double> >::const_iterator> ret;
      ret = contacting_atoms_mm.equal_range(static_cast<int>(breakdata[i][8]));
      for (std::multimap<int,std::vector<double> >::const_iterator it = ret.first; it != ret.second; ++it) {
        ext_atoms.push_back(it->second);
      }
    }

    // overlapping primitive walls
    std::vector<PrimitiveWall*> prim_walls;
    {
      std::pair<std::multimap<int,PrimitiveWall*>::const_iterator, std::multimap<int,PrimitiveWall*>::const_iterator> ret;
      ret = prim_walls_mm.equal_range(static_cast<int>(breakdata[i][8]));
      for (std::multimap<int,PrimitiveWall*>::const_iterator it = ret.first; it != ret.second; ++it) {
        prim_walls.push_back(it->second);
      }
    }

    // overlapping meshes
    std::vector<TriMeshContacts*> meshes;
    {
      std::pair<std::multimap<int,TriMeshContacts*>::const_iterator, std::multimap<int,TriMeshContacts*>::const_iterator> ret;
      ret = meshes_mm.equal_range(static_cast<int>(breakdata[i][8]));
      for (std::multimap<int,TriMeshContacts*>::const_iterator it = ret.first; it != ret.second; ++it) {
        meshes.push_back(it->second);
      }
    }

    ParticleSpatialDistribution pxd(random, breakdata[i][10]);
    pxd.randomInsertion(&breakdata[i][0], breakdata[i][6], radii, x_sphere_list[i], ext_atoms, prim_walls, meshes);

    double energy = elastic_energy(r_sphere_list[i], x_sphere_list[i]);
    double CF = breakdata[i][9]/energy;
    CF = cbrt(CF*CF);
    collision_factor.push_back(CF);
  }
}


// geometrically correct elastic energy
double FixTemplateFragments::elastic_energy(const std::vector<double> &radius, const std::vector<std::vector<double> > &x)
{
  const unsigned int fragments = radius.size();
  double energy = 0.0;
  int overlaps = 0;

  if(fragments > 0) {
    const int64_t normalmodel = pair_gran->hashcode() & 0x0f;
    const double density = expectancy(pdf_density);
    const double *Y, *nu;
    const int itype = atom_type - 1;
    Y  = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0,style))->get_values();
    nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0,style))->get_values();
    const double Yeff = Y[itype]/( 2. * (1. - nu[itype]*nu[itype]) ); // fragments have same type

    for(unsigned int pi = 0; pi < fragments-1; ++pi) {
      for(unsigned int pj = pi+1; pj < fragments; ++pj) {
        const double delx = x[pj][0] - x[pi][0];
        const double dely = x[pj][1] - x[pi][1];
        const double delz = x[pj][2] - x[pi][2];
        const double rsq = delx * delx + dely * dely + delz * delz;
        const double radsum = radius[pi] + radius[pj];

        if (rsq < radsum * radsum) {
          ++overlaps;
          const double rij = sqrt(rsq);
          const double deltan = radsum - rij;
          double reff = radius[pi] * radius[pj] / radsum;
          double kn = 0.0;

          if (normalmodel == 6) { // hertz/break
            kn = 4./3.*Yeff*sqrt(reff * deltan);
            energy += 0.4 * kn * deltan * deltan; // 0.4 = (2./5.)
          } else if (false) { // hooke/break
            const double sqrtval = sqrt(reff);
            const double mi = MY_4PI3 * radius[pi] * radius[pi] * radius[pi] * density;
            const double mj = MY_4PI3 * radius[pj] * radius[pj] * radius[pj] * density;
            const double meff = mi * mj / (mi + mj);
            const double charVel = static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0,style))->compute_scalar();
            kn = (16./15.) * sqrtval * Yeff * pow(15. * meff * charVel * charVel /(16. * sqrtval * Yeff), 0.2);
            energy += 0.5 * kn * deltan * deltan;
          }
        }
      }
    }
  }

  return energy;
}

