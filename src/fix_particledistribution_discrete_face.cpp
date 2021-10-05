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

#include "fix_particledistribution_discrete_face.h"
#include "atom.h"
#include "modify.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "comm.h"
#include "container.h"
#include "fix_massflow_mesh_face.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscreteFace::FixParticledistributionDiscreteFace(LAMMPS *lmp, int narg, char **arg) :
  FixParticledistribution(lmp, narg, arg)
{
  restart_global = 1;

  if(narg < 4)
    error->all(FLERR,"Illegal fix particledistribution/discrete/face command, not enough arguments");
  seed = atoi(arg[3]) + comm->me;
  random = new RanPark(lmp,seed);

  mintype = 1;
  maxtype = atom->ntypes;

  maxrad = maxrbound = 0.0;
  n_pti_max = 0;
}

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscreteFace::~FixParticledistributionDiscreteFace()
{
  delete_pti_list_face_local();
}

/* ----------------------------------------------------------------------*/

void FixParticledistributionDiscreteFace::delete_pti_list_face_local()
{
  if(!pti_list_face_local.empty())
  {
    std::vector<pti_list_type>::iterator it_face = pti_list_face_local.begin();
    for(; it_face != pti_list_face_local.end(); it_face++)
    {
      pti_list_type::iterator it_pti = it_face->begin();
      for(; it_pti != it_face->end(); it_pti++)
      {
        delete *it_pti;
      }
    }
    pti_list_face_local.clear();
  }
}

/* ----------------------------------------------------------------------*/

void FixParticledistributionDiscreteFace::set_distribution_local(const std::vector<DiscreteParticleDistribution>& distributions, const std::vector<std::vector<int> > & distributions_face_local, double cg, int type_offset)
{
  // TODO: cg = force->cg() can be used when differently resolved levels are separate simulations -> remove parameter
  //       also, type_offset will be 0 in that case
  delete_pti_list_face_local();

  n_pti_max = 0;
  volexpect=0.;
  massexpect=0.;
  double masstotal = 0.;
  maxrad = maxrbound = 0.;
  minrad = 1000.;
  maxnspheres = 1;

  if(distributions.empty()) return;

  int n_face_ids = distributions.size();

  pti_list_face_local.resize(n_face_ids);
  for(int iface=0; iface<n_face_ids; ++iface)
  {
    DiscreteParticleDistribution::const_iterator it_dist = distributions[iface].begin();
    for(int idist=0; it_dist!=distributions[iface].end(); ++it_dist, ++idist)
    {
      int type = it_dist->first.atomtype_ + type_offset;

      double radius = it_dist->first.radius_ * cg;
      if(radius > maxrad) maxrad = radius;
      else if(radius < minrad)  minrad = radius;

      int ntemplates_to_insert = distributions_face_local[iface][idist];
      for(int itemplate=0; itemplate<ntemplates_to_insert; ++itemplate)
      {
        ParticleToInsert *pti = new ParticleToInsert(lmp);
        pti->atom_type = type;
        pti->radius_ins[0] =  pti->r_bound_ins = radius;
        pti->density_ins = it_dist->first.density_;
        pti->volume_ins = radius * radius * radius * 4.*M_PI/3.;
        pti->mass_ins = it_dist->first.mass_ * cg * cg *cg;
        // init insertion position
        vectorZeroize3D(pti->x_ins[0]);
        vectorZeroize3D(pti->v_ins);
        vectorZeroize3D(pti->omega_ins);
        pti->groupbit = groupbit;
        pti_list_face_local[iface].push_back(pti);
        ++n_pti_max;
        volexpect  += pti->mass_ins*pti->volume_ins;
        massexpect += pti->mass_ins*pti->mass_ins;
        masstotal  += pti->mass_ins;
      }
    }
  }

  if(masstotal > 0.)
  {
    volexpect  /= masstotal;
    massexpect /= masstotal;
  }

  maxrbound = maxrad;
}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_list() command
   typically only called once before first insertion step
------------------------------------------------------------------------- */

void FixParticledistributionDiscreteFace::random_init_list(int ntotal)
{
}

/* ----------------------------------------------------------------------
   returns number of particles to be inserted.
   typically called once per insertion step
------------------------------------------------------------------------- */

int FixParticledistributionDiscreteFace::randomize_list(int ntotal,int insert_groupbit,int exact_number)
{
  if(ntotal > n_pti_max)
  {
    error->one(FLERR,"Faulty implementation: FixParticledistributionDiscreteFace::randomize_list() called for more particles than defined in random_init_list()");
  }

  groupbit = insert_groupbit;
  ninsert = n_pti_max;
  ninserted = ninsert;
  return ninsert;
}

/* ----------------------------------------------------------------------
   preparations before insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscreteFace::pre_insert(int n, FixPropertyAtom *fp, double val, int idx, int ival, int iidx)
{
  FixParticledistribution::pre_insert();

  // set fix property as desired by fix insert
  if(fp)
  {
    std::vector<pti_list_type>::iterator it_face = pti_list_face_local.begin();
    for(; it_face!=pti_list_face_local.end(); it_face++)
    {
      pti_list_type::iterator it_pti = it_face->begin();
      for(; it_pti!=it_face->end(); it_pti++)
      {
        (*it_pti)->fix_properties.push_back(fp);
        std::vector<double> value(1,val);
        (*it_pti)->fix_property_values.push_back(value);
      }
    }
  }
  else if(idx >= 0)
  {
    std::vector<pti_list_type>::iterator it_face = pti_list_face_local.begin();
    for(; it_face!=pti_list_face_local.end(); it_face++)
    {
      pti_list_type::iterator it_pti = it_face->begin();
      for(; it_pti!=it_face->end(); it_pti++)
      {
        (*it_pti)->property_index = idx;
        (*it_pti)->fix_property_dvalue = val;
      }
    }
  }
  else if(iidx >= 0)
  {
    std::vector<pti_list_type>::iterator it_face = pti_list_face_local.begin();
    for(; it_face!=pti_list_face_local.end(); it_face++)
    {
      pti_list_type::iterator it_pti = it_face->begin();
      for(; it_pti!=it_face->end(); it_pti++)
      {
        (*it_pti)->property_iindex = iidx;
        (*it_pti)->fix_property_ivalue = ival;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   set particle properties - only pti needs to know which properties to set
   loop over all particles that have been inserted
------------------------------------------------------------------------- */

int FixParticledistributionDiscreteFace::insert(int n)
{
  int ninserted_spheres_local = 0;

  std::vector<pti_list_type>::iterator it_face = pti_list_face_local.begin();
  for(; it_face!=pti_list_face_local.end(); it_face++)
  {
    pti_list_type::iterator it_pti = it_face->begin();
    for(; it_pti!=it_face->end(); it_pti++)
    {
       ninserted_spheres_local += (*it_pti)->insert();
    }
  }

  return n_pti_max;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscreteFace::min_rad(int type) const
{
  // get minrad
  double minrad_type = 1000.;

  std::vector<pti_list_type>::const_iterator it_face = pti_list_face_local.begin();
  for(; it_face!=pti_list_face_local.end(); it_face++)
  {
    pti_list_type::const_iterator it_pti = it_face->begin();
    for(; it_pti!=it_face->end(); it_pti++)
    {
      if((*it_pti)->atom_type != type)
        continue;
      if((*it_pti)->radius_ins[0] < minrad_type)
        minrad_type = (*it_pti)->radius_ins[0];
    }
  }

  return minrad_type;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscreteFace::max_rad(int type) const
{
  // get maxrad
  double maxrad_type = 0.;

  std::vector<pti_list_type>::const_iterator it_face = pti_list_face_local.begin();
  for(; it_face!=pti_list_face_local.end(); it_face++)
  {
    pti_list_type::const_iterator it_pti = it_face->begin();
    for(; it_pti!=it_face->end(); it_pti++)
    {
      if((*it_pti)->atom_type != type)
        continue;
      if((*it_pti)->radius_ins[0] > maxrad_type)
        maxrad_type = (*it_pti)->radius_ins[0];
    }
  }

  return maxrad_type;
}

