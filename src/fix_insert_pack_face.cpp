/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2015-     JKU Linz

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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_insert_pack_face.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "region_mesh_hex.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_particledistribution_discrete_face.h"
#include "fix_template_sphere.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "particleToInsert.h"
#include "fix_multisphere.h"
#include "math_const.h"
#include "math_extra_liggghts.h"
#include "fix_massflow_mesh_face.h"

#define SEED_OFFSET 12

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixInsertPackFace::FixInsertPackFace(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while (iarg < narg && hasargs) {
    hasargs = false;
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      ins_region = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"ntry_mc") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      ntry_mc = atoi(arg[iarg+1]);
      if (ntry_mc < 1000) error->fix_error(FLERR,this,"ntry_mc must be > 1000");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"warn_region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      if (strcmp(arg[iarg+1],"yes") == 0)
        warn_region = true;
      else if (strcmp(arg[iarg+1],"no") == 0)
        warn_region = false;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after 'warn_region'");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"massflow_face") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0 || strcmp(modify->fix[ifix]->style,"massflow/mesh/face"))
        error->fix_error(FLERR,this,"Fix insert/pack/face requires you to define a valid ID for a fix of type massflow/mesh/face");
      massflowface = static_cast<FixMassflowMeshFace*>(modify->fix[ifix]);
      int n = strlen(arg[iarg+1]) + 1;
      idmassflowface = new char[n];
      strcpy(idmassflowface,arg[iarg+1]);

      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"cg") == 0) { // TODO: remove when cg and fg are separate simulations -> cg_ = force->cg();
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      cg_ = atof(arg[iarg+1]);
      cg3_ = cg_*cg_*cg_;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(style,"insert/pack/face") == 0) {
      error->fix_error(FLERR,this,"unknown keyword");
    }
  }

  // no fixed total number of particles inserted by this fix exists
  if (strcmp(style,"insert/pack/face") == 0)
    ninsert_exists = 0;
}

/* ---------------------------------------------------------------------- */

FixInsertPackFace::~FixInsertPackFace()
{
  delete [] idmassflowface;
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

void FixInsertPackFace::post_create()
{
  int nfaces = massflowface->get_face_ids_size();
  if (nfaces > 0) {
    fraction_face_local.resize(nfaces);
    min_face_extent_local.resize(nfaces);
    volume_face_absolut.resize(nfaces);
  }
  rounding_all.resize(comm->nprocs, 0.0);
  rounding_all[0] = 0.5;

  if (!fix_distribution) {
    char * id_name = new char[38+strlen(id)];
    sprintf(id_name,"insert_particle_distribution_per_face%s",id);

    char **fixarg = new char*[4];
    fixarg[0] = id_name;
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "particledistribution/discrete/face";
    fixarg[3] = (char *) "5531"; // TODO: set RNG seed dynamically

    modify->add_fix(7,fixarg);
    fix_distribution = static_cast<FixParticledistribution*>(modify->find_fix_id(id_name));
    delete []fixarg;
    delete []id_name;
  }

  FixInsert::post_create();
}

/* ---------------------------------------------------------------------- */

// do NOT call FixInsert::init_defaults() from here
// this would overwrite settings that were parsed in FixInsert constructor
// since init_defaults() is called from constructor of both classes, both
// FixInsert::init_defaults() and FixInsertPack::init_defaults() are called
// at the correct point

void FixInsertPackFace::init_defaults()
{
  massflowface = NULL;
  idmassflowface = NULL;
  ins_region = NULL;
  ins_region_mesh_hex = NULL;
  idregion = 0;
  ntry_mc = 100000;

  region_volume = region_volume_local = 0.;

  insertion_ratio = 0.;

  warn_region = true;
  cg_ = 1.0;
  cg3_ = 1.0;
}

/* ---------------------------------------------------------------------- */

void FixInsertPackFace::init()
{
  FixInsert::init();

  if (ins_region) {
    int iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->fix_error(FLERR,this,"region ID does not exist");
    ins_region = domain->regions[iregion];
  }

  if (strcmp(ins_region->style, "mesh/hex") == 0) {
    ins_region_mesh_hex = static_cast<RegHexMesh*>(ins_region);
  }
}

/* ----------------------------------------------------------------------
   perform error checks
------------------------------------------------------------------------- */

void FixInsertPackFace::calc_insertion_properties()
{
  // error check on region
  if (!ins_region)
    error->fix_error(FLERR,this,"must define an insertion region");
  ins_region->reset_random(seed + SEED_OFFSET);

  calc_region_volume_local();
  if (region_volume <= 0. || region_volume_local < 0. || (region_volume_local - region_volume)/region_volume > 1e-3 )
    error->one(FLERR,"Fix insert: Region volume calculation with MC failed");

  if (ins_region->dynamic_check())
    error->fix_error(FLERR,this,"dynamic regions are not allowed");

  // error check on insert_every
  if (insert_every < 0)
    error->fix_error(FLERR,this,"must define 'insert_every'");

  // error checks to disallow args from FixInsert
  if (ninsert > 0 || massinsert > 0.)
    error->fix_error(FLERR,this,"specifying 'nparticles' or 'mass' not allowed");
  if (nflowrate > 0. || massflowrate > 0.)
    error->fix_error(FLERR,this,"specifying 'nflowrate' or 'massflowrate' not allowed");
}

/* ----------------------------------------------------------------------
   calculate volume of region on my subbox
   has to be called at initialization and before every insertion in case
   box is changing
------------------------------------------------------------------------- */

void FixInsertPackFace::calc_region_volume_local()
{
  ins_region->volume_mc(ntry_mc,all_in_flag==0?false:true,fix_distribution->max_r_bound(),
                        region_volume,region_volume_local);
}

/* ----------------------------------------------------------------------
   number of particles to insert this timestep
   depends on number of particles in region already
------------------------------------------------------------------------- */

int FixInsertPackFace::calc_ninsert_this()
{
  int ninsert_this = 0;

  // check if region extends outside simulation box
  // if so, throw error if boundary setting is "f f f"

  if (warn_region && ins_region->bbox_extends_outside_box()) {
    for (int idim = 0; idim < 3; idim++)
      for (int iface = 0; iface < 2; iface++)
        if (domain->boundary[idim][iface] == 1)
          error->fix_error(FLERR,this,"Insertion region extends outside simulation box and a fixed boundary is used. "
                            "Please use non-fixed boundaries in this case only");
  }

  int ifix = modify->find_fix(idmassflowface);
  if (ifix < 0 || strcmp(modify->fix[ifix]->style,"massflow/mesh/face")) {
    massflowface = NULL;
    //error->fix_error(FLERR,this,"Fix insert/pack/face requires you to define a valid ID for a fix of type massflow/mesh/face");
  } else {
    massflowface = static_cast<FixMassflowMeshFace*>(modify->fix[ifix]);
  }

  if (massflowface) {
      ninsert_this = static_cast<int>(massflowface->compute_vector(6));
  }
  insertion_ratio = 0.;

  // can be < 0 due to overflow, round-off etc
  if (ninsert_this < -200000)
    error->fix_error(FLERR,this,"overflow in particle number calculation: inserting too many particles in one step");
  if (ninsert_this < 0) ninsert_this = 0;

  if (     insertion_ratio < 0.) insertion_ratio = 0.;
  else if (insertion_ratio > 1.) insertion_ratio = 1.;

  return ninsert_this;
}

/* ---------------------------------------------------------------------- */

double FixInsertPackFace::insertion_fraction()
{
  // have to re-calculate region_volume_local in case simulation box is changing
  if (domain->box_change)
    calc_region_volume_local();

  const std::map<int,int>& faceids = massflowface->get_face_ids();
  for (std::map<int,int>::const_iterator it=faceids.begin(); it!=faceids.end(); ++it) {
    fraction_face_local[it->second] = insertion_fraction_face(it->first);
  }

  // mc method may be too vague -> normalize fraction
  int nfaces = massflowface->get_face_ids_size();
  std::vector<double> fraction_face_all(nfaces, 1.0);
  MPI_Allreduce(&fraction_face_local[0], &fraction_face_all[0], nfaces, MPI_DOUBLE, MPI_SUM, world);

  for (int i=0; i<nfaces; ++i) {
    fraction_face_local[i] /= fraction_face_all[i];
  }

    return region_volume_local/region_volume;
}

/* ---------------------------------------------------------------------- */

double FixInsertPackFace::insertion_fraction_face(int face_id)
{
  // check face_id
  if (massflowface->get_face_ids().find(face_id) == massflowface->get_face_ids().end())
    error->fix_error(FLERR,this,"invalid face id");

  // generate random points in hexahedral cell (hex_randpos), then test how many of them are in the local domain
  // do not: generate random points in domain and test if they are inside the hexahedral cell
  // since point in cell test may run into numerical problems!

  int inside = 0;
  int iHex = ins_region_mesh_hex->get_hex("face_id", face_id);
  const int index = massflowface->get_face_ids().at(face_id);
  min_face_extent_local[index] = 0.;
  volume_face_absolut[index] = 0.;
  if (iHex >= 0 && iHex < ins_region_mesh_hex->n_hex()) {
    volume_face_absolut[index] = ins_region_mesh_hex->hex_vol(iHex);

    double bounds[6] = {};
    ins_region_mesh_hex->hex_bounds(iHex, bounds);
    BoundingBox hexbb(bounds);
    hexbb.shrinkToSubbox(domain->sublo,domain->subhi);
    double extent[3] = {};
    hexbb.getExtent(extent);
    min_face_extent_local[index] = std::min(extent[0], std::min(extent[1], extent[2]));

    double pos[3] = {};
    for (int i = 0; i < ntry_mc; ++i) {
      ins_region_mesh_hex->hex_randpos(iHex, pos);

      if (!domain->is_in_subdomain(pos))
        continue;

      // TEST cut???
      ++inside;
    }
  }

  return static_cast<double>(inside)/static_cast<double>(ntry_mc);
}

/* ----------------------------------------------------------------------
   distribute insertions across processors
------------------------------------------------------------------------- */

int FixInsertPackFace::distribute_ninsert_this(int ninsert_this)
{
  int me, nprocs, ninsert_this_local=0;

  me = comm->me;
  nprocs = comm->nprocs;

  // for exact_number==1 have to allgather to exactly match ninsert_this

  int nfaces = fraction_face_local.size();
  double *fraction_face_local_all = NULL;
  double *min_face_extent_local_all = NULL;
  MPI_Allgather_Vector(&fraction_face_local[0], nfaces, fraction_face_local_all, world);
  MPI_Allgather_Vector(&min_face_extent_local[0], nfaces, min_face_extent_local_all, world);

  // all procs calculate local insertion distribution because:
  // a) synchronization of distribution data would be pretty complicated
  // b) if only proc 0 does the work, all other procs would be idle/waiting anyway
  // we just have to make sure that all procs arrive at the same distribution in the end

  const std::vector<DiscreteParticleDistribution>& distributions_face = massflowface->get_distributions();

  // NOTE: fraction_face_local not sorted by face ID! --> fraction_face_local[faceid2index_[face_id]];
  //       distributions_face not sorted by face ID!
  std::vector<std::vector<std::vector<int> > > distributions_face_local_all(nprocs); //procs - faces - distributions:#particles
  for (int iproc = 0; iproc < nprocs; ++iproc) {
    distributions_face_local_all[iproc].resize(nfaces);
  }

  // loop over faces
  for (int iface=0; iface<nfaces; ++iface) {
    std::vector<double> accumulatedParticleVolume_facei(nprocs, 0.);

    // loop over templates
    DiscreteParticleDistribution::const_iterator it_dist = distributions_face[iface].begin();
    for (int idist=0; it_dist!=distributions_face[iface].end(); ++it_dist, ++idist) {
      std::vector<int> accumulatedTemplateParticles_facei(distributions_face[iface].size(), 0);

      // TODO: scale number of particles by cg3_ factor!!! (if this is not the finest grained level)
      double nparticles = it_dist->second;
      double diameter = 2. * it_dist->first.radius_;
      double volume_single = (1./6.)*MY_PI*diameter*diameter*diameter;

      for (int iproc = 0; iproc < nprocs; ++iproc) {
        double procvolfracface = fraction_face_local_all[iproc*nfaces+iface]; // local volume fraction
        double procminfaceextent = min_face_extent_local_all[iproc*nfaces+iface];
        // round to integer number of particles
        // rounding_all tries to avoid that first proc is preferred for insertion
        int nparts = diameter < procminfaceextent? static_cast<int>(procvolfracface*nparticles + rounding_all[iproc]) : 0;
        while(nparts > 0 && nparts + accumulatedTemplateParticles_facei[idist] > nparticles) {
          --nparts; // just in case some wierd rounding issues occured
        }
        double volume = nparts*volume_single;
        while (nparts > 0 && volume + accumulatedParticleVolume_facei[iproc] > procvolfracface*volume_face_absolut[iface]) {
          // volume of selected particles larger than local fraction of cell!
          // that may happen if e.g. local cell volume fraction is smaller than larger particles
          // but large enough to hold a couple of small particles
          volume -= volume_single;
          --nparts;
        }
        accumulatedParticleVolume_facei[iproc] += volume; // volume of particles selected to insert locally this far
        accumulatedTemplateParticles_facei[idist] += nparts;
        distributions_face_local_all[iproc][iface].push_back(nparts);
        if (iproc == me)
          ninsert_this_local += nparts;
      }

      // now fix up number distribution to add up to original sum
      {
        // is there a mismatch we need to correct?
        while (nparticles - static_cast<double>(accumulatedTemplateParticles_facei[idist]) > 0.5) {
          // get a list of potential proc candidates for additional insertion
          // sorted by their mismatch between desired and assigned fraction of particles
          std::map<double, int> potential_procs;
          for (int iproc = 0; iproc < nprocs; ++iproc) {
            double procvolfracface = fraction_face_local_all[iproc*nfaces+iface];
            double procminfaceextent = min_face_extent_local_all[iproc*nfaces+iface];
            // is there enough local space?
            if (diameter < procminfaceextent &&
                accumulatedParticleVolume_facei[iproc] + volume_single < procvolfracface * volume_face_absolut[iface]) {
                double desired = procvolfracface*nparticles;
                double assigned = static_cast<double>(distributions_face_local_all[iproc][iface].back());
                double mismatch = desired - assigned;
                potential_procs[mismatch] = iproc;
            }
          }

          if (potential_procs.empty()) {
            // failed to correct particle distribution: not enough local volume!
            break;
          }

          int correcting_proc = potential_procs.rbegin()->second;
          distributions_face_local_all[correcting_proc][iface].back()++;
          if (correcting_proc == me)
            ninsert_this_local++;
          accumulatedTemplateParticles_facei[idist]++;
        }
      }
    }
  }

  // send result to each proc
  FixParticledistributionDiscreteFace *fix_pddf = (FixParticledistributionDiscreteFace*)fix_distribution;
  fix_pddf->set_distribution_local(massflowface, distributions_face_local_all[me]);
  maxrad = std::max(maxrad, fix_pddf->max_rad());

  delete [] fraction_face_local_all;
  delete [] min_face_extent_local_all;

  rounding_all.push_front(rounding_all.back());
  rounding_all.pop_back();

   return ninsert_this_local;
}


/* ---------------------------------------------------------------------- */

inline int FixInsertPackFace::is_nearby(int i)
{
  double pos[3], rad, cut;

  vectorCopy3D(atom->x[i],pos);
  rad = atom->radius[i];

  // choose right distance depending on all_in_flag

  if (all_in_flag) cut = maxrad;
  else cut = rad + maxrad;

  if (ins_region->match_expandby_cut(pos,cut))
    return 1;

  return 0;
}

/* ---------------------------------------------------------------------- */

BoundingBox FixInsertPackFace::getBoundingBox() const
{
  BoundingBox bb(ins_region->extent_xlo, ins_region->extent_xhi,
                 ins_region->extent_ylo, ins_region->extent_yhi,
                 ins_region->extent_zlo, ins_region->extent_zhi);

  const double cut = 2*maxrad;
  bb.extendByDelta(cut);
  bb.shrinkToSubbox(domain->sublo, domain->subhi);

  return bb;
}

/* ----------------------------------------------------------------------
   calc # of maximum tries
   propertional to total desired # of particles to insert on this
   subdomain to ensure insertion "does not give up too early" if a low
   remaining # of insertions is to be performed
------------------------------------------------------------------------- */

int FixInsertPackFace::calc_maxtry(int ninsert_this_local)
{
  if (insertion_ratio >= 1.) return ninsert_this_local * maxattempt;
  else return static_cast<int>(static_cast<double>(ninsert_this_local*maxattempt) / (1.-insertion_ratio));
}

/* ----------------------------------------------------------------------
   generate random positions within insertion volume
   perform overlap check via xnear
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertPackFace::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
  ninserted_this_local = ninserted_spheres_this_local = 0;
  mass_inserted_this_local = 0.;

  double pos[3];
  ParticleToInsert *pti;

  int ntry = 0;
  // int maxtry = calc_maxtry(ninsert_this_local);

  double v_toInsert[3];
  vectorZeroize3D(v_toInsert);

  // always do overlap check
  {
    int ifix = modify->find_fix(idmassflowface);
    if (ifix < 0 || strcmp(modify->fix[ifix]->style,"massflow/mesh/face"))
      return;

    massflowface = static_cast<FixMassflowMeshFace*>(modify->fix[ifix]);
    if (!massflowface)
      error->fix_error(FLERR, this, "Failed to locate fix massflow/mesh/face");

    const std::map<int,int>& faceids = massflowface->get_face_ids();

    int iface = 0;
    FixParticledistributionDiscreteFace *fix_pddf = (FixParticledistributionDiscreteFace*)fix_distribution;
    for (std::map<int,int>::const_iterator it=faceids.begin(); it!=faceids.end(); ++it, ++iface) {
      int particle_this_face = 0;
      int particle_this_face_local = fix_pddf->pti_list_face_local[iface].size();
      int ninserted_this_face_local = 0;
      if (massflowface) {
        if (cg_ > 1.0) {
          double particle_this_face_real = massflowface->compute_array_by_id(it->first, 6) + remaining_ptis[it->first];
          double remainder = fmod(particle_this_face_real, cg3_);
          particle_this_face = static_cast<int>(particle_this_face_real / cg3_);
          if ((remainder/cg3_) > 0.5) {
            remainder -= cg3_;
            ++particle_this_face;
          }
          remaining_ptis[it->first] = remainder;
        } else {
          particle_this_face = massflowface->compute_array_by_id(it->first, 6);
        }
        v_insert[0] = massflowface->compute_array_by_id(it->first, 7);
        v_insert[1] = massflowface->compute_array_by_id(it->first, 8);
        v_insert[2] = massflowface->compute_array_by_id(it->first, 9);
      }

      ntry = 0;

      int maxtry = maxattempt*particle_this_face_local;
      while (ntry < maxtry && ninserted_this_face_local < particle_this_face_local) {
        pti = fix_pddf->pti_list_face_local[iface][ninserted_this_face_local];
        double rbound = pti->r_bound_ins;

        int nins = 0;
        while (nins == 0 && ntry < maxattempt) {
          do {
            // generate a point in my subdomain
            int iHex = ins_region_mesh_hex->get_hex("face_id", it->first);
            if (iHex >= 0 && iHex < ins_region_mesh_hex->n_hex()) {
              do {
                do {
                  ins_region_mesh_hex->hex_randpos(iHex, pos);
                } while (!domain->is_in_subdomain(pos));
              } while (all_in_flag && ins_region_mesh_hex->match_hex_cut(iHex,pos,rbound));
            } else {
              error->fix_error(FLERR, this, "Failed to find hex cell by face_id!");
            }
            ++ntry;
          } while (ntry < maxattempt /* *particle_this_face_local*/ && domain->dist_subbox_borders(pos) < rbound);

          // randomize vel, omega, quat here
          vectorCopy3D(v_insert,v_toInsert);

          // could randomize vel, omega, quat here
          generate_random_velocity(v_toInsert);

          if (quat_random_)
            MathExtraLiggghts::random_unit_quat(random,quat_insert);

          nins = pti->check_near_set_x_v_omega(pos,v_toInsert,omega_insert,quat_insert,neighList);
        }

        if (nins > 0) {
          ninserted_spheres_this_local += nins;
          mass_inserted_this_local += pti->mass_ins;
          ninserted_this_local++;
          ninserted_this_face_local++;
        } else {
          fix_pddf->pti_list_face_local[iface].erase(fix_pddf->pti_list_face_local[iface].begin()+ninserted_this_face_local);
          --particle_this_face_local;
          // failed to insert particle on this proc
        }
      }
    }
  }
}



/* ---------------------------------------------------------------------- */

void FixInsertPackFace::restart(char *buf)
{
  FixInsert::restart(buf);

  ins_region->reset_random(seed + SEED_OFFSET);
}
