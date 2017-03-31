/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2017-     JKU Linz

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_insert_pack_face_universe.h"
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
#include "universe.h"

#define SEED_OFFSET 12

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtraLiggghts;

/* ---------------------------------------------------------------------- */

FixInsertPackFaceUniverse::FixInsertPackFaceUniverse(LAMMPS *lmp, int narg, char **arg) :
  FixInsertPackFace(lmp, narg, arg)
{
  if (universe->nworlds == 1)
    error->all(FLERR,"Must have more than one processor partition for fix insert/pack/face/universe");

  cg_ = force->cg();
  cg3_ = cg_*cg_*cg_;
  receive_from_world_ = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"receive_from_partition") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      receive_from_world_ = atoi(arg[iarg+1])-1;
      if(receive_from_world_ < 0 || receive_from_world_ >= universe->nworlds)
        error->fix_error(FLERR,this,"receive_from_world_ must be a valid communicator");
      iarg += 2;
      break;
    }
    ++iarg;
  }
  idmassflowface_hash = JSHash(idmassflowface, mpi_tag_upper_bound(universe->uworld));
  // no fixed total number of particles inserted by this fix exists
  if (strcmp(style,"insert/pack/face/universe") == 0)
    ninsert_exists = 0;
}

/* ---------------------------------------------------------------------- */

FixInsertPackFaceUniverse::~FixInsertPackFaceUniverse()
{
}

/* ---------------------------------------------------------------------- */

void FixInsertPackFaceUniverse::post_create_per_face_data()
{
  int nfaces = 0;
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&nfaces, 1, MPI_INT, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&nfaces, 1, MPI_INT, 0, world);

  std::vector<int> faceids(nfaces,0);
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&faceids[0], nfaces, MPI_INT, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&faceids[0], nfaces, MPI_INT, 0, world);

  for(int i=0; i<nfaces; ++i)
    faceid2index_[faceids[i]] = i;

  if (nfaces > 0) {
    nparticles_face_.resize(nfaces, 0.);
    average_vx_face_out_.resize(nfaces, 0.);
    average_vy_face_out_.resize(nfaces, 0.);
    average_vz_face_out_.resize(nfaces, 0.);
    average_omegax_face_out_.resize(nfaces, 0.);
    average_omegay_face_out_.resize(nfaces, 0.);
    average_omegaz_face_out_.resize(nfaces, 0.);

    fraction_face_local.resize(nfaces);
    min_face_extent_local.resize(nfaces);
    volume_face_absolut.resize(nfaces);
  }
}

/* ----------------------------------------------------------------------
   number of particles to insert this timestep
------------------------------------------------------------------------- */
void FixInsertPackFaceUniverse::receive_ninsert_this(int& ninsert_this)
{
  if(comm->me == 0 && receive_from_world_ >= 0) {
    MPI_Recv(&ninsert_this, 1, MPI_INT, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
  }
  MPI_Bcast(&ninsert_this, 1, MPI_INT, 0, world);
}

/* ---------------------------------------------------------------------- */

double FixInsertPackFaceUniverse::insertion_fraction()
{
  // have to re-calculate region_volume_local in case simulation box is changing
  if (domain->box_change)
    calc_region_volume_local();

  for (std::map<int,int>::const_iterator it=faceid2index_.begin(); it!=faceid2index_.end(); ++it) {
    fraction_face_local[it->second] = insertion_fraction_face(it->first);
  }

  // mc method may be too vague -> normalize fraction
  int nfaces = faceid2index_.size();
  std::vector<double> fraction_face_all(nfaces, 1.0);
  MPI_Allreduce(&fraction_face_local[0], &fraction_face_all[0], nfaces, MPI_DOUBLE, MPI_SUM, world);

  for (int i=0; i<nfaces; ++i) {
    fraction_face_local[i] /= fraction_face_all[i];
  }

  return region_volume_local/region_volume;
}

/* ---------------------------------------------------------------------- */

int FixInsertPackFaceUniverse::get_face_index_check(int face_id)
{
  // check face_id
  if (faceid2index_.find(face_id) == faceid2index_.end())
    error->fix_error(FLERR,this,"invalid face id");

  return faceid2index_.at(face_id);
}

/* ----------------------------------------------------------------------
   distribute insertions across processors
------------------------------------------------------------------------- */

int FixInsertPackFaceUniverse::distribute_ninsert_this(int ninsert_this)
{
  int me, nprocs, ninsert_this_local=0;

  me = comm->me;
  nprocs = comm->nprocs;

  insertion_fraction();

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

  std::vector<DiscreteParticleDistribution> distributions_face(nfaces);

  // NOTE: fraction_face_local not sorted by face ID! --> fraction_face_local[faceid2index_[face_id]];
  //       distributions_face not sorted by face ID!
  std::vector<std::vector<std::vector<int> > > distributions_face_local_all(nprocs); //procs - faces - distributions:#particles
  for (int iproc = 0; iproc < nprocs; ++iproc) {
    distributions_face_local_all[iproc].resize(nfaces);
  }

  if(ninsert_this <= 0) {
    std::fill(nparticles_face_.begin(),nparticles_face_.end(), 0.);
    std::fill(average_vx_face_out_.begin(),average_vx_face_out_.end(), 0.);
    std::fill(average_vy_face_out_.begin(),average_vy_face_out_.end(), 0.);
    std::fill(average_vz_face_out_.begin(),average_vz_face_out_.end(), 0.);
    std::fill(average_omegax_face_out_.begin(),average_omegax_face_out_.end(), 0.);
    std::fill(average_omegay_face_out_.begin(),average_omegay_face_out_.end(), 0.);
    std::fill(average_omegaz_face_out_.begin(),average_omegaz_face_out_.end(), 0.);
  } else {
    // receive distribution data size
    int recv_data_size = 0;
    if(comm->me == 0 && receive_from_world_ >= 0) {
      MPI_Recv(&recv_data_size, 1, MPI_INT, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&recv_data_size, 1, MPI_INT, 0, world);

    // receive distribution data
    std::vector<double> recv_data(recv_data_size, 0.);
    if(comm->me == 0 && receive_from_world_ >= 0) {
      MPI_Recv(&recv_data[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&recv_data[0], recv_data_size, MPI_DOUBLE, 0, world);

    // build distributions for each face
    int ncpts = recv_data_size/(3+nfaces);
    for (int icpts=0; icpts<ncpts; ++icpts) {
      double rad = recv_data[icpts*(3+nfaces)];
      double mass = recv_data[icpts*(3+nfaces)+1];
      int type = static_cast<int>(recv_data[icpts*(3+nfaces)+2]);

      for (int iface=0; iface<nfaces; ++iface) {
        double nparticles = recv_data[icpts*(3+nfaces)+3+iface];
        if (nparticles > 0.) {
          distributions_face[iface][ConstantParticleTemplateSphere(rad, mass, type)] = nparticles;
        }
      }
    }

    // receive per face data
    recv_data_size = faceid2index_.size();
    if(comm->me == 0 && receive_from_world_ >= 0) {
      MPI_Recv(&nparticles_face_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_vx_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_vy_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_vz_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_omegax_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_omegay_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
      MPI_Recv(&average_omegaz_face_out_[0], recv_data_size, MPI_DOUBLE, universe->root_proc[receive_from_world_], idmassflowface_hash, universe->uworld, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(&nparticles_face_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_vx_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_vy_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_vz_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_omegax_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_omegay_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);
    MPI_Bcast(&average_omegaz_face_out_[0], recv_data_size, MPI_DOUBLE, 0, world);

    // loop over faces
    for (int iface=0; iface<nfaces; ++iface) {
      std::vector<double> accumulatedParticleVolume_facei(nprocs, 0.);

      // loop over templates
      DiscreteParticleDistribution::const_iterator it_dist = distributions_face[iface].begin();
      for (int idist=0; it_dist!=distributions_face[iface].end(); ++it_dist, ++idist) {
        std::vector<int> accumulatedTemplateParticles_facei(distributions_face[iface].size(), 0);

        // scale number of particles by cg3_ factor (massflowface gives numbers for fully resolved level)
        double nparticles = it_dist->second / cg3_;
        double diameter = 2. * it_dist->first.radius_ * cg_;
        double volume_single = (1./6.)*MY_PI*diameter*diameter*diameter;

        for (int iproc = 0; iproc < nprocs; ++iproc) {
          double procvolfracface = fraction_face_local_all[iproc*nfaces+iface]; // local volume fraction
          double procminfaceextent = min_face_extent_local_all[iproc*nfaces+iface];
          // round to integer number of particles
          // rounding_all tries to avoid that first proc is preferred for insertion
          int nparts = diameter < procminfaceextent? static_cast<int>(procvolfracface*nparticles + rounding_all[iproc]) : 0;
          while(nparts > 0 && nparts + accumulatedTemplateParticles_facei[idist] > nparticles) {
            --nparts; // just in case some weird rounding issues occured
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
  fix_pddf->set_distribution_local(distributions_face, distributions_face_local_all[me], cg_, type_offset);
  maxrad = std::max(maxrad, fix_pddf->max_rad());

  delete [] fraction_face_local_all;
  delete [] min_face_extent_local_all;

  rounding_all.push_front(rounding_all.back());
  rounding_all.pop_back();

  return ninsert_this_local;
}

/* ----------------------------------------------------------------------
   generate random positions within insertion volume
   perform overlap check via xnear
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertPackFaceUniverse::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
  ninserted_this_local = ninserted_spheres_this_local = 0;
  mass_inserted_this_local = 0.;

  double pos[3];
  ParticleToInsert *pti;

  int ntry = 0;
  double v_toInsert[3];
  vectorZeroize3D(v_toInsert);

  // always do overlap check
  {
    int iface = 0;
    FixParticledistributionDiscreteFace *fix_pddf = (FixParticledistributionDiscreteFace*)fix_distribution;
    for (std::map<int,int>::const_iterator it=faceid2index_.begin(); it!=faceid2index_.end(); ++it, ++iface) {
      int iHex = ins_region_mesh_hex->get_hex("face_id", it->first);
      if (iHex < 0 || iHex >= ins_region_mesh_hex->n_hex())
        error->fix_error(FLERR, this, "Failed to find hex cell by face_id!");

      int particle_this_face = 0;
      int particle_this_face_local = fix_pddf->pti_list_face_local[iface].size();
      int ninserted_this_face_local = 0;

      if (cg_ > 1.0) {
        double particle_this_face_real = nparticles_face_[iface] + remaining_ptis[it->first];
        double remainder = fmod(particle_this_face_real, cg3_);
        particle_this_face = static_cast<int>(particle_this_face_real / cg3_);
        if ((remainder/cg3_) > 0.5) {
          remainder -= cg3_;
          ++particle_this_face;
        }
        remaining_ptis[it->first] = remainder;
      } else {
        particle_this_face = nparticles_face_[iface];
      }

      v_insert[0] = average_vx_face_out_[iface];
      v_insert[1] = average_vy_face_out_[iface];
      v_insert[2] = average_vz_face_out_[iface];
      omega_insert[0] = average_omegax_face_out_[iface];
      omega_insert[1] = average_omegay_face_out_[iface];
      omega_insert[2] = average_omegaz_face_out_[iface];

      ntry = 0;

      int maxtry = calc_maxtry(particle_this_face_local);
      while (ntry < maxtry && ninserted_this_face_local < particle_this_face_local) {
        pti = fix_pddf->pti_list_face_local[iface][ninserted_this_face_local];
        double rbound = pti->r_bound_ins;

        int nins = 0;
        while (nins == 0 && ntry < maxattempt) {
          do {
            // generate a point in my subdomain
            do {
              do {
                ins_region_mesh_hex->hex_randpos(iHex, pos);
              } while (!domain->is_in_subdomain(pos));
            } while (all_in_flag && ins_region_mesh_hex->match_hex_cut(iHex,pos,rbound));
            ++ntry;
          } while (ntry < maxattempt && domain->dist_subbox_borders(pos) < rbound);

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

