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

#include <algorithm>
#include <functional>
#include <stdlib.h>
#include <string.h>
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "container.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_massflow_mesh_face_universe.h"
#include "universe.h"

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace MathConst;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassflowMeshFaceUniverse::FixMassflowMeshFaceUniverse(LAMMPS *lmp, int narg, char **arg) :
  FixMassflowMeshFace(lmp, narg, arg)
{
    if (universe->nworlds == 1)
      error->all(FLERR,"Must have more than one processor partition for fix massflow/mesh/face/universe");

    cg_ = force->cg();
    cg3_ = cg_*cg_*cg_;

    id_hash_ = JSHash(id);
    couple_every_ = 0;
    send_to_world_ = -1;

    // parse args for this class
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if (strcmp(arg[iarg],"couple_every") == 0 || strcmp(arg[iarg],"every") == 0) {
          if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
          couple_every_ = atoi(arg[iarg+1]);
          if(couple_every_ < 0) error->fix_error(FLERR,this,"couple_every must be >= 0");
          iarg += 2;
          hasargs = true;
        } else if (strcmp(arg[iarg],"send_to_partition") == 0) {
          if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
          send_to_world_ = atoi(arg[iarg+1])-1;
          if(send_to_world_ < 0 || send_to_world_ >= universe->nworlds) error->fix_error(FLERR,this,"send_to_world must be a valid communicator");
          iarg += 2;
          hasargs = true;
        } else {
          ++iarg;
          hasargs = true;
        }
    }
}

/* ---------------------------------------------------------------------- */

FixMassflowMeshFaceUniverse::~FixMassflowMeshFaceUniverse()
{
}

/* ---------------------------------------------------------------------- */

template <typename K, typename V> V get_key(const std::pair<K,V>& p) { return p.first; }

template <typename K, typename V> void mapkeys2vector(const std::map<K, V> &m, std::vector<K>& vkey)
{
    vkey.clear();
    vkey.reserve(m.size());
    std::transform (m.begin(), m.end(), std::back_inserter(vkey), get_key<K,V>);
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFaceUniverse::send_post_create_data()
{
    // send information to child
    if (comm->me == 0 && send_to_world_ >= 0)
    {
        std::vector<int> keys;
        int nfaceids = faceid2index_.size();
        mapkeys2vector(faceid2index_, keys);
        MPI_Send(&nfaceids, 1, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
        MPI_Send(&keys[0], nfaceids, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
    }
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFaceUniverse::send_coupling_data()
{
    if((couple_every_ > 0) &&  (update->ntimestep > 0) && ((update->ntimestep-1) % couple_every_ == 0))
    {
        // only proc 0 sends data
        if(comm->me == 0)
        {
            int particles_since_last = (nparticles_-nparticles_last_)*cg3_;

            MPI_Send(&particles_since_last, 1, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);

            if(particles_since_last > 0)
            {
                int send_data_size = distributions_faces_.size()*(3+faceid2index_.size());
                std::vector<double> send_data;
                send_data.reserve(send_data_size);
                for(std::map<ConstantParticleTemplateSphere, std::vector<double> >::iterator it=distributions_faces_.begin(); it!=distributions_faces_.end(); ++it)
                {
                    send_data.push_back(it->first.radius_);
                    send_data.push_back(it->first.mass_);
                    send_data.push_back(it->first.atomtype_);
                    send_data.insert(send_data.end(), it->second.begin(), it->second.end());
                }

                MPI_Send(&send_data_size, 1, MPI_INT, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
                MPI_Send(&send_data[0], send_data_size, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);

                // send per face data
                send_data_size = faceid2index_.size();
                send_data.clear();
                send_data.resize(send_data_size);

                // send # particles per face
                std::transform (nparticles_face_.begin(), nparticles_face_.end(), nparticles_face_last_.begin(), send_data.begin(), std::minus<int>());
                std::transform (send_data.begin(), send_data.end(), send_data.begin(), std::bind1st(std::multiplies<double>(),cg3_));
                MPI_Send(&send_data[0], send_data_size, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);

                std::vector<double> mass_face_since_last_(faceid2index_.size(), 0.);
                std::transform (mass_face_.begin(), mass_face_.end(), mass_face_last_.begin(), mass_face_since_last_.begin(), std::minus<double>());

                // send vx per face
                for(int i=0; i<send_data_size; ++i)
                    send_data[i] = (mass_face_since_last_[i] > 0.) ? average_vx_face_out_[i]/mass_face_since_last_[i] : 0.;
                MPI_Send(&send_data[0], send_data_size, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
                // send vy per face
                for(int i=0; i<send_data_size; ++i)
                    send_data[i] = (mass_face_since_last_[i] > 0.) ? average_vy_face_out_[i]/mass_face_since_last_[i] : 0.;
                MPI_Send(&send_data[0], send_data_size, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
                // send vz per face
                for(int i=0; i<send_data_size; ++i)
                    send_data[i] = (mass_face_since_last_[i] > 0.) ? average_vz_face_out_[i]/mass_face_since_last_[i] : 0.;
                MPI_Send(&send_data[0], send_data_size, MPI_DOUBLE, universe->root_proc[send_to_world_], id_hash_, universe->uworld);
            }
        }

        if(reset_t_count_)
        {
            delta_t_ = t_count_;
            t_count_ = 0.;
            reset_t_count_ = false;
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFaceUniverse::reset_distributions(int /*size*/)
{
    distributions_faces_.clear();
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFaceUniverse::increment_distribution(const ConstantParticleTemplateSphere& cpts, int iface)
{
    distributions_faces_[cpts].resize(faceid2index_.size(), 0.);
    distributions_faces_[cpts][iface] += cg3_;
}

