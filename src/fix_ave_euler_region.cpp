/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Department for Particule Flow Modelling
   Copyright 2016- JKU Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:
   Daniel Queteschiner <daniel.queteschiner@jku.at> (JKU Linz)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi_liggghts.h"
#include "fix_ave_euler_region.h"
#include "compute_stress_atom.h"
#include "math_extra_liggghts.h"
#include "atom.h"
#include "force.h"
#include "domain.h"
#include "modify.h"
#include "neighbor.h"
#include "region.h"
#include "region_mesh_hex.h"
#include "update.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

#define INVOKED_PERATOM 8

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAveEulerRegion::FixAveEulerRegion(LAMMPS *lmp, int narg, char **arg) :
  FixAveEuler(lmp, narg, arg),
  idregion_grid_(NULL),
  region_grid_(NULL),
  region_grid_mesh_hex_(NULL),
  v_min_(NULL),
  v_max_(NULL)
{
  // parse args
  if (narg < 4) error->all(FLERR,"Illegal fix ave/euler/region command");
  int iarg = 3;

  if(strcmp(arg[iarg++],"nevery"))
    error->fix_error(FLERR,this,"expecting keyword 'nevery'");
  exec_every_ = force->inumeric(FLERR,arg[iarg++]);
  if(exec_every_ < 1)
    error->fix_error(FLERR,this,"'nevery' > 0 required");
  nevery = exec_every_;

  parallel_ = false;

  while(iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion_grid_ = new char[n];
      strcpy(idregion_grid_,arg[iarg+1]);
      region_grid_ = domain->regions[iregion];
      if (strcmp(region_grid_->style, "mesh/hex") == 0) {
        region_grid_mesh_hex_ = static_cast<RegHexMesh*>(region_grid_);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"basevolume_region") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->fix_error(FLERR,this,"region ID does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion_ = new char[n];
      strcpy(idregion_,arg[iarg+1]);
      region_ = domain->regions[iregion];
      iarg += 2;
    } else if (strcmp(style,"fix ave/euler/region") == 0) {
      char *errmsg = new char[strlen(arg[iarg])+50];
      sprintf(errmsg,"unknown keyword or wrong keyword order: %s", arg[iarg]);
      error->fix_error(FLERR,this,errmsg);
      delete []errmsg;
    } else {
      ++iarg;
    }
  }
}

/* ---------------------------------------------------------------------- */

FixAveEulerRegion::~FixAveEulerRegion()
{
  delete [] idregion_grid_;
  memory->destroy(v_min_);
  memory->destroy(v_max_);
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegion::post_create()
{
  //  stress computation, just for pairwise contribution
  if (!compute_stress_) {
    const char* arg[5];
    arg[0]="stress_faveu";
    arg[1]="all";
    arg[2]="stress/atom";
    arg[3]="pair";
    arg[4]="fix";

    // create compute if it doesn't exist otherwise reuse
    if(modify->find_compute(arg[0]) < 0)
      modify->add_compute(5,(char**)arg);

    compute_stress_ = static_cast<ComputeStressAtom*>(modify->compute[modify->find_compute(arg[0])]);
  }


  if (region_grid_mesh_hex_) {
    ncells_ = region_grid_mesh_hex_->n_hex();
    ScalarContainer<int> *values = region_grid_mesh_hex_->prop().getElementProperty<ScalarContainer<int> >("cell_id");
    if (values) {
      for (int iHex=0; iHex<region_grid_mesh_hex_->n_hex(); ++iHex) {
        int cell_id = (*values)(iHex);
        if (cellid2index_.find(cell_id) != cellid2index_.end())
          error->warning(FLERR, "Cell ids in hexahedral mesh are not unique!");
        cellid2index_[cell_id] = iHex;
        cellid_.push_back(cell_id);
      }
    } else error->fix_error(FLERR, this, "requires cell data 'cell_id'");
  } else {
    error->fix_error(FLERR, this, "requires region of style mesh/hex");
  }

  send_post_create_data();
}

/* ---------------------------------------------------------------------- */

int FixAveEulerRegion::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegion::init()
{
  if (!atom->radius_flag)
    error->fix_error(FLERR,this,"requires atom attribute radius");
  if (!atom->rmass_flag)
    error->fix_error(FLERR,this,"requires atom attribute mass");

  if (region_) {
    int iregion = domain->find_region(idregion_);
    if (iregion == -1)
      error->fix_error(FLERR,this,"regions used by this command must not be deleted");

    region_ = domain->regions[iregion];
  }

  if (region_grid_) {
    int iregion = domain->find_region(idregion_grid_);
    if (iregion == -1)
      error->fix_error(FLERR,this,"regions used by this command must not be deleted");

    region_grid_ = domain->regions[iregion];

    if (strcmp(region_grid_->style, "mesh/hex") == 0) {
      region_grid_mesh_hex_ = static_cast<RegHexMesh*>(region_grid_);
    }
  }

  // error checks

  if (domain->triclinic == 1)
    error->fix_error(FLERR,this,"triclinic boxes not supported");
}

/* ----------------------------------------------------------------------
   setup 3d bins and their extent and coordinates
   called at setup() and when averaging occurs if box size changes
   similar to FixAveSpatial::setup_bins() and PairDSMC::init_style()

   bins are subbox - skin/2 so owned particles cannot move outside
   bins - so do not have to extrapolate
------------------------------------------------------------------------- */

void FixAveEulerRegion::setup_bins()
{
  ncells_ = region_grid_mesh_hex_->n_hex(); // global number of cells

  // (re) allocate spatial bin arrays
  if (ncells_ > ncells_max_) {
    ncells_max_ = ncells_;
    memory->grow(cellhead_,ncells_max_,"ave/euler:cellhead_");
    memory->grow(v_av_,  ncells_max_,3,"ave/euler:v_av_");
    memory->grow(vol_fr_,ncells_max_,  "ave/euler:vol_fr_");
    memory->grow(weight_,ncells_max_,  "ave/euler:vol_fr_");
    memory->grow(radius_,ncells_max_,  "ave/euler:radius_");
    memory->grow(ncount_,ncells_max_,  "ave/euler:ncount_");
    memory->grow(mass_,ncells_max_,    "ave/euler:mass_");
    memory->grow(stress_,ncells_max_,7,"ave/euler:stress_");
    memory->grow(v_min_,  ncells_max_,3,"ave/euler:v_min_");
    memory->grow(v_max_,  ncells_max_,3,"ave/euler:v_max_");
  }

  // calculate weight_[icell]
  if (!region_) {
    for (int icell = 0; icell < ncells_max_; ++icell)
      weight_[icell] = 1.;
  }

  // MC calculation if region_ exists
  if (region_) {
    double x_try[3];

    double contribution = 1./static_cast<double>(ntry_per_cell());  // contrib of each try

    for (int icell = 0; icell < ncells_max_; ++icell) {
      weight_[icell] = 0.;

      for (int itry=0; itry<ntry_per_cell(); ++itry) {
        region_grid_mesh_hex_->hex_randpos(icell, x_try);
        if (!domain->is_in_subdomain(x_try)) continue;
        if (region_->match(x_try[0],x_try[1],x_try[2])) {
          weight_[icell] += contribution;
        }
      }
    }

    // allreduce weights
    MPI_Sum_Vector(weight_,ncells_,world);

    // limit weight to 1
    for (int icell = 0; icell < ncells_max_; ++icell)
        if (weight_[icell] > 1.) weight_[icell] = 1.;
  }
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
   this also implies we do not need to wrap around PBCs
   bin ghost atoms only if inside my grid
------------------------------------------------------------------------- */

void FixAveEulerRegion::bin_atoms()
{
  int i,ibin;
  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < ncells_max_; ++i)
    cellhead_[i] = -1;

  // re-alloc cellptr_ if necessary
  if (nall > ncellptr_max_) {
    ncellptr_max_ = nall;
    memory->grow(cellptr_,ncellptr_max_,"ave/pic:cellptr_");
  }

  // bin in reverse order so linked list will be in forward order
  // also use ghost atoms
  // skip if any atom is out of (sub) box

  for (i = nall-1; i >= 0; --i) {
    if (!(mask[i] & groupbit)) continue;

    // skip particles outside my subdomain
    if (!domain->is_in_subdomain(x[i])) continue;

    ibin = region_grid_mesh_hex_->get_hex(x[i]);

    // particles outside grid may return values ibin < 0 || ibin >= ncells_
    // these are ignored

    if (ibin < 0 || ibin >= ncells_) continue;

    cellptr_[i] = cellhead_[ibin];
    cellhead_[ibin] = i;
  }
  dirty_ = false;
}

/* ---------------------------------------------------------------------- */

bool FixAveEulerRegion::is_inside_bounds(double bounds[6], double *pos)
{
  if(pos[0] < bounds[0] || pos[0] > bounds[1]) return false; // outside bb of hex
  if(pos[1] < bounds[2] || pos[1] > bounds[3]) return false; // outside bb of hex
  if(pos[2] < bounds[4] || pos[2] > bounds[5]) return false; // outside bb of hex
  return true;
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegion::lazy_bin_atoms(int icell)
{
  // skip if no particles in cell
  if (cellhead_[icell] == -1)
    return;

  double **x = atom->x;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;
  double bounds[6];

  region_grid_mesh_hex_->hex_bounds(icell,bounds);
  int previous_idx = -1;
  for (int iatom = cellhead_[icell]; iatom >= 0; ) {

    if (iatom > nall || // particle does no longer exist
        !(mask[iatom] & groupbit) || // particle does not belong to our group (e.g. due to deletions & copying inside atom array)
        !domain->is_in_subdomain(x[iatom]) || // particle no longer in subdomain
        !is_inside_bounds(bounds, x[iatom]) || // particle no longer inside AABB of this cell
        (region_grid_mesh_hex_->is_inside_hex(icell,x[iatom]) <= 0) ) { // particle no longer inside this cell

      if (iatom == cellhead_[icell]) {
        cellhead_[icell] = cellptr_[iatom];
        cellptr_[iatom] = -1;
        iatom = cellhead_[icell];
      } else if (previous_idx < 0) {
        cellptr_[cellhead_[icell]] = cellptr_[iatom];
        cellptr_[iatom] = -1;
        iatom = cellptr_[cellhead_[icell]];
      } else {
        cellptr_[previous_idx] = cellptr_[iatom];
        cellptr_[iatom] = -1;
        iatom = cellptr_[previous_idx];
      }
    } else {
      previous_idx = iatom;
      iatom = cellptr_[iatom];
    }
  }
}

/* ----------------------------------------------------------------------
   calculate Eulerian data, use interpolation function
------------------------------------------------------------------------- */

void FixAveEulerRegion::calculate_eu()
{
  double **v = atom->v;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double vel_x_mass[3] = {};

  // wrap compute with clear/add
  modify->clearstep_compute();

  // invoke compute if not previously invoked
  //NP also mess around with invoked_flag b/c addstep_compute
  //NP will only add a new step to those computes which were
  //NP invoked this time-step
  if (!(compute_stress_->invoked_flag & INVOKED_PERATOM)) {
    compute_stress_->compute_peratom();
    compute_stress_->invoked_flag |= INVOKED_PERATOM;
  }

  // forward comm per-particle stress from compute so neighs have it
  comm->forward_comm_compute(compute_stress_);

  // need to get pointer here since compute_peratom() may realloc
  double **stress_atom = compute_stress_->array_atom;

  // loop all binned particles
  // each particle can contribute to the cell that it has been binned to
  // optionally plus its 26 neighs

  for (int icell = 0; icell < ncells_; ++icell) {
    ncount_[icell] = 0;
    vectorZeroize3D(v_av_[icell]);
    vol_fr_[icell] = 0.;
    radius_[icell] = 0.;
    mass_[icell] = 0.;
    vectorZeroizeN(stress_[icell],7);
    vectorZeroize3D(v_min_[icell]);
    vectorZeroize3D(v_max_[icell]);

    // skip if no particles in cell

    if (-1 == cellhead_[icell])
      continue;

    // add contributions of particle - v and volume fraction
    // v is favre-averaged (mass-averaged)
    // radius is number-averaged

    for (int iatom = cellhead_[icell]; iatom >= 0; iatom = cellptr_[iatom]) {
      vectorScalarMult3D(v[iatom],rmass[iatom],vel_x_mass);
      vectorAdd3D(v_av_[icell],vel_x_mass,v_av_[icell]);
      vol_fr_[icell] += radius[iatom]*radius[iatom]*radius[iatom];
      radius_[icell] += radius[iatom];
      mass_[icell] += rmass[iatom];
      ncount_[icell]++;
      vectorComponentMin3D(v_min_[icell],v[iatom],v_min_[icell]);
      vectorComponentMax3D(v_max_[icell],v[iatom],v_max_[icell]);
    }
  }

  // allreduce contributions so far if not parallel
  if (ncells_ > 0) {
    MPI_Sum_Vector(&(v_av_[0][0]),3*ncells_,world);
    MPI_Sum_Vector(vol_fr_,ncells_,world);
    MPI_Sum_Vector(radius_,ncells_,world);
    MPI_Sum_Vector(mass_,ncells_,world);
    MPI_Sum_Vector(ncount_,ncells_,world);
    MPI_Min_Vector(&(v_min_[0][0]),3*ncells_,world);
    MPI_Max_Vector(&(v_max_[0][0]),3*ncells_,world);
  }

  // perform further calculations

  double eps_ntry = 1./static_cast<double>(ntry_per_cell());
  for (int icell = 0; icell < ncells_; ++icell) {

    // calculate average vel and radius
    if (ncount_[icell]) vectorScalarDiv3D(v_av_[icell],mass_[icell]);
    if (ncount_[icell]) radius_[icell]/=static_cast<double>(ncount_[icell]);

    // calculate volume fraction
    //safety check, add an epsilon to weight if any particle ended up in that cell
    if (vol_fr_[icell] > 0. && MathExtraLiggghts::compDouble(weight_[icell],0.,1e-6))
      weight_[icell] = eps_ntry;

    cell_volume_ = region_grid_mesh_hex_->hex_vol(icell);
    double prefactor_stress = 1./(cell_volume_*weight_[icell]);
    double prefactor_vol_fr = (4./3.)*M_PI*prefactor_stress;

    if (weight_[icell] < eps_ntry)
      vol_fr_[icell] = 0.;
    else
      vol_fr_[icell] *= prefactor_vol_fr;

    // add contribution of particle - stress
    // need v before can calculate stress
    // stress is molecular diffusion + contact forces

    for (int iatom = cellhead_[icell/*+stencil*/]; iatom >= 0; iatom = cellptr_[iatom]) {
      stress_[icell][1] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][0]-v_av_[icell][0]) + stress_atom[iatom][0];
      stress_[icell][2] += -rmass[iatom]*(v[iatom][1]-v_av_[icell][1])*(v[iatom][1]-v_av_[icell][1]) + stress_atom[iatom][1];
      stress_[icell][3] += -rmass[iatom]*(v[iatom][2]-v_av_[icell][2])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][2];
      stress_[icell][4] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][1]-v_av_[icell][1]) + stress_atom[iatom][3];
      stress_[icell][5] += -rmass[iatom]*(v[iatom][0]-v_av_[icell][0])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][4];
      stress_[icell][6] += -rmass[iatom]*(v[iatom][1]-v_av_[icell][1])*(v[iatom][2]-v_av_[icell][2]) + stress_atom[iatom][5];
    }

    stress_[icell][0] = -0.333333333333333*(stress_[icell][1]+stress_[icell][2]+stress_[icell][3]);
    if (weight_[icell] < eps_ntry)
      vectorZeroizeN(stress_[icell],7);
    else
      vectorScalarMultN(7,stress_[icell],prefactor_stress);
  }

  // allreduce stress if not parallel
  if (ncells_ > 0) {
    MPI_Sum_Vector(&(stress_[0][0]),7*ncells_,world);

    // recalc pressure based on allreduced stress
    for (int icell = 0; icell < ncells_; ++icell)
      stress_[icell][0] = -0.333333333333333*(stress_[icell][1]+stress_[icell][2]+stress_[icell][3]);
  }

  // wrap with clear/add
  int nextstep = (update->ntimestep/nevery)*nevery + nevery;
  modify->addstep_compute(nextstep);
  send_coupling_data();
}

/* ---------------------------------------------------------------------- */

double FixAveEulerRegion::cell_volume(int i)
{
  return region_grid_mesh_hex_->hex_vol(i);
}

/* ---------------------------------------------------------------------- */

double FixAveEulerRegion::cell_center(int i, int j)
{
  return region_grid_mesh_hex_->hex_center(i)[j];
}

/* ---------------------------------------------------------------------- */

void FixAveEulerRegion::cell_bounds(int i, double bounds[6])
{
  return region_grid_mesh_hex_->hex_bounds(i, bounds);
}

/* ---------------------------------------------------------------------- */

double* FixAveEulerRegion::cell_vector_property(int i, const char* property)
{
  return (*region_grid_mesh_hex_->prop().getElementProperty<VectorContainer<double,3> >(property))(i);
}

/* ----------------------------------------------------------------------
   return I,J array value
   if I exceeds current bins, return 0.0 instead of generating an error
   column 1,2,3 = bin coords, next column = vol fr,
   remaining columns = vel, stress, radius
------------------------------------------------------------------------- */

double FixAveEulerRegion::compute_array(int i, int j)
{
  if (i >= ncells_) return 0.0;
  else if(j < 3) return cell_center(i, j);
  else if(j == 3) return vol_fr_[i];
  else if(j < 7) return v_av_[i][j-4];
  else if(j == 7) return stress_[i][0];   // 0
  else if(j < 14) return stress_[i][j-7]; // 1 - 6
  else if(j < 15) return radius_[i];
  else if(j < 18) return v_min_[i][j-15];
  else if(j < 21) return v_max_[i][j-18];
  else if(j == 21) return ncount_[i];

  else return 0.0;
}

/* ---------------------------------------------------------------------- */

double FixAveEulerRegion::compute_array_by_id(int cell_id, int j)
{
  if (cellid2index_.find(cell_id) == cellid2index_.end())
      error->fix_error(FLERR, this, "Invalid cell id!");

  return compute_array(cellid2index_[cell_id], j);
}
