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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include "fix_mesh_surface.h"
#include "fix_neighlist_mesh.h"
#include "fix_multisphere.h"
#include "fix_property_atom.h"
#include "mpi_liggghts.h"
#include "math_extra_liggghts.h"
#include "fix_massflow_mesh.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMassflowMesh::FixMassflowMesh(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  once_(false),
  delete_atoms_(false),
  mass_deleted_(0.),
  nparticles_deleted_(0),
  fix_orientation_(0),
  fix_mesh_(0),
  fix_counter_(0),
  fix_neighlist_(0),
  havePointAtOutlet_(false),
  insideOut_(false),
  mass_(0.),
  nparticles_(0),
  fix_property_(0),
  property_sum_(0.),
  verbose_(false),
  screenflag_(false),
  fp_(0),
  writeTime_(false),
  mass_last_(0.),
  nparticles_last_(0.),
  t_count_(0.),
  delta_t_(0.),
  reset_t_count_(true),
  fix_ms_(0),
  ms_(0),
  ms_counter_(0)
{
    vectorZeroize3D(nvec_);
    vectorZeroize3D(pref_);
    vectorZeroize3D(sidevec_);
    vectorZeroize3D(pointAtOutlet_);

    // parse args for this class

    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"vec_side") == 0) {
            if(narg < iarg+4)
                error->fix_error(FLERR,this,"not enough arguments for 'vec_side'");
            iarg++;
            sidevec_[0] = atof(arg[iarg++]);
            sidevec_[1] = atof(arg[iarg++]);
            sidevec_[2] = atof(arg[iarg++]);
            if(vectorMag3D(sidevec_) == 0.)
                error->fix_error(FLERR,this,"vec_side > 0 required");
            hasargs = true;
        } else if(strcmp(arg[iarg],"mesh") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'mesh'");
            iarg++;
            fix_mesh_ = static_cast<FixMeshSurface*>(modify->find_fix_id_style(arg[iarg++],"mesh/surface"));
            if(!fix_mesh_)
                error->fix_error(FLERR,this,"fix mesh ID does not exist");
            hasargs = true;
        } else if(strcmp(arg[iarg],"sum_property") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'sum_property'");
            iarg++;
            fix_property_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(arg[iarg++],"property/atom","scalar",0,0,style));
            hasargs = true;
        } else if(strcmp(arg[iarg],"count") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'count'");
            iarg++;
            if(strcmp(arg[iarg],"once") == 0)
                once_ = true;
            else if(strcmp(arg[iarg],"multiple") == 0)
                once_ = false;
            else
                error->fix_error(FLERR,this,"expecting 'once' or 'multiple' after 'count'");
            iarg++;
            hasargs = true;
        } else if( strcmp(arg[iarg],"writeTime") == 0) {
            writeTime_ = true;
            iarg++;
            hasargs = true;
        } else if(strcmp(arg[iarg],"point_at_outlet") == 0) {
            if(narg < iarg+4)
                error->fix_error(FLERR,this,"not enough arguments for 'point_at_outlet'");
            havePointAtOutlet_ = true;
            iarg++;
            pointAtOutlet_[0] = atof(arg[iarg++]);
            pointAtOutlet_[1] = atof(arg[iarg++]);
            pointAtOutlet_[2] = atof(arg[iarg++]);
            hasargs = true;
        } else if(strcmp(arg[iarg],"inside_out") == 0) {
            insideOut_ = true;
            iarg++;
            if(!havePointAtOutlet_)
                error->fix_error(FLERR,this,"the setting 'inside_out' has no meaning in case you do not use 'point_at_outlet'");
            hasargs = true;
        } else if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0)
          {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");

            char* filecurrent = new char[strlen(arg[iarg+1]) + 8];
            if (1 < comm->nprocs) //open a separate file for each processor
                sprintf(filecurrent,"%s.%d",arg[iarg+1],comm->me);
            else  //open one file for proc 0
                sprintf(filecurrent,"%s",arg[iarg+1]);

            if (strcmp(arg[iarg],"file") == 0)
                fp_ = fopen(filecurrent,"w");
            else
                fp_ = fopen(filecurrent,"a");
            delete [] filecurrent;
            if (fp_ == NULL) {
                char str[512];
                sprintf(str,"Cannot open file %s",arg[iarg+1]);
                error->fix_error(FLERR,this,str);
            }
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"verbose") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'verbose'");
            if(strcmp(arg[iarg+1],"yes") == 0)
                verbose_ = true;
            else if(strcmp(arg[iarg+1],"no") != 0)
                error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'verbose'");
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"screen") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'screen'");
            if (strcmp(arg[iarg+1],"yes") == 0) screenflag_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) screenflag_ = false;
            else error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'screen'");
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"delete_atoms") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"not enough arguments for 'delete_atoms'");
            if (strcmp(arg[iarg+1],"yes") == 0) delete_atoms_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) delete_atoms_ = false;
            else error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'delete_atoms'");
            iarg += 2;
            hasargs = true;
        } else if(strcmp(style,"massflow/mesh") == 0)
            error->fix_error(FLERR,this,"unknown keyword");
    }

    if(fp_ && comm->nprocs > 1 && comm->me == 0 && screen)
      fprintf(screen,"**FixMassflowMesh: > 1 process - "
                     " will write to multiple files\n");

    if(fp_)
    {
        //write header
        fprintf(fp_,"# ID");

        if(writeTime_)
          fprintf(fp_," time ");

        fprintf(fp_," diameter x y z u v w");

        fprintf(fp_,"  (ex ey ez, color)\n");

        fflush(fp_);
    }

    // error checks on necessary args

    if(!once_ && delete_atoms_)
           error->fix_error(FLERR,this,"using 'delete_atoms yes' requires 'count once'");
    if( !once_ && havePointAtOutlet_)
        error->fix_error(FLERR,this,"setting 'point_at_outlet' requires 'count once'");
    if( vectorMag3D(sidevec_)==0. && !havePointAtOutlet_)
        error->fix_error(FLERR,this,"expecting keyword 'vec_side'");
    if (!fix_mesh_)
        error->fix_error(FLERR,this,"expecting keyword 'mesh'");

    // get reference point on face
    // calculate normalvec

    fix_mesh_->triMesh()->node(0,0,pref_);
    fix_mesh_->triMesh()->surfaceNorm(0,nvec_);
    double dot = vectorDot3D(nvec_,sidevec_);

    /*NL*/ //if (screen) printVec3D(screen,"nvec_",nvec_);
    /*NL*/ //if (screen) printVec3D(screen,"sidevec_",sidevec_);

    if(fabs(dot) < 1e-6 && !havePointAtOutlet_ )
        error->fix_error(FLERR,this,"need to change 'vec_side', it is currently in or to close to the mesh plane");
    else if(dot < 0.)
        vectorFlip3D(nvec_);

    restart_global = 1;

    vector_flag = 1;
    size_vector = 6;
    if(fix_property_)
        size_vector = 7;
    global_freq = 1; //NP available always
}

/* ---------------------------------------------------------------------- */

FixMassflowMesh::~FixMassflowMesh()
{
   if(fp_) fclose(fp_);
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::post_create()
{
    // add per-particle count flag

    const char * fixarg[9];

    sprintf(fixid_,"massflow_%s",id);
    fixarg[0]=fixid_;
    fixarg[1]="all";
    fixarg[2]="property/atom";
    fixarg[3]=fixid_;
    fixarg[4]="scalar"; //NP 1 scalar per particle to be registered
    fixarg[5]="yes";    //NP restart yes
    fixarg[6]="no";    //NP communicate ghost yes
    fixarg[7]="no";    //NP communicate rev no
    fixarg[8]="-1.";     //NP take 0 as default
    modify->add_fix(9,const_cast<char**>(fixarg));

    fix_counter_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid_,"property/atom","scalar",0,0,style));

    // add neighbor list

    fix_neighlist_ = fix_mesh_->createOtherNeighList(igroup,id);

    // need to find multisphere here to be able to add property

    fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0));
    if(fix_ms_)
    {
        ms_ = &fix_ms_->data();
        ms_counter_ = ms_->prop().addElementProperty< ScalarContainer<int> >("counter_ms","comm_exchange_borders","frame_invariant", "restart_yes");
        ms_counter_->setDefaultValue(-1);

        if(delete_atoms_)
            error->fix_error(FLERR,this,"can not use 'delete_atoms' with fix multisphere");
        if(!once_)
            error->fix_error(FLERR,this,"must use 'count once' with fix multisphere");
    }
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::pre_delete(bool unfixflag)
{
    if (unfixflag) modify->delete_fix(fixid_);
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixMassflowMesh::init()
{
    if (atom->rmass_flag == 0)
        error->fix_error(FLERR,this,"requires atoms have mass");

    if(delete_atoms_ && 1 != atom->map_style)
        error->fix_error(FLERR,this,"requires an atom map of type 'array', via an 'atom_modify map array' command");

    if(!fix_ms_ && static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0)))
        error->fix_error(FLERR,this,"fix multisphere must come before fix massflow/mesh in input script");
}

/* ---------------------------------------------------------------------- */

void FixMassflowMesh::setup(int vflag)
{
    // check if face planar
    //NP do this here since may only do this after mesh setup
    if(!fix_mesh_->triMesh()->isPlanar() && !havePointAtOutlet_)
       error->fix_error(FLERR,this,"requires a planar face mass flow measurement or using 'point_at_outlet'");
}

/* ---------------------------------------------------------------------- */

int FixMassflowMesh::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    if(delete_atoms_) mask |= PRE_EXCHANGE;
    return mask;
}

/* ----------------------------------------------------------------------
   evaluate mass which went thru face
   all nearby particles
  -1 = default value for counter, set if particle is not
       a neighbor of the TriMesh
   0 = was not on nvec_ last step and is thus eligible
   1 = was on nvec_ side last step
   2 = do not re-count, was counted already
------------------------------------------------------------------------- */

void FixMassflowMesh::post_integrate()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    int *mask = atom->mask;
    double *counter = fix_counter_->vector_atom;
    double dot,delta[3]={};
    double mass_this = 0.;
    int nparticles_this = 0.;
    double property_this = 0.;
    double deltan;
    int *tag = atom->tag;

    class FixPropertyAtom* fix_color=static_cast<FixPropertyAtom*>(modify->find_fix_property("color","property/atom","scalar",0,0,style,false));
    bool fixColFound = false;
    if (fix_color) fixColFound=true;

    fix_orientation_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("ex",
        "property/atom","vector",0,0,style,false));

    TriMesh *mesh = fix_mesh_->triMesh();
    int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

    // update time for counter
    // also store values for last invokation
    t_count_ += update->dt;
    if(!reset_t_count_)
    {
        nparticles_last_ = nparticles_;
        mass_last_ = mass_;
        reset_t_count_ = true;
    }

    /*NL*///if(fix_property_ && screen)
    /*NL*///    fprintf(screen,"FOUNDE PROP, id %s style %s\n",fix_property_->id,fix_property_->style);

    // loop owned and ghost triangles
    // count only if owned particle

    for(int iTri = 0; iTri < nTriAll; iTri++)
    {
        //NP    TODO.... maybe some particles are double-counted if they are in the
        //NP    neighlist multiple times (for different triangles)?? (for case multiple)

        const std::vector<int> & neighborList = fix_neighlist_->get_contact_list(iTri);
        const int numneigh = neighborList.size();
        for(int iNeigh = 0; iNeigh < numneigh; iNeigh++)
        {
            const int iPart = neighborList[iNeigh];

            /*NL*/ //if(screen && 523 == atom->tag[iPart]) fprintf(screen,"step " BIGINT_FORMAT ": checking particle tag %d\n",update->ntimestep,atom->tag[iPart]);

            // skip ghost particles
            if(iPart >= nlocal)
                continue;

            // skip particles not in fix group
            if (!(mask[iPart] & groupbit))
                continue;

            const int ibody = fix_ms_ ? ( (fix_ms_->belongs_to(iPart) > -1) ? (ms_->map(fix_ms_->belongs_to(iPart))) : -1 ) : -1;

            // in case of once_ == true, ignore everything which has been already counted
            if((ibody > -1) ?  ((*ms_counter_)(ibody) == 2) : (compDouble(counter[iPart],2.)) ) continue;

            if(havePointAtOutlet_)
            {
                //get the vector from the particle center
                //to the next triangle
                deltan = fix_mesh_->triMesh()->resolveTriSphereContact(iPart,iTri,radius[iPart],x[iPart],delta);
                if(deltan < radius[iPart])
                {
                    vectorSubtract3D(x[iPart],pointAtOutlet_,nvec_); //vector pointing to the particle location
                    dot = vectorDot3D(delta,nvec_);
                }
                else //particle is not overlapping with mesh, so continue
                    continue;

                if(insideOut_) dot =-dot;
            }
            else
            {
                vectorSubtract3D(x[iPart],pref_,delta);
                dot = vectorDot3D(delta,nvec_);
            }

            // first-time, just set 0 or 1 depending on what side of the mesh
            if((ibody > -1) ? ((*ms_counter_)(ibody) == -1) : (compDouble(counter[iPart],-1.)) )
            {
                if(ibody > -1)
                    (*ms_counter_)(ibody) = (dot <= 0.) ? 0. : 1.;
                else
                    counter[iPart] = (dot <= 0.) ? 0. : 1.;
                /*NL*/ //if(screen && 523 == atom->tag[iPart]) fprintf(screen,"    1 counter set to %f\n",counter[iPart]);
                continue;
            }

            // particle is now on nvec_ side
            if(dot > 0.)
            {
                //particle was not on nvec_ side before
                if((ibody > -1) ? ((*ms_counter_)(ibody) == 0) : (compDouble(counter[iPart],0.)) ) // compDouble(counter[iPart],0.))
                {
                    //NP each particle holds whole MS mass
                    mass_this += rmass[iPart];
                    nparticles_this ++;
                    if(fix_property_)
                    {
                        /*NL*/ //if (screen) fprintf(screen,"adding %e\n",fix_property_->vector_atom[iPart]);
                        property_this += fix_property_->vector_atom[iPart];
                    }

                    if(delete_atoms_)
                    {
                        //reset counter to avoid problems with other fixes & mark to be deleted
                        counter[iPart] = -1.0;
                        atom_tags_delete_.push_back(atom->tag[iPart]);
                    }

                    if (screenflag_ && screen)
                        fprintf(screen," %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g \n ",
                                       tag[iPart],2.*radius[iPart]/force->cg(),
                                       x[iPart][0],x[iPart][1],x[iPart][2],
                                       v[iPart][0],v[iPart][1],v[iPart][2]);
                    if(fp_)
                    {
                        fprintf(fp_,"%d", tag[iPart]);

                        if(writeTime_)
                            fprintf(fp_,"  %4.4g ", update->dt*update->ntimestep);

                        fprintf(fp_," %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g",
                                   2.*radius[iPart]/force->cg(),
                                   x[iPart][0],x[iPart][1],x[iPart][2],
                                   v[iPart][0],v[iPart][1],v[iPart][2]);

                        if(fix_orientation_)
                        {
                            double **orientation = NULL;
                            orientation = fix_orientation_->array_atom;
                            fprintf(fp_,"    %4.4g %4.4g %4.4g ",
                                    orientation[iPart][0], orientation[iPart][1], orientation[iPart][2]);
                        }
                        if (fixColFound)
                            fprintf(fp_,"    %4.0g ", fix_color->vector_atom[iPart]);
                        fprintf(fp_,"\n");
                        fflush(fp_);
                    }
                }

                if(ibody > -1)
                    (*ms_counter_)(ibody) = once_ ? 2. : 1.;
                else if(!delete_atoms_) //only set if not marked for deletion
                    counter[iPart] = once_ ? 2. : 1.;
                /*NL*/ //if(screen && 523 == atom->tag[iPart]) fprintf(screen,"    2 counter set to %f\n",counter[iPart]);
            }
            else // dot <= 0
            {
                if(ibody > -1)
                    (*ms_counter_)(ibody) = 0;
                else
                    counter[iPart] = 0.;
                /*NL*/ //if(screen && 523 == atom->tag[iPart]) fprintf(screen,"    3 counter set to %f\n",counter[iPart]);
            }
        }
    }

    MPI_Sum_Scalar(mass_this,world);
    MPI_Sum_Scalar(nparticles_this,world);
    if(fix_property_) MPI_Sum_Scalar(property_this,world);

    /*NL*/ //if (screen) fprintf(screen,"property_this %e\n",property_this);

    mass_ += mass_this;
    nparticles_ += nparticles_this;
    property_sum_ += property_this;

}

/* ----------------------------------------------------------------------
   perform particle deletion of marked particles
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixMassflowMesh::pre_exchange()
{
    if (delete_atoms_)
    {
        double mass_deleted_this_ = 0.;
        int nparticles_deleted_this_ = 0.;
        int tag_max = atom->tag_max();

        // delete particles

        while (atom_tags_delete_.size() > 0)
        {
            if(atom_tags_delete_[0] <= tag_max)
            {
                int iPart = atom->map(atom_tags_delete_[0]);

                if(iPart >= 0)
                {
                    mass_deleted_this_ += atom->rmass[iPart];
                    nparticles_deleted_this_++;

                    atom->avec->copy(atom->nlocal-1,iPart,1);

                    // update atom map
                    // need to do this since atom map is needed to get local index for deletion
                    atom->map_one(atom->tag[atom->nlocal-1], iPart); // update map for copied particle

                    atom->nlocal--;
                }
                else if (verbose_)
                {
                    // particle may have been removed already by a different deleting command
                    // e.g. if the particle is in the neighbor list of the meshes of multiple massflow/mesh fixes
                    error->fix_warning(FLERR, this, "failed to find atom for deletion (possibly already deleted by another deleting command)");
                }
            }
            else if (verbose_)
            {
                // particle may have been removed already by a different deleting command
                // e.g. if the particle is in the neighbor list of the meshes of multiple massflow/mesh fixes
                error->fix_warning(FLERR, this, "failed to find atom for deletion (possibly already deleted by another deleting command)");
            }

            atom_tags_delete_.erase(atom_tags_delete_.begin());
        }

        atom_tags_delete_.clear();

        MPI_Sum_Scalar(mass_deleted_this_,world);
        MPI_Sum_Scalar(nparticles_deleted_this_,world);

        mass_deleted_ += mass_deleted_this_;
        nparticles_deleted_ += nparticles_deleted_this_;

        //NP tags and maps
        if(nparticles_deleted_this_)
        {
            bigint nblocal = atom->nlocal;
            MPI_Allreduce(&nblocal,&atom->natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

            if (atom->tag_enable) {
              if (atom->map_style) {
                atom->nghost = 0;
                atom->map_init();
                atom->map_set();
              }
            }
        }
    }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMassflowMesh::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = mass_;
  list[n++] = t_count_;
  list[n++] = mass_last_;
  list[n++] = static_cast<double>(nparticles_last_);
  list[n++] = mass_deleted_;
  list[n++] = static_cast<double>(nparticles_deleted_);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMassflowMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  mass_ = list[n++];
  t_count_ = list[n++];
  mass_last_ = list[n++];
  nparticles_last_ = static_cast<int>(list[n++]);
  nparticles_ = nparticles_last_;
  mass_deleted_ = list[n++];
  nparticles_deleted_ = static_cast<int>(list[n++]);
}

/* ----------------------------------------------------------------------
   output mass
------------------------------------------------------------------------- */

double FixMassflowMesh::compute_vector(int index)
{
    if(reset_t_count_)
    {
        delta_t_ = t_count_;
        t_count_ = 0.;
        reset_t_count_ = false;
    }

    if(index == 0)
        return mass_;
    if(index == 1)
        return static_cast<double>(nparticles_);
    if(index == 2)
        return delta_t_ == 0. ? 0. : (mass_-mass_last_)/delta_t_;
    if(index == 3)
        return delta_t_ == 0. ? 0. : static_cast<double>(nparticles_-nparticles_last_)/delta_t_;
    if(index == 4)
        return mass_deleted_;
    if(index == 5)
        return static_cast<double>(nparticles_deleted_);
    if(index == 6 && fix_property_)
        return property_sum_;

    return 0.;
}
