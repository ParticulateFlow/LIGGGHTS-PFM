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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include "fix_move_mesh.h"
#include "fix_mesh.h"
#include "tri_mesh.h"
#include "modify.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "mesh_mover.h"
#include "container.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMoveMesh::FixMoveMesh(LAMMPS *lmp, int narg, char **arg)
  : Fix(lmp,narg,arg),
    fix_mesh_(0),
    mesh_(0),
    move_(0),
    time_(0),
    time_since_setup_(0)
{
    vectorZeroize3D(reference_point_);

    if(narg < 6)
      error->all(FLERR,"Illegal fix move/mesh command, you need to specify a mesh");

    int iarg = 3;

    if(strcmp(arg[iarg++],"mesh"))
      error->all(FLERR,"Illegal fix move/mesh command, expecting keyword 'mesh'");

    fix_mesh_ = static_cast<FixMesh*>(modify->find_fix_id(arg[iarg++]));
    if(fix_mesh_ == 0)
        error->all(FLERR,"Illegal fix move/mesh command, illegal mesh ID provided");

    mesh_ = fix_mesh_->mesh();
    move_ = createMeshMover(lmp,mesh_,this,&arg[iarg],narg-iarg);

    if(move_ == 0)
      error->all(FLERR,"Illegal fix move/mesh command, illegal arguments");

    // not compatible because surface velocity is just set once
    // and for moving mesh it is set every step
    if(fix_mesh_->surfaceVel())
      error->all(FLERR,"Illegal fix move/mesh command, cannot apply move to a mesh using keywords 'velocity' or 'angular_velocity'");

    //NP indicate that restart data is stored
    restart_global = 1;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::pre_delete(bool unfixflag)
{
    // check if another fix move operates on the same mesh
    // which came after me in the imput script
    // if so, it is illegal to delete this command
    // without first deleting the other

    //NP this is because node_orig_ would be messed up
    //NP this way, user is forced to unfix and re-instantiate
    //NP all commands, so node_orig_ is recreated

    //NP can throw errors only for case unfixflag = true
    //NP otherwise error would produce segfault

    if(unfixflag)
    {
        int nmove = modify->n_fixes_style("move/mesh");

        for(int imove = 0; imove < nmove; imove++)
        {
            FixMoveMesh* fix_move_mesh = static_cast<FixMoveMesh*>(modify->find_fix_style("move/mesh",imove));
            if(fix_move_mesh != this && fix_move_mesh->fixMesh() == fixMesh() && move_->isFirst())
                error->all(FLERR,"Illegal deletion of a fix move/mesh. There is another fix move/mesh command active on the same mesh. "
                           "Superposed fix move/mesh commands must be unfixed in reverse order of creation");
        }

        //NP have MeshMover unregister with mesh
        move_->pre_delete();

        // do not delete property v, as a dump command may still refer to it
        // set velocity to zero
        //NP if other fix move/mesh is still active, it will correctly set v on next step
        MultiVectorContainer<double,3,3> *v;
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        if(v) v->setAll(0.);

        // remove reference point if have one
        char refpt_id[200];
        sprintf(refpt_id, "REFPT_%s",id);

        if(mesh_->prop().getGlobalProperty<   VectorContainer<double,3> >(refpt_id))
           mesh_->prop().removeGlobalProperty(refpt_id);

    }

    delete move_;
}

/* ---------------------------------------------------------------------- */

FixMoveMesh::~FixMoveMesh()
{
    //NP do not need to delete v here, destructor is called in AssociativePtrArray
}

/* ---------------------------------------------------------------------- */

int FixMoveMesh::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::setup(int vflag)
{
    /*NL*///  fprintf(screen,"FixMoveMesh::setup()\n");
    time_since_setup_ = 0.;
    /*NL*/// fprintf(screen,"time_since_setup_ %f\n",time_since_setup_);

    reset_reference_point();

    if(!mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v"))
    {
        mesh_->prop().addElementProperty<MultiVectorContainer<double,3,3> >("v","comm_none","frame_invariant","restart_no");
        /*NL*/// int size = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v")->size();
        /*NL*/// fprintf(screen,"proc %d, FixMoveMesh::setup size %d\n",comm->me,size);
    }

    move_->setup();
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::initial_integrate(int dummy)
{
    MultiVectorContainer<double,3,3> *v;

    double dt = update->dt;
    time_ += update->dt;
    time_since_setup_ += update->dt;

    //NP reset velocity if I am first fix move/mesh on this mesh
    //NP first fix will always be first as it cant be deleted without
    //NP deleting all other fixes

    if(move_->isFirst())
    {
        v = mesh_->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");
        v->setAll(0.);
    }

    // integration
    //NP need to pass time since set up b/c node is copied to node_orig_ upon setup
    move_->initial_integrate(time_,time_since_setup_,dt);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::final_integrate()
{
    double dt = update->dt;

    // useful only if accelerations are known
    //NP need to pass time since set up b/c node is copied to node_orig_ upon setup
    move_->final_integrate(time_,time_since_setup_,dt);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixMoveMesh::write_restart(FILE *fp)
{
  int n = 0;
  double * list = new double[1 + move_->n_restart()];
  list[n++] = time_;
  /*NL*///fprintf(screen,"time %f\n",time_);

  move_->write_restart(&(list[n]));
  n += move_->n_restart();

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  delete[] list;
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixMoveMesh::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  time_ = static_cast<double> (list[n++]);
  move_->read_restart(&(list[n]));

  /*NL*///fprintf(screen,"time %f\n",time_);
}

/* ----------------------------------------------------------------------
   called by mesh mover
------------------------------------------------------------------------- */

void FixMoveMesh::add_reference_point(double *point)
{
    char refpt_id[200];
    sprintf(refpt_id, "REFPT_%s",id);

    if(mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id))
        error->fix_error(FLERR,this,"only one reference point allowed");

    /*NL*/ //printVec3D(screen,"adding point",point);
    vectorCopy3D(point,reference_point_);

    mesh_->prop().addGlobalProperty<VectorContainer<double,3> >(refpt_id,"comm_none","frame_general","restart_no");
    mesh_->prop().setGlobalProperty<VectorContainer<double,3> >(refpt_id,point);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::get_reference_point(double *point)
{
    VectorContainer<double,3> *refpt;
    char refpt_id[200];

    sprintf(refpt_id, "REFPT_%s",id);
    refpt = mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id);

    if(!refpt)
        error->fix_error(FLERR,this,"internal error");

    //NP need to explicitly reset reference point
    //NP otherwise would be too late since reset is called in rotate()
    //NP or move() only and first mesh mover needs resetted reference point
    //NP so only do this for first mesh mover
    if(move_->isFirst())
        mesh_->prop().resetGlobalPropToOrig(refpt_id);

    refpt->get(0,point);
    vectorCopy3D(point,reference_point_);
    /*NL*/ //printVec3D(screen,"getting point",point);
}

/* ---------------------------------------------------------------------- */

void FixMoveMesh::reset_reference_point()
{
    //NP need to re-set reference point from local copy upon setup
    //NP this ensures orig value of reference_point for fix i has been
    //NP handled by fixes i-1, i-2,... only (not by i, i+1,...)

    VectorContainer<double,3> *refpt;
    char refpt_id[200];

    sprintf(refpt_id, "REFPT_%s",id);
    refpt = mesh_->prop().getGlobalProperty<VectorContainer<double,3> >(refpt_id);

    // no error since not all moves have reference points
    if(!refpt)
        return;

    // set value for property
    refpt->set(0,reference_point_);

    // set orig value for property
    mesh_->prop().storeGlobalPropOrig(refpt_id);

    /*NL*/ //printVec3D(screen,"re-setting point",reference_point_);
}
