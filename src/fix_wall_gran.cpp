/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   This file was modified with respect to the release in LAMMPS
   Modifications are Copyright 2009-2012 JKU Linz
                     Copyright 2012-     DCS Computing GmbH, Linz

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_wall_gran.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair_gran.h"
#include "fix_rigid.h"
#include "fix_mesh.h"
#include "fix_contact_history_mesh.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_contact_property_atom_wall.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "compute_pair_gran_local.h"
#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface_stress.h"
#include "tri_mesh.h"
#include "primitive_wall.h"
#include "primitive_wall_definitions.h"
#include "mpi_liggghts.h"
#include "neighbor.h"
#include "contact_interface.h"
#include "fix_property_global.h"
#include <vector>
#include "granular_wall.h"
#include <assert.h>
#include <string>
#include <sstream>

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  #include "math_extra_liggghts_nonspherical.h"
#endif

using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::PRIMITIVE_WALL_DEFINITIONS;
using namespace LIGGGHTS::Walls;
using namespace LIGGGHTS::ContactModels;

/*NL*/ #define DEBUGMODE_LMP_FIX_WALL_GRAN false //(20950 < update->ntimestep) //(comm->me ==1)
/*NL*/ #define DEBUG_LMP_FIX_FIX_WALL_GRAN_M_ID 774
/*NL*/ #define DEBUG_LMP_FIX_FIX_WALL_GRAN_P_ID 206

const double SMALL = 1e-12;

  // modes for conduction contact area calaculation
  // same as in fix_heat_gran_conduction.cpp

  enum{ CONDUCTION_CONTACT_AREA_OVERLAP,
        CONDUCTION_CONTACT_AREA_CONSTANT,
        CONDUCTION_CONTACT_AREA_PROJECTION};

/* ---------------------------------------------------------------------- */

FixWallGran::FixWallGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    // wall/gran requires gran properties
    // sph not
    if (strncmp(style,"wall/gran",9) == 0 && (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag))
        error->fix_error(FLERR,this,"requires atom attributes radius, omega, torque");

    // defaults
    store_force_ = false;
    store_force_contact_ = false;
    stress_flag_ = false;
    n_FixMesh_ = 0;
    dnum_ = 0;
    skinDistance_ = 0.0;

    r0_ = 0.;

    shear_ = 0;
    shearDim_ = shearAxis_ = -1;
    vectorZeroize3D(shearAxisVec_);

    atom_type_wall_ = 1; // will be overwritten during execution, but other fixes require a value here

    // initializations
    fix_wallforce_ = 0;
    fix_wallforce_contact_ = 0;
    fix_rigid_ = NULL;
    heattransfer_flag_ = false;

    FixMesh_list_ = NULL;

    rebuildPrimitiveNeighlist_ = false;

    addflag_ = 0;
    cwl_ = NULL;

    computeflag_ = 1;

    meshwall_ = -1;

    Temp_wall = -1.;
    fixed_contact_area_ = 0.;
    Q = Q_add = 0.;

    area_calculation_mode_ = CONDUCTION_CONTACT_AREA_OVERLAP;

    // parse args
    //style = new char[strlen(arg[2])+2];
    //strcpy(style,arg[2]);

    iarg_ = 3;
    narg_ = narg;

    int nremaining = narg - 3;
    char ** remaining_args = &arg[3];

    int64_t variant = Factory::instance().selectVariant(style, nremaining, remaining_args);
    impl = Factory::instance().create(style, variant, lmp, this);

    iarg_ = narg - nremaining;

    bool hasargs = true;
    while(iarg_ < narg && hasargs)
    {
        hasargs = false;

        /*NL*/// if (screen) fprintf(screen,"arg %s\n",arg[iarg_]);

        if (strcmp(arg[iarg_],"primitive") == 0) {
           iarg_++;
           meshwall_ = 0;

           if (meshwall_ == 1)
             error->fix_error(FLERR,this,"'mesh' and 'primitive' are incompatible, choose either of them");

           if (strcmp(arg[iarg_++],"type"))
             error->fix_error(FLERR,this,"expecting keyword 'type'");
           atom_type_wall_ = force->inumeric(FLERR,arg[iarg_++]);
           if (atom_type_wall_ < 1 || atom_type_wall_ > atom->ntypes)
             error->fix_error(FLERR,this,"1 <= type <= max type as defined in create_box'");

           int n_primitives_ = 1;

           if (strcmp(arg[iarg_],"num") == 0) {
             ++iarg_;
             n_primitives_ = atoi(arg[iarg_++]);
           }

           for(int iWall = 0; iWall < n_primitives_; ++iWall) {

             char *wallstyle = arg[iarg_++];
             int nPrimitiveArgs = PRIMITIVE_WALL_DEFINITIONS::numArgsPrimitiveWall(wallstyle);
             /*NL*/// if (screen) fprintf(screen,"nPrimitiveArgs %d\n",nPrimitiveArgs);

             if(narg-iarg_ < nPrimitiveArgs)
              error->fix_error(FLERR,this,"not enough arguments for primitive wall");

             double * argVec = new double[nPrimitiveArgs];
             for(int i=0;i<nPrimitiveArgs;i++)
             {
               /*NL*/ //if (screen) fprintf(screen,"primitive arg %s\n",arg[iarg_]);
               argVec[i] = force->numeric(FLERR,arg[iarg_++]);
             }

             bool setflag = false;
             for(int w=0;w<(int)PRIMITIVE_WALL_DEFINITIONS::NUM_WTYPE;w++)
             {
               /*NL*///if (screen) fprintf(screen,"WallString %s %s\n",PRIMITIVE_WALL_DEFINITIONS::wallString[w],wallstyle);
               if(strcmp(wallstyle,PRIMITIVE_WALL_DEFINITIONS::wallString[w]) == 0)
               {
                 primitiveWalls_.push_back(new PrimitiveWall(lmp,(PRIMITIVE_WALL_DEFINITIONS::WallType)w,nPrimitiveArgs,argVec));
                 setflag = true;
                 break;
               }
             }
             if(!setflag) error->fix_error(FLERR,this,"unknown primitive wall style");
             delete[] argVec;
           }
           hasargs = true;
        } else if (strcmp(arg[iarg_],"mesh") == 0) {
           hasargs = true;
           meshwall_ = 1;
           iarg_ += 1;
        } else if (strcmp(arg[iarg_],"store_force") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"store_force_contact") == 0) {
           if (iarg_+2 > narg)
              error->fix_error(FLERR,this," not enough arguments");
           if (strcmp(arg[iarg_+1],"yes") == 0) store_force_contact_ = true;
           else if (strcmp(arg[iarg_+1],"no") == 0) store_force_contact_ = false;
           else error->fix_error(FLERR,this,"expecting 'yes' or 'no' after keyword 'store_force_contact_'");
           hasargs = true;
           iarg_ += 2;
        } else if (strcmp(arg[iarg_],"n_meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'n_meshes'");
          if (iarg_+2 > narg)
             error->fix_error(FLERR,this,"not enough arguments");
          n_FixMesh_ = atoi(arg[iarg_+1]);
          if(n_FixMesh_ < 1)
              error->fix_error(FLERR,this,"'n_meshes' > 0 required");
          hasargs = true;
          iarg_ += 2;
        } else if (strcmp(arg[iarg_],"meshes") == 0) {
          if (meshwall_ != 1)
             error->fix_error(FLERR,this,"have to use keyword 'mesh' before using 'meshes'");
          if(n_FixMesh_ == 0)
              error->fix_error(FLERR,this,"have to define 'n_meshes' before 'meshes'");
          if (narg < iarg_+1+n_FixMesh_)
              error->fix_error(FLERR,this,"not enough arguments");

          FixMesh_list_ = new FixMeshSurface*[n_FixMesh_];
          for(int i = 1; i <= n_FixMesh_; i++)
          {
              int f_i = modify->find_fix(arg[iarg_+i]);
              if (f_i == -1)
                  error->fix_error(FLERR,this,"could not find fix mesh id you provided");
              if (strncmp(modify->fix[f_i]->style,"mesh/surface",12))
                  error->fix_error(FLERR,this,"the fix belonging to the id you provided is not of type mesh");
              FixMesh_list_[i-1] = static_cast<FixMeshSurface*>(modify->fix[f_i]);

              if(FixMesh_list_[i-1]->trackStress())
                stress_flag_ = true;
              //NP create neighbor list for each mesh
          }
          hasargs = true;
          iarg_ += 1+n_FixMesh_;
        } else if (strcmp(arg[iarg_],"shear") == 0) {
          if (iarg_+3 > narg)
            error->fix_error(FLERR,this,"not enough arguments for 'shear'");
          if(primitiveWalls_.size() == 0)
            error->fix_error(FLERR,this,"have to define primitive wall before 'shear'. For mesh walls, please use fix move/mesh");

          if (strcmp(arg[iarg_+1],"x") == 0) shearDim_ = 0;
          else if (strcmp(arg[iarg_+1],"y") == 0) shearDim_ = 1;
          else if (strcmp(arg[iarg_+1],"z") == 0) shearDim_ = 2;
          else error->fix_error(FLERR,this,"illegal 'shear' dim");
          vshear_ = force->numeric(FLERR,arg[iarg_+2]);
          shear_ = 1;

          // update axis for cylinder etc if needed
          /* FIXME if(shearDim_ != primitiveWall_->axis())
          {
            shearAxis_ = primitiveWall_->axis();
            shearAxisVec_[shearAxis_] = vshear_;
          }*/

          hasargs = true;
          iarg_ += 3;
        } else if (strcmp(arg[iarg_],"temperature") == 0) {
            if (iarg_+1 >= narg)
              error->fix_error(FLERR,this,"not enough arguments for 'temperature'");
            Temp_wall = force->numeric(FLERR,arg[iarg_+1]);
            hasargs = true;
            iarg_ += 2;
        } else if(strcmp(arg[iarg_],"contact_area") == 0) {

          if(strcmp(arg[iarg_+1],"overlap") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_OVERLAP;
          else if(strcmp(arg[iarg_+1],"projection") == 0)
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_PROJECTION;
          else if(strcmp(arg[iarg_+1],"constant") == 0)
          {
            if (iarg_+3 > narg)
                error->fix_error(FLERR,this,"not enough arguments for keyword 'contact_area constant'");
            area_calculation_mode_ =  CONDUCTION_CONTACT_AREA_CONSTANT;
            fixed_contact_area_ = force->numeric(FLERR,arg[iarg_+2]);
            if (fixed_contact_area_ <= 0.)
                error->fix_error(FLERR,this,"'contact_area constant' value must be > 0");
            iarg_++;
          }
          else error->fix_error(FLERR,this,"expecting 'overlap', 'projection' or 'constant' after 'contact_area'");
          iarg_ += 2;
          hasargs = true;
        }
    }

    if(impl)
      impl->settings(narg - iarg_, &arg[iarg_]);

    // error checks

    if(meshwall_ == -1 && primitiveWalls_.size() == 0)
        error->fix_error(FLERR,this,"Need to use define style 'mesh' or 'primitive'");

    if(meshwall_ == 1 && !FixMesh_list_)
        error->fix_error(FLERR,this,"Need to provide the number and a list of meshes by using 'n_meshes' and 'meshes'");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_create()
{
    if(strncmp(style,"wall/gran",9) != 0)
    {
      // case non-granular (sph)
      dnum_ = 0;
    }

    // register storage for wall force if required
    if(store_force_)
    {
          char *wallforce_name = new char[strlen(id)+1+6];
          strcpy(wallforce_name,"force_");
          strcat(wallforce_name,id);
          const char *fixarg[11];
          fixarg[0] = wallforce_name;
          fixarg[1] = "all";
          fixarg[2] = "property/atom";
          fixarg[3] = wallforce_name;
          fixarg[4] = "vector";
          fixarg[5] = "no";    // restart
          fixarg[6] = "no";    // communicate ghost
          fixarg[7] = "no";    // communicate rev
          fixarg[8] = "0.";
          fixarg[9] = "0.";
          fixarg[10] = "0.";
          modify->add_fix(11,const_cast<char**>(fixarg));
          fix_wallforce_ =
              static_cast<FixPropertyAtom*>(modify->find_fix_property(wallforce_name,"property/atom","vector",3,0,style));
          delete []wallforce_name;
   }

   //NP case primitve wall - register here
   if(store_force_contact_ && 0 == meshwall_)
   {
        const char *fixarg[19];
        char fixid[200],ownid[200];
        sprintf(fixid,"contactforces_%s",id);
        sprintf(ownid,"%s",id);
        fixarg[0]=fixid;
        fixarg[1]="all";
        fixarg[2]="contactproperty/atom/wall";
        fixarg[3]=fixid;
        fixarg[4]="6";
        fixarg[5]="fx";
        fixarg[6]="0";
        fixarg[7]="fy";
        fixarg[8]="0";
        fixarg[9]="fz";
        fixarg[10]="0";
        fixarg[11]="tx";
        fixarg[12]="0";
        fixarg[13]="ty";
        fixarg[14]="0";
        fixarg[15]="tz";
        fixarg[16]="0";
        fixarg[17]="primitive";
        fixarg[18]=ownid;
        modify->add_fix(19,const_cast<char**>(fixarg));
        fix_wallforce_contact_ = static_cast<FixContactPropertyAtomWall*>(modify->find_fix_id(fixid));
   }

   // create neighbor list for each mesh
   //NP also create contact history in case of dnum == 0
   //NP so to be able to exclude duplicates contacts for coplanar faces
   for(int i=0;i<n_FixMesh_;i++)
   {
       //NP order is important!!
       //NP neigh list must be built before contact history is refreshed
       FixMesh_list_[i]->createWallNeighList(igroup);
       FixMesh_list_[i]->createContactHistory(dnum());

       if(store_force_contact_)
         FixMesh_list_[i]->createMeshforceContact();
   }

   // contact history for primitive wall
   if(meshwall_ == 0 && dnum_ > 0)
   {
     assert(primitiveWallsHistory_.size() == 0);

     for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        ostringstream os;
        os << "history_" << id << "_" << i;
        string hist_name = os.str();
        const char **fixarg = new const char*[8+dnum_];
        fixarg[0] = hist_name.c_str();
        fixarg[1] = "all";
        fixarg[2] = "property/atom";
        fixarg[3] = hist_name.c_str();
        fixarg[4] = "vector";
        fixarg[5] = "yes"; // restart
        fixarg[6] = "no";  // communicate ghost
        fixarg[7] = "no";  // communicate rev
        for(int j = 8; j < 8+dnum_; ++j)
          fixarg[j] = "0.";
        modify->add_fix(8+dnum_,const_cast<char**>(fixarg));
        FixPropertyAtom* history = static_cast<FixPropertyAtom*>(modify->find_fix_property(hist_name.c_str(),"property/atom","vector",dnum_,0,style));
        primitiveWallsHistory_.push_back(history);
        delete []fixarg;
     }
   }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::pre_delete(bool unfixflag)
{
    if(unfixflag && store_force_)
        modify->delete_fix(fix_wallforce_->id);
    if(unfixflag && primitiveWallsHistory_.size() > 0) {
      for(size_t i = 0; i < primitiveWallsHistory_.size(); ++i){
        modify->delete_fix(primitiveWallsHistory_[i]->id);
        primitiveWallsHistory_[i] = NULL;
      }
      primitiveWallsHistory_.clear();
    }

   if(unfixflag && store_force_contact_) {
     modify->delete_fix(fix_wallforce_contact_->id);
   }

    if(unfixflag)
    {
       for(int i=0;i<n_FixMesh_;i++)
       {
           //NP remove contact history as new fix wall gran might have a
           //NP different contact history needs
           FixMesh_list_[i]->deleteWallNeighList();
           FixMesh_list_[i]->deleteContactHistory();
       }
    }
}

/* ---------------------------------------------------------------------- */

FixWallGran::~FixWallGran()
{
    if(primitiveWalls_.size() > 0) {
      for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        delete primitiveWalls_[i];
        primitiveWalls_[i] = NULL;
      }
      primitiveWalls_.clear();
    }
    if(FixMesh_list_) delete []FixMesh_list_;
    delete impl;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::setmask()
{
    int mask = 0;
    mask |= PRE_NEIGHBOR;
    mask |= PRE_FORCE;
    mask |= POST_FORCE;
    mask |= POST_FORCE_RESPA;
    return mask;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::min_type() const
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::max_type() const
{
    return atom_type_wall_;
}

/* ---------------------------------------------------------------------- */

PrimitiveWall* FixWallGran::primitiveWall()
{ return primitiveWalls_[0]; } // FIXME

/* ---------------------------------------------------------------------- */

void FixWallGran::init()
{
    dt_ = update->dt;

    // case granular
    if(strncmp(style,"wall/gran",9) == 0)
    {
        // check if a fix rigid is registered - important for damp
        fix_rigid_ = static_cast<FixRigid*>(modify->find_fix_style_strict("rigid",0));

        if (strcmp(update->integrate_style,"respa") == 0)
          nlevels_respa_ = ((Respa *) update->integrate)->nlevels;

        if(impl)
          impl->init_granular();
        else
        {
          // init for derived classes
          init_granular();
        }

        // disallow more than one wall of non-rimitive style
        //NP because contact history is stored in the mesh
        //NP and mesh can hold only one history
        if(is_mesh_wall())
        {
            int nfix = modify->n_fixes_style("wall/gran");
            for (int ifix = 0; ifix < nfix; ifix++)
            {
                FixWallGran *fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
                if (fwg == this) continue;
                if (fwg->is_mesh_wall())
                    error->fix_error(FLERR,this,"More than one wall of type 'mesh' is not supported");
            }
        }
    }
    //NP out = fopen("shearhistory","w");
}

/* ---------------------------------------------------------------------- */

void FixWallGran::setup(int vflag)
{
    if (strstr(update->integrate_style,"verlet"))
    {
      pre_neighbor();
      pre_force(vflag);
      post_force(vflag);
    }
    else
    {
      ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa_-1);
      post_force_respa(vflag,nlevels_respa_-1,0);
      ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa_-1);
    }

    //NP doing this here because deltan_ratio is set in init() of fix heat/gran
    init_heattransfer();
}

/* ----------------------------------------------------------------------
   neighbor list via fix wall/gran is only relevant for primitive walls
------------------------------------------------------------------------- */

void FixWallGran::pre_neighbor()
{
    rebuildPrimitiveNeighlist_ = (primitiveWalls_.size() > 0);
}

void FixWallGran::pre_force(int vflag)
{
    const double halfskin = neighbor->skin*0.5;
    const int nlocal = atom->nlocal;

    x_ = atom->x;
    radius_ = atom->radius;
    cutneighmax_ = neighbor->cutneighmax;
#ifdef SUPERQUADRIC_ACTIVE_FLAG
    if(atom->superquadric_flag) {
      quat_ = atom->quaternion;
      shape_ = atom->shape;
      blockiness_ = atom->blockiness;
    }
#endif

    /*NL*/// if (screen) fprintf(screen,"cutneighmax_ %f, radius is %s, rebuild primitive %s\n",cutneighmax_,radius_?"notnull":"null",rebuildPrimitiveNeighlist_?"yes":"no");
    /*NL*/// error->all(FLERR,"end");

    // build neighlist for primitive walls
    //NP as in neighlist/mesh - hack for sph
    if(rebuildPrimitiveNeighlist_) {
      for(size_t i = 0; i < primitiveWalls_.size(); ++i) {
        primitiveWalls_[i]->buildNeighList(radius_ ? halfskin:(r0_+halfskin),x_,radius_,nlocal);
      }
    }

    rebuildPrimitiveNeighlist_ = false;
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via verlet
------------------------------------------------------------------------- */

void FixWallGran::post_force(int vflag)
{
    computeflag_ = 1;
    shearupdate_ = 1;
    if (update->setupflag) shearupdate_ = 0;
    addflag_ = 0;

    post_force_wall(vflag);
}

/* ----------------------------------------------------------------------
   force on each atom calculated via post_force
   called via compute wall/gran
------------------------------------------------------------------------- */

void FixWallGran::post_force_pgl()
{
    computeflag_ = 0;
    shearupdate_ = 0;
    addflag_ = 1;

    post_force_wall(0);
}

/* ----------------------------------------------------------------------
   post_force
------------------------------------------------------------------------- */

void FixWallGran::post_force_wall(int vflag)
{
  // set pointers and values appropriately
  nlocal_ = atom->nlocal;
  x_ = atom->x;
  f_ = atom->f;
  radius_ = atom->radius;
  rmass_ = atom->rmass;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
  if(atom->superquadric_flag) {
    quat_ = atom->quaternion;
    shape_ = atom->shape;
    blockiness_ = atom->blockiness;
  }
#endif

  if(fix_rigid_)
  {
      body_ = fix_rigid_->body;
      masstotal_ = fix_rigid_->masstotal;
  }

  if(fix_wallforce_)
    wallforce_ = fix_wallforce_->array_atom;

  cutneighmax_ = neighbor->cutneighmax;

  if(nlocal_ && !radius_ && r0_ == 0.)
    error->fix_error(FLERR,this,"need either per-atom radius or r0_ being set");

  if(store_force_)
  {
      for(int i = 0; i < nlocal_; i++)
      {
          vectorZeroize3D(wallforce_[i]);
      }
  }

  if(meshwall_ == 1)
    post_force_mesh(vflag);
  else
    post_force_primitive(vflag);

  if(meshwall_ == 0 && store_force_contact_)
    fix_wallforce_contact_->do_forward_comm();

  if(meshwall_ == 1 && store_force_contact_)
  {
    for(int imesh = 0; imesh < n_FixMesh_; imesh++)
        FixMesh_list_[imesh]->meshforceContact()->do_forward_comm();
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGran::post_force_respa(int vflag, int ilevel, int iloop)
{
    if (ilevel == nlevels_respa_-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   post_force for mesh wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_mesh(int vflag)
{
    //NP groupbit accounted for via neighlist/mesh

    // contact properties
    double v_wall[3],bary[3];
    double delta[3],deltan;
    const int nlocal = atom->nlocal;
    const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;

    CollisionData cdata;
    cdata.is_wall = true;
    cdata.computeflag = computeflag_;
    cdata.shearupdate = shearupdate_;

    /*NL*/// if(comm->me == 3 && update->ntimestep == 3735 && screen)
    /*NL*///   fprintf(screen,"proc 3 start\n");

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();
      // mark all contacts for deletion at this point
      //NP all detected contacts will be un-marked by fix_contact->handleContact()

      if(fix_contact) {
        FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
        fix_contact->resetDeletionPage(0);

        std::vector<int> & particle_indices = meshNeighlist->particle_indices;

        // mark all contacts for deletion at this point
        //NP all detected contacts will be un-marked by fix_contact->handleContact()
        for(std::vector<int>::iterator it = particle_indices.begin(); it != particle_indices.end(); ++it) {
          fix_contact->markForDeletion(0, *it);
        }
      }

      //NP extremely dirty; this is that FixWallGranBase can use fix_wallforce_contact_
      //NP pointer for both primitive and mesh case
      if(store_force_contact_)
        fix_wallforce_contact_ = FixMesh_list_[iMesh]->meshforceContact();
    }
    //if (screen) fprintf(screen, "[%d] Mark All took %g seconds.\n", tid, MPI_Wtime() - markAllTime);

    for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
    {
      TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
      const int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();
      FixContactHistoryMesh * const fix_contact = FixMesh_list_[iMesh]->contactHistory();

      // get neighborList and numNeigh
      FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
      std::vector<int>::iterator b = meshNeighlist->particle_indices.begin();
      std::vector<int>::iterator e = meshNeighlist->particle_indices.end();
      if(b == e) continue; // nothing to do (no neighbors)

      vectorZeroize3D(v_wall);
      MultiVectorContainer<double,3,3> *vMeshC = mesh->prop().getElementProperty<MultiVectorContainer<double,3,3> >("v");

      cdata.jtype = FixMesh_list_[iMesh]->atomTypeWall();

      // moving mesh
      if(vMeshC)
      {
        double ***vMesh = vMeshC->begin();

        // loop owned and ghost triangles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(atom->superquadric_flag)
            {
              Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);

              if(mesh->sphereTriangleIntersection(iTri, radius_[iPart], x_[iPart])) //check for Bounding Sphere-triangle intersection
              {
                deltan = mesh->resolveTriSuperquadricContact(iTri, delta, cdata.contact_point, particle, bary);
              }
              else
              {
                deltan = LARGE_TRIMESH;
              }

              cdata.is_non_spherical = true; //by default it is false
            }
            else
            {
              deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta,bary);
            }
#else
            deltan = mesh->resolveTriSphereContactBary(iPart,iTri,radius_ ? radius_[iPart]:r0_ ,x_[iPart],delta,bary);
#endif

            if(deltan > cutneighmax_) continue;

            bool intersectflag = (deltan <= 0);

            if(intersectflag || (radius_ && deltan < contactDistanceMultiplier*radius_[iPart]))
            {
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,cdata.contact_history)) continue;

              for(int i = 0; i < 3; i++)
                v_wall[i] = (bary[0]*vMesh[iTri][0][i] + bary[1]*vMesh[iTri][1][i] + bary[2]*vMesh[iTri][2][i]);

              cdata.i = iPart;
              cdata.deltan = -deltan;
              cdata.delta[0] = -delta[0];
              cdata.delta[1] = -delta[1];
              cdata.delta[2] = -delta[2];
              post_force_eval_contact(cdata, intersectflag, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }

          }
        }
      }
      // non-moving mesh - do not calculate v_wall, use standard distance function
      else
      {
        // loop owned and ghost particles
        for(int iTri = 0; iTri < nTriAll; iTri++)
        {
          const std::vector<int> & neighborList = meshNeighlist->get_contact_list(iTri);
          const int numneigh = neighborList.size();

          for(int iCont = 0; iCont < numneigh; iCont++)
          {
            const int iPart = neighborList[iCont];

            // do not need to handle ghost particles
            if(iPart >= nlocal) continue;

            int idTri = mesh->id(iTri);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
            if(atom->superquadric_flag)
            {
              Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);

              if(mesh->sphereTriangleIntersection(iTri, radius_[iPart], x_[iPart])) //check for Bounding Sphere-triangle intersection
              {
                deltan = mesh->resolveTriSuperquadricContact(iTri, delta, cdata.contact_point, particle);
              }
              else
              {
                deltan = LARGE_TRIMESH;
              }

              cdata.is_non_spherical = true; //by default it is false
            }
            else
            {
              deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
            }
#else
            deltan = mesh->resolveTriSphereContact(iPart,iTri,radius_ ? radius_[iPart]:r0_,x_[iPart],delta);
#endif

            if(deltan > cutneighmax_) continue;

            bool intersectflag = (deltan <= 0);

            //NP hack for SPH
            if(intersectflag || (radius_ && deltan < contactDistanceMultiplier*radius_[iPart]))
            {
              //NP continue in case already have a contact with a coplanar face
              if(fix_contact && ! fix_contact->handleContact(iPart,idTri,cdata.contact_history)) continue;

              cdata.i = iPart;
              cdata.deltan = -deltan;
              cdata.delta[0] = -delta[0];
              cdata.delta[1] = -delta[1];
              cdata.delta[2] = -delta[2];
              post_force_eval_contact(cdata, intersectflag, v_wall,iMesh,FixMesh_list_[iMesh],mesh,iTri);
            }
          }
        }
      }
    }

  // clean-up contacts
  //NP i.e. delete all which have not been detected this time-step
  //NP have to do this here: values cannot be deleted before
  //NP since must be available to copy contact history for coplanar neighs
  //double cleanupTime = MPI_Wtime();
  for(int iMesh = 0; iMesh < n_FixMesh_; iMesh++)
  {
    FixContactHistoryMesh *fix_contact = FixMesh_list_[iMesh]->contactHistory();
      // clean-up contacts
      //NP i.e. delete all which have not been detected this time-step
      if(fix_contact) {
        FixNeighlistMesh * meshNeighlist = static_cast<FixNeighlistMesh*>(FixMesh_list_[iMesh]->meshNeighlist());
        std::vector<int> & particle_indices = meshNeighlist->particle_indices;

        for(std::vector<int>::iterator it = particle_indices.begin(); it != particle_indices.end(); ++it) {
          fix_contact->cleanUpContact(*it);
        }
      }
  }
}

/* ----------------------------------------------------------------------
   post_force for primitive wall
------------------------------------------------------------------------- */

void FixWallGran::post_force_primitive(int vflag)
{
  int *mask = atom->mask;

  CollisionData cdata;
  cdata.is_wall = true;
  cdata.computeflag = computeflag_;
  cdata.shearupdate = shearupdate_;
  cdata.jtype = atom_type_wall_;
  const double contactDistanceMultiplier = neighbor->contactDistanceFactor - 1.0;

  // contact properties
  double v_wall[] = {0.,0.,0.};

  // if shear, set velocity accordingly
  if (shear_) v_wall[shearDim_] = vshear_; // FIXME

  for(size_t iWall = 0; iWall < primitiveWalls_.size(); ++iWall) {
    PrimitiveWall * primitiveWall = primitiveWalls_[iWall];
    double **c_history = 0;

    if(dnum() > 0) {
      c_history = primitiveWallsHistory_[iWall]->array_atom;
    }

    // loop neighbor list
    int *neighborList;
    const int nNeigh = primitiveWall->getNeighbors(neighborList);

    for (int iCont = 0; iCont < nNeigh ; iCont++, neighborList++)
    {
      int iPart = *neighborList;

      if(!(mask[iPart] & groupbit)) continue;

      double delta[3]={};
      double deltan = primitiveWall->resolveContact(x_[iPart],radius_?radius_[iPart]:r0_,delta);

      if(deltan > cutneighmax_) continue;

      if(deltan <= 0 || deltan < contactDistanceMultiplier*radius_[iPart])
      {
        // spheres
        bool intersectflag = (deltan <= 0);

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(atom->superquadric_flag) {
          double sphere_contact_point[3];
          vectorAdd3D(x_[iPart], delta, sphere_contact_point);
          double closestPoint[3], closestPointProjection[3], point_of_lowest_potential[3];

          Superquadric particle(x_[iPart], quat_[iPart], shape_[iPart], blockiness_[iPart]);
          intersectflag = particle.plane_intersection(delta, sphere_contact_point, closestPoint, point_of_lowest_potential);
          deltan = -MathExtraLiggghtsNonspherical::point_wall_projection(delta, sphere_contact_point, closestPoint, closestPointProjection);

          vectorCopy3D(closestPoint, cdata.contact_point);
          cdata.is_non_spherical = true; // by default it is false
        }
#endif

        if(shear_ && shearAxis_ >= 0)
        {
          double rdist[3];
          primitiveWall->calcRadialDistance(x_[iPart],rdist);
          vectorCross3D(shearAxisVec_,rdist,v_wall);
        }
        cdata.i = iPart;
        cdata.contact_history = c_history ? c_history[iPart] : NULL;

#ifdef SUPERQUADRIC_ACTIVE_FLAG
        if(atom->superquadric_flag && deltan > 0.0)
        {
          if(c_history)
            vectorZeroizeN(c_history[iPart],dnum_);
          break;
        }

        if(!cdata.is_non_spherical || atom->superquadric_flag)
#endif
        cdata.deltan = -deltan;

        cdata.delta[0] = -delta[0];
        cdata.delta[1] = -delta[1];
        cdata.delta[2] = -delta[2];
        post_force_eval_contact(cdata,intersectflag,v_wall);
      }
      else
      {
        if(c_history) vectorZeroizeN(c_history[iPart],dnum_);
      }
    }
  }
}

void FixWallGran::compute_force(CollisionData &, double *)
{
}

/* ----------------------------------------------------------------------
   actually calculate force, called for both mesh and primitive
------------------------------------------------------------------------- */

inline void FixWallGran::post_force_eval_contact(CollisionData & cdata, bool intersectflag, double * v_wall, int iMesh, FixMeshSurface *fix_mesh, TriMesh *mesh, int iTri)
{
  const int iPart = cdata.i;

  // deltan > 0 in compute_force
  // but negative in distance algorithm
  cdata.radi = radius_ ? radius_[iPart] : r0_;
  cdata.r = cdata.radi - cdata.deltan; // sign of corrected, because negative value is passed
  cdata.rsq = cdata.r*cdata.r;
  cdata.meff = rmass_ ? rmass_[iPart] : atom->mass[atom->type[iPart]];
  cdata.area_ratio = 1.;

  ForceData i_forces;
  ForceData j_forces;

  // add to cwl
  if(cwl_ && addflag_)
  {
      double contactPoint[3];
      vectorSubtract3D(x_[cdata.i],cdata.delta,contactPoint);
#ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(atom->superquadric_flag) {
          vectorCopy3D(cdata.contact_point, contactPoint);
      }
#endif
      cwl_->add_wall_1(iMesh,mesh->id(iTri),iPart,contactPoint,v_wall);
  }

  if(impl) {
    impl->compute_force(this, cdata, intersectflag, v_wall, i_forces, j_forces);
  } else {
    double force_old[3]={};

    // if force should be stored - remember old force
    if(store_force_ || stress_flag_)
      vectorCopy3D(f_[iPart],force_old);

    compute_force(cdata, v_wall); // LEGACY CODE (SPH)

    if(store_force_ || stress_flag_)
    {
      vectorSubtract3D(f_[iPart], force_old, i_forces.delta_F);
    }
  }

  // if force should be stored or evaluated
  if(cdata.has_force_update && (store_force_ || stress_flag_))
  {
    if(store_force_)
        vectorAdd3D (wallforce_[iPart], i_forces.delta_F, wallforce_[iPart]);

    if(stress_flag_ && fix_mesh->trackStress())
    {
        double delta[3];
        delta[0] = -cdata.delta[0];
        delta[1] = -cdata.delta[1];
        delta[2] = -cdata.delta[2];
        static_cast<FixMeshSurfaceStress*>(fix_mesh)->add_particle_contribution
        (
           iPart,i_forces.delta_F,delta,iTri,v_wall
        );
    }
  }

  // add heat flux
  if(heattransfer_flag_)
    addHeatFlux(mesh,iPart,cdata.deltan,1.);
}

/* ---------------------------------------------------------------------- */

int FixWallGran::is_moving()
{
    if(is_mesh_wall())
    {
        for(int i = 0; i < n_FixMesh_; i++) {
            if(FixMesh_list_[i]->mesh()->isMoving())
               return 1;
        }
        return 0;
    }
    return shear_;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local()
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts();

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all()
{
    int ncontacts = n_contacts_local();
    MPI_Sum_Scalar(ncontacts,world);
    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_local(int contact_groupbit)
{
    if (!is_mesh_wall() || dnum() == 0) return 0;

    int ncontacts = 0;
    for(int i = 0; i < n_FixMesh_; i++)
        ncontacts += FixMesh_list_[i]->contactHistory()->n_contacts(contact_groupbit);

    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_contacts_all(int contact_groupbit)
{
    int ncontacts = n_contacts_local(contact_groupbit);
    MPI_Sum_Scalar(ncontacts,world);
    return ncontacts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_ghosts_local()
{
    if (!is_mesh_wall()) return 0;

    int nghosts = 0;
    for(int iMesh = 0; iMesh < n_FixMesh_; ++iMesh) {
        TriMesh *mesh = FixMesh_list_[iMesh]->triMesh();
        nghosts += mesh->sizeGhost();
    }

    return nghosts;
}

/* ---------------------------------------------------------------------- */

int FixWallGran::n_ghosts_all()
{
    int nghosts = n_ghosts_local();
    MPI_Sum_Scalar(nghosts,world);
    return nghosts;
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixWallGran::register_compute_wall_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cwl_ != NULL)
     error->fix_error(FLERR,this,"Fix wall/gran allows only one compute of type wall/gran/local");
   cwl_ = ptr;
   dnum_compute = dnum_; //history values
}

void FixWallGran::unregister_compute_wall_local(ComputePairGranLocal *ptr)
{
   if(cwl_ != ptr)
     error->fix_error(FLERR,this,"Illegal situation in FixWallGran::unregister_compute_wall_local");
   cwl_ = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGran::init_heattransfer()
{
    fppa_T = NULL;
    fppa_hf = NULL;
    deltan_ratio = NULL;

    // decide if heat transfer is to be calculated

    if (!is_mesh_wall() && Temp_wall < 0.) return;
    else if (is_mesh_wall())
    {
        int heatflag = 0;
        for(int imesh = 0; imesh < n_meshes(); imesh++)
        {
            heatflag = heatflag || mesh_list()[imesh]->mesh()->prop().getGlobalProperty<ScalarContainer<double> >("Temp") != NULL;
        }

        if(!heatflag) return;
    }

    // heat transfer is to be calculated - continue with initializations

    // set flag so addHeatFlux function is called
    heattransfer_flag_ = true;

    // if(screen && comm->me == 0) fprintf(screen,"Initializing wall/gran heat transfer model\n");
    fppa_T = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",1,0,style));
    fppa_hf = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",1,0,style));

    th_cond = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",0,0,style))->get_values();

    // if youngsModulusOriginal defined, get deltan_ratio
    Fix* ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,style,false);
    // deltan_ratio is defined by heat transfer fix, see if there is one
    int n_htf = modify->n_fixes_style("heat/gran/conduction");
    if (n_htf < 1) n_htf = modify->n_fixes_style_strict("heat/gran");

    // get deltan_ratio set by the heat transfer fix
    if(ymo_fix && n_htf) deltan_ratio = static_cast<FixPropertyGlobal*>(ymo_fix)->get_array_modified();
}

/* ---------------------------------------------------------------------- */

//NP modified C.K.
void FixWallGran::addHeatFlux(TriMesh *mesh,int ip, double delta_n, double area_ratio)
{
    //r is the distance between the sphere center and wall
    double tcop, tcowall, hc, Acont=0.0, r;
    double reff_wall = atom->radius[ip];
    int itype = atom->type[ip];
    double ri = atom->radius[ip];

    if(mesh)
        Temp_wall = (*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp"))(0);

    double *Temp_p = fppa_T->vector_atom;
    double *heatflux = fppa_hf->vector_atom;

    /*NL*/ //if (screen) fprintf(screen,"size %d Temp %f\n",(*mesh->prop().getGlobalProperty< ScalarContainer<double> >("Temp")).size(),Temp_wall);

    if(CONDUCTION_CONTACT_AREA_OVERLAP == area_calculation_mode_)
    {
        //NP adjust overlap that may be superficially large due to softening
        if(deltan_ratio)
           delta_n *= deltan_ratio[itype-1][atom_type_wall_-1];

        r = ri - delta_n;

        Acont = (reff_wall*reff_wall-r*r)*M_PI*area_ratio; //contact area sphere-wall
    }
    else if (CONDUCTION_CONTACT_AREA_CONSTANT == area_calculation_mode_)
        Acont = fixed_contact_area_;
    else if (CONDUCTION_CONTACT_AREA_PROJECTION == area_calculation_mode_)
    {
        Acont = M_PI*ri*ri;
    }

    tcop = th_cond[itype-1]; //types start at 1, array at 0
    tcowall = th_cond[atom_type_wall_-1];

    if ((fabs(tcop) < SMALL) || (fabs(tcowall) < SMALL)) hc = 0.;
    else hc = 4.*tcop*tcowall/(tcop+tcowall)*sqrt(Acont);

    if(computeflag_)
    {
        heatflux[ip] += (Temp_wall-Temp_p[ip]) * hc;
        Q_add += (Temp_wall-Temp_p[ip]) * hc * update->dt;
    }
    if(cwl_ && addflag_)
        cwl_->add_heat_wall(ip,(Temp_wall-Temp_p[ip]) * hc);
    /*NL*/ //if (screen) fprintf(screen,"adding heat flux of %g to particle %d, wall %f Part %f hc %f\n",(Temp_wall-Temp_p[ip]) * hc,ip,Temp_wall,Temp_p[ip],hc);
    /*NL*/ //if (screen) fprintf(screen,"tcop %f tcowall %f Acont %f reff_wall %f r %f delta_n %f\n",tcop,tcowall,Acont,reff_wall,r,delta_n);
}

int64_t FixWallGran::hashcode() {
  return impl->hashcode();
}
