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

#include "math.h"
#include "stdlib.h"
#include "string.h"
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
#include "fix_massflow_mesh_face.h"
#include <algorithm>

using namespace LAMMPS_NS;
using namespace MathExtraLiggghts;
using namespace MathConst;
using namespace FixConst;

enum {UNDEFINED=0, INSIDE, OUTSIDE, IGNORE};

#define FAVRE_AVERAGED
/* ---------------------------------------------------------------------- */

FixMassflowMeshFace::FixMassflowMeshFace(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  once_(false),
  delete_atoms_(false),
  mass_deleted_(0.),
  nparticles_deleted_(0),
  fix_orientation_(0),
  fix_mesh_(0),
  fix_counter_(0),
  fix_prevx_(0),
  fix_neighlist_(0),
  havePointAtOutlet_(false),
  insideOut_(false),
  mass_(0.),
  nparticles_(0),
  average_vx_out_(0.),
  average_vy_out_(0.),
  average_vz_out_(0.),
  fix_property_(0),
  property_sum_(0.),
  screenflag_(false),
  fp_(NULL),
  mass_last_(0.),
  nparticles_last_(0.),
  t_count_(0.),
  delta_t_(0.),
  reset_t_count_(true),
  fix_ms_(0),
  ms_(0),
  ms_counter_(0),
  cg_(1.),
  cg3_(1.)
{
    vectorZeroize3D(nvec_);
    vectorZeroize3D(pref_);
    vectorZeroize3D(pointAtOutlet_);

    // parse args for this class
    int iarg = 3;

    bool hasargs = true;
    while(iarg < narg && hasargs)
    {
        hasargs = false;

        if(strcmp(arg[iarg],"mesh") == 0) {
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
            hasargs = true;
        } else if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");

            char* filecurrent = new char[strlen(arg[iarg+1]) + 8];
            if (1 < comm->nprocs) //open a separate file for each processor
                 sprintf(filecurrent,"%s%s%d",arg[iarg+1],".",comm->me);
            else  //open one file for proc 0
                 sprintf(filecurrent,"%s",arg[iarg+1]);

            if (strcmp(arg[iarg],"file") == 0)
                fp_ = fopen(filecurrent,"w");
            else
                fp_ = fopen(filecurrent,"a");
            if (fp_ == NULL) {
                char str[128];
                sprintf(str,"Cannot open file %s",arg[iarg+1]);
                error->fix_error(FLERR,this,str);
            }
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"screen") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg+1],"yes") == 0) screenflag_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) screenflag_ = false;
            else error->fix_error(FLERR,this,"Illegal screen option");
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"delete_atoms") == 0) {
            if(narg < iarg+2)
                error->fix_error(FLERR,this,"Illegal keyword entry");
            if (strcmp(arg[iarg+1],"yes") == 0) delete_atoms_ = true;
            else if (strcmp(arg[iarg+1],"no") == 0) delete_atoms_ = false;
            else error->fix_error(FLERR,this,"Illegal delete_atoms option");
            iarg += 2;
            hasargs = true;
        } else if (strcmp(arg[iarg],"cg") == 0) { // TODO: remove if cg-fg are in separate simulations
            if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
            cg_ = atof(arg[iarg+1]);
            cg3_ = cg_*cg_*cg_;
            iarg += 2;
            hasargs = true;
        } else if(strcmp(style,"massflow/mesh/face") == 0)
            error->fix_error(FLERR,this,"unknown keyword");
    }

    if(fp_ && 1 < comm->nprocs && 0 == comm->me)
      fprintf(screen,"**FixMassflowMeshFace: > 1 process - will write to multiple files\n");

    // error checks on necessary args

    if(!once_ && delete_atoms_)
        error->fix_error(FLERR,this,"using 'delete_atoms yes' requires 'count once'");
    if(!once_ && havePointAtOutlet_)
        error->fix_error(FLERR,this,"setting 'point_at_outlet' requires 'count once'");
    if(!fix_mesh_)
        error->fix_error(FLERR,this,"expecting keyword 'mesh'");

    restart_global = 1;

    vector_flag = 1;
    size_vector = 10;
    if(fix_property_)
        ++size_vector;
    global_freq = 1; // available always

    array_flag = 1;
    size_array_rows = 1; // rows in global array
    size_array_cols = 10; // columns in global array
}

/* ---------------------------------------------------------------------- */

FixMassflowMeshFace::~FixMassflowMeshFace()
{
   if(fp_) fclose(fp_);
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFace::post_create()
{
    // add per-particle count flag

    if(!fix_counter_)
    {
        const char * fixarg[9];

        sprintf(fixid_,"massflowface_%s",id);
        fixarg[0]=fixid_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=fixid_;
        fixarg[4]="scalar"; // 1 scalar per particle to be registered
        fixarg[5]="yes";    // restart yes
        fixarg[6]="no";     // communicate ghost no
        fixarg[7]="no";     // communicate rev no
        fixarg[8]="0";      // take 0 (UNDEFINED) as default
        modify->add_fix(9,const_cast<char**>(fixarg));

        fix_counter_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid_,"property/atom","scalar",0,0,style));
    }

    if(!fix_prevx_)
    {
        char * var_name = new char[16+strlen(id)];
        sprintf(var_name,"massflowface_x_%s",id);

        const char* fixarg[11];
        fixarg[0]=var_name;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=var_name;
        fixarg[4]="vector";
        fixarg[5]="yes";
        fixarg[6]="no";
        fixarg[7]="no";
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fix_prevx_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
        delete [] var_name;
    }

    // add neighbor list

    fix_neighlist_ = fix_mesh_->createOtherNeighList(igroup,id);

    // need to find multisphere here to be able to add property

    fix_ms_ = static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0));
    if(fix_ms_)
    {
        ms_ = &fix_ms_->data();
        ms_counter_ = ms_->prop().addElementProperty< ScalarContainer<int> >("counter_ms","comm_exchange_borders","frame_invariant", "restart_yes");
        ms_counter_->setDefaultValue(UNDEFINED);

        if(delete_atoms_)
            error->fix_error(FLERR,this,"can not use 'delete_atoms' with fix multisphere");
        if(!once_)
            error->fix_error(FLERR,this,"must use 'count once' with fix multisphere");
    }


    // get face IDs from mesh
    TriMesh *mesh = fix_mesh_->triMesh();

    int nTriLocal = mesh->sizeLocal();
    int *tri_face_ids_local = new int[nTriLocal];

    if(mesh->prop().getElementProperty<ScalarContainer<int> >("face_id"))
    {
        ScalarContainer<int> *face_id = mesh->prop().getElementProperty<ScalarContainer<int> >("face_id");
        for(int iTri = 0; iTri < nTriLocal; ++iTri)
        {
            tri_face_ids_local[iTri] = face_id->get(iTri);
        }
    }
    else
    {
        error->fix_error(FLERR,this,"Mesh element property 'face_id' required");
    }

    int nTriGlobal = nTriLocal;
    int *tri_face_ids_recv = NULL;

    if(mesh->isParallel())
    {
        nTriGlobal = mesh->sizeGlobal();
        MPI_Allgather_Vector(tri_face_ids_local, nTriLocal, tri_face_ids_recv, world);
        delete [] tri_face_ids_local;
    }
    else // every proc has full mesh information
    {
        // NOTE: mesh->sizeGlobal() not initialized yet;
        tri_face_ids_recv = tri_face_ids_local;
    }


    // faceid2index_->first:  unique face id
    // faceid2index_->second: index in std::vector member variables
    for(int iTri = 0; iTri < nTriGlobal; ++iTri)
    {
        faceid2index_[tri_face_ids_recv[iTri]] = 0;
    }

    int i = 0;
    for(std::map<int,int>::iterator it=faceid2index_.begin(); it!=faceid2index_.end(); ++it)
    {
        it->second = i++;
    }

    int nfaceids = faceid2index_.size();
    mass_face_.resize(nfaceids);
    nparticles_face_.resize(nfaceids);
    average_vx_face_out_.resize(nfaceids);
    average_vy_face_out_.resize(nfaceids);
    average_vz_face_out_.resize(nfaceids);
    mass_face_last_.resize(nfaceids);
    nparticles_face_last_.resize(nfaceids);
    distributions_face_.resize(nfaceids);
    // shouldn't be necessary (?):
    std::fill_n(mass_face_.begin(),            nfaceids, 0.);
    std::fill_n(nparticles_face_.begin(),      nfaceids, 0 );
    std::fill_n(average_vx_face_out_.begin(),  nfaceids, 0.);
    std::fill_n(average_vy_face_out_.begin(),  nfaceids, 0.);
    std::fill_n(average_vz_face_out_.begin(),  nfaceids, 0.);
    std::fill_n(mass_face_last_.begin(),       nfaceids, 0.);
    std::fill_n(nparticles_face_last_.begin(), nfaceids, 0 );

    delete [] tri_face_ids_recv;

    size_array_rows = nfaceids; // rows in global array
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFace::pre_delete(bool unfixflag)
{
    if (unfixflag)
    {
        modify->delete_fix(fixid_);
        modify->delete_fix(fix_prevx_->id);
    }
}

/* ----------------------------------------------------------------------
   initialize this fix
------------------------------------------------------------------------- */

void FixMassflowMeshFace::init()
{
    if (atom->rmass_flag == 0)
        error->fix_error(FLERR,this,"requires atoms have mass");

    if(delete_atoms_ && atom->map_style != 1)
        error->fix_error(FLERR,this,"requires an atom map of type 'array', via an 'atom_modify map array' command");

    if(!fix_ms_ && static_cast<FixMultisphere*>(modify->find_fix_style_strict("multisphere",0)))
        error->fix_error(FLERR,this,"fix multisphere must come before fix massflow/mesh in input script");


    // TODO: enable if cg and fg are in separate simulations
    //cg_ = force->cg();
    //cg3_ = cg_*cg_*cg_;
}

/* ---------------------------------------------------------------------- */

void FixMassflowMeshFace::setup(int /*vflag*/)
{
}

/* ---------------------------------------------------------------------- */

int FixMassflowMeshFace::setmask()
{
    int mask = 0;
    mask |= POST_INTEGRATE;
    if(delete_atoms_) mask |= PRE_EXCHANGE;
    return mask;
}

/* ----------------------------------------------------------------------
   evaluate mass which went through face
   all nearby particles
   0 = default value for counter, set if particle is not
       a neighbor of the TriMesh
   1 = was not on nvec_ last step and is thus eligible
   2 = was on nvec_ side last step
   3 = do not re-count, was counted already
------------------------------------------------------------------------- */

void FixMassflowMeshFace::post_integrate()
{
    int nlocal = atom->nlocal;
    double **x = atom->x;
    double **v = atom->v;
    double *radius = atom->radius;
    double *rmass = atom->rmass;
    double dot,delta[3]={};
    double mass_this = 0.;
    int nparticles_this = 0.;
    double average_vx_out_this = 0.0;
    double average_vy_out_this = 0.0;
    double average_vz_out_this = 0.0;
    double property_this = 0.;
    double deltan;
    int *tag = atom->tag;
    double **prevx = fix_prevx_->array_atom;

    class FixPropertyAtom* fix_color=static_cast<FixPropertyAtom*>(modify->find_fix_property("color","property/atom","scalar",0,0,style,false));
    bool fixColFound = false;
    if (fix_color) fixColFound=true;

    fix_orientation_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("ex",
        "property/atom","vector",0,0,style,false));

    TriMesh *mesh = fix_mesh_->triMesh();
    int nTriAll = mesh->sizeLocal() + mesh->sizeGhost();

    ScalarContainer<int> *tri_face_ids = mesh->prop().getElementProperty<ScalarContainer<int> >("face_id");
    int nfaceids = faceid2index_.size();

    // update time for counter
    // also store values for last invokation
    t_count_ += update->dt;

    if(!reset_t_count_)
    {
        average_vx_out_ = 0.;
        average_vy_out_ = 0.;
        average_vz_out_ = 0.;
        nparticles_last_ = nparticles_;
        mass_last_ = mass_;
        reset_t_count_ = true;

        std::fill_n(average_vx_face_out_.begin(),  nfaceids, 0.);
        std::fill_n(average_vy_face_out_.begin(),  nfaceids, 0.);
        std::fill_n(average_vz_face_out_.begin(),  nfaceids, 0.);
        nparticles_face_last_ = nparticles_face_;
        mass_face_last_ = mass_face_;
        distributions_face_.clear();
        distributions_face_.resize(nfaceids);
    }

    std::vector<std::vector<double> > radius_dist_face_local_this(nfaceids);
    std::vector<std::vector<double> > mass_dist_face_local_this(nfaceids);
    std::vector<std::vector<int> > atomtype_dist_face_local_this(nfaceids);
    std::vector<double> mass_face_this(nfaceids, 0.0);
    std::vector<int> nparticles_face_this(nfaceids, 0);
    std::vector<double> average_vx_face_out_this(nfaceids, 0.);
    std::vector<double> average_vy_face_out_this(nfaceids, 0.);
    std::vector<double> average_vz_face_out_this(nfaceids, 0.);

    std::map<int, int> classified_particles_this;
    std::map<int, int> crossing_particles_this;
    std::set<int> ignore_this;
    std::set<int> once_this;
    std::vector<bool> defined_this(nlocal,false); // are in neighlist of a triangle

    // loop owned and ghost triangles
    // count only if owned particle

    for(int iTri = 0; iTri < nTriAll; iTri++)
    {
        int face_id = tri_face_ids->get(iTri);
        mesh->surfaceNorm(iTri,nvec_);
        mesh->node(iTri,0,pref_);

        const std::vector<int> & neighborList = fix_neighlist_->get_contact_list(iTri);
        const int numneigh = neighborList.size();

        for(int iNeigh = 0; iNeigh < numneigh; iNeigh++)
        {
            const int iPart = neighborList[iNeigh];

            // skip ghost particles
            if(iPart >= nlocal) continue;

            defined_this[iPart] = true;

            const int ibody = fix_ms_ ? ( (fix_ms_->belongs_to(iPart) > -1) ? (ms_->map(fix_ms_->belongs_to(iPart))) : -1 ) : -1;

            // in case of once_ == true, ignore everything which has been already counted
            if((ibody > -1) ?  ((*ms_counter_)(ibody) == IGNORE) : (fix_counter_->get_vector_atom_int(iPart) == IGNORE))
                continue;

            if(ignore_this.find(iPart) != ignore_this.end())
              continue;

            if(havePointAtOutlet_)
            {
                // get the vector from the particle center to the next triangle
                deltan = fix_mesh_->triMesh()->resolveTriSphereContact(iPart,iTri,radius[iPart],x[iPart],delta);
                if(deltan < radius[iPart])
                {
                    vectorSubtract3D(x[iPart],pointAtOutlet_,nvec_); // vector pointing to the particle location
                    dot = vectorDot3D(delta,nvec_);
                }
                else // particle is not overlapping with mesh, so continue
                    continue;

                if(insideOut_) dot =-dot;
            }
            else
            {
                vectorSubtract3D(x[iPart],pref_,delta);
                dot = vectorDot3D(delta,nvec_);
                if(insideOut_) dot =-dot;
            }

            // first-time, just set INSIDE or OUTSIDE depending on what side of the mesh
            if(ibody > -1)
            {
                if((*ms_counter_)(ibody) == UNDEFINED)
                    (*ms_counter_)(ibody) = (dot <= 0.) ? INSIDE : OUTSIDE;
                classified_particles_this[iPart] = iTri;
                prevx[iPart][0] = x[iPart][0];
                prevx[iPart][1] = x[iPart][1];
                prevx[iPart][2] = x[iPart][2];
                continue;
            }
            else if(fix_counter_->get_vector_atom_int(iPart) == UNDEFINED)
            {
                fix_counter_->set_vector_atom_int(iPart, (dot <= 0.) ? INSIDE : OUTSIDE);
                classified_particles_this[iPart] = iTri;
                prevx[iPart][0] = x[iPart][0];
                prevx[iPart][1] = x[iPart][1];
                prevx[iPart][2] = x[iPart][2];
                continue;
            }

            // A particle might be in the neighbor list of more than one triangle,
            // for complex geometries it might get tricky to correctly classifiy the particle position
            // the following code is not perfect but covers a broad range of situations
            if(classified_particles_this.find(iPart) != classified_particles_this.end()) // already classified by other triangle
            {
                if( confirm_classification(ibody, iPart, classified_particles_this[iPart], dot > 0. ? OUTSIDE : INSIDE) )
                    ignore_this.insert(iPart);

                continue;
            }

            // WE'RE JUST COUNTING PARTICLES CROSSING BORDERS FROM INSIDE TO OUTSIDE
            if(dot > 0.) // particle is now OUTSIDE
            {
                // particle was not OUTSIDE before
                if((ibody > -1) ? ((*ms_counter_)(ibody) == INSIDE) : (fix_counter_->get_vector_atom_int(iPart) == INSIDE))
                {
                    // make sure that particle really crossed the tri
                    {
                        double dir[3]; MathExtra::sub3(prevx[iPart], x[iPart], dir);
                        double v0[3];  mesh->node(iTri,0,v0);
                        double v1[3];  mesh->node(iTri,1,v1);
                        double v2[3];  mesh->node(iTri,2,v2);

                        if(!MathExtraLiggghts::line_triangle_intersect(x[iPart], dir, v0, v1, v2))
                            continue;
                        ignore_this.insert(iPart);
                    }
                    // each particle holds whole MS mass
                    mass_this += rmass[iPart]; // total mass crossing any mesh face
                    nparticles_this++; // total number of cg particles crossing any mesh face
#ifdef FAVRE_AVERAGED
                    average_vx_out_this += rmass[iPart]*v[iPart][0];
                    average_vy_out_this += rmass[iPart]*v[iPart][1];
                    average_vz_out_this += rmass[iPart]*v[iPart][2];
#else
                    average_vx_out_this += v[iPart][0];
                    average_vy_out_this += v[iPart][1];
                    average_vz_out_this += v[iPart][2];
#endif
                    mass_face_this[faceid2index_[face_id]] += rmass[iPart]; // total mass crossing current face

                    radius_dist_face_local_this  [faceid2index_[face_id]].push_back(radius[iPart]);
                    mass_dist_face_local_this    [faceid2index_[face_id]].push_back(rmass[iPart]);
                    atomtype_dist_face_local_this[faceid2index_[face_id]].push_back(atom->type[iPart]);

                    nparticles_face_this[faceid2index_[face_id]]++; // number of cg particles crossing current face
#ifdef FAVRE_AVERAGED
                    average_vx_face_out_this[faceid2index_[face_id]] += rmass[iPart]*v[iPart][0];
                    average_vy_face_out_this[faceid2index_[face_id]] += rmass[iPart]*v[iPart][1];
                    average_vz_face_out_this[faceid2index_[face_id]] += rmass[iPart]*v[iPart][2];
#else
                    average_vx_face_out_this[faceid2index_[face_id]] += v[iPart][0];
                    average_vy_face_out_this[faceid2index_[face_id]] += v[iPart][1];
                    average_vz_face_out_this[faceid2index_[face_id]] += v[iPart][2];
#endif
                    crossing_particles_this[iPart] = iTri;

                    if(fix_property_)
                    {
                        property_this += fix_property_->vector_atom[iPart];
                    }

                    if(delete_atoms_)
                    {
                        atom_tags_delete_.push_back(atom->tag[iPart]);
                    }

                    if(ibody > -1)
                        (*ms_counter_)(ibody) = OUTSIDE;
                    else
                        fix_counter_->set_vector_atom_int(iPart, OUTSIDE);

                    if(once_) once_this.insert(iPart);

                }
            }
            else // dot <= 0 INSIDE
            {
                if((ibody > -1) ? ((*ms_counter_)(ibody) == OUTSIDE) : (fix_counter_->get_vector_atom_int(iPart) == OUTSIDE))
                {
                    // make sure that particle really crossed the tri
                    {
                        double dir[3]; MathExtra::sub3(prevx[iPart], x[iPart], dir);
                        double v0[3];  mesh->node(iTri,0,v0);
                        double v1[3];  mesh->node(iTri,1,v1);
                        double v2[3];  mesh->node(iTri,2,v2);

                        if(!MathExtraLiggghts::line_triangle_intersect(x[iPart], dir, v0, v1, v2))
                            continue;
                    }

                    ignore_this.insert(iPart);

                    if(ibody > -1)
                        (*ms_counter_)(ibody) = INSIDE;
                    else
                        fix_counter_->set_vector_atom_int(iPart, INSIDE);
                }
            }
        }
    }

    if(once_)
    {
        for(std::set<int>::iterator it=once_this.begin(); it!=once_this.end(); ++it)
        {
            const int ibody = fix_ms_ ? ( (fix_ms_->belongs_to(*it) > -1) ? (ms_->map(fix_ms_->belongs_to(*it))) : -1 ) : -1;
            if(ibody > -1)
                (*ms_counter_)(ibody) = IGNORE;
            else
                fix_counter_->set_vector_atom_int(*it, IGNORE);
        }
    }

    for(int iPart = 0; iPart < nlocal; ++iPart)
    {
        prevx[iPart][0] = x[iPart][0];
        prevx[iPart][1] = x[iPart][1];
        prevx[iPart][2] = x[iPart][2];

        if(!defined_this[iPart])
        {
            const int ibody = fix_ms_ ? ( (fix_ms_->belongs_to(iPart) > -1) ? (ms_->map(fix_ms_->belongs_to(iPart))) : -1 ) : -1;
            if(ibody > -1)
            {
                if((*ms_counter_)(ibody) != IGNORE)
                    (*ms_counter_)(ibody) = UNDEFINED;
            }
            else
            {
                if(fix_counter_->get_vector_atom_int(iPart) != IGNORE)
                {
                    fix_counter_->set_vector_atom_int(iPart, UNDEFINED);
                }
            }
        }
    }

    if(screenflag_ && screen)
    {
        for(std::map<int,int>::iterator it=crossing_particles_this.begin(); it!=crossing_particles_this.end(); ++it)
        {
            int iPart = it->first;
            fprintf(screen,"%ld %d %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g \n ",
                       update->ntimestep, tri_face_ids->get(it->second), tag[iPart],
                       2.*radius[iPart]/force->cg(),
                       x[iPart][0],x[iPart][1],x[iPart][2],
                       v[iPart][0],v[iPart][1],v[iPart][2]);
        }
    }

    if(fp_)
    {
        for(std::map<int,int>::iterator it=crossing_particles_this.begin(); it!=crossing_particles_this.end(); ++it)
        {
            int iPart = it->first;
            fprintf(fp_," %ld %d %d %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g %4.4g",
                       update->ntimestep, tri_face_ids->get(it->second), tag[iPart],
                       2.*radius[iPart]/force->cg(),
                       x[iPart][0],x[iPart][1],x[iPart][2],
                       v[iPart][0],v[iPart][1],v[iPart][2]);
            if (fixColFound)
                fprintf(fp_,"    %4.0g ", fix_color->vector_atom[iPart]);

            if(fix_orientation_)
            {
                double **orientation = NULL;
                orientation = fix_orientation_->array_atom;
                fprintf(fp_,"    %4.4g %4.4g %4.4g ",
                        orientation[iPart][0], orientation[iPart][1], orientation[iPart][2]);
            }
            fprintf(fp_,"\n");
            fflush(fp_);
        }
    }

    // sum global values from all processes
    MPI_Sum_Scalar(mass_this,world);
    MPI_Sum_Scalar(nparticles_this,world);
    MPI_Sum_Scalar(average_vx_out_this,world);
    MPI_Sum_Scalar(average_vy_out_this,world);
    MPI_Sum_Scalar(average_vz_out_this,world);
    if(fix_property_) MPI_Sum_Scalar(property_this,world);


    if(nparticles_this > 0) // some particles crossed the mesh faces
    {
        double *radius_dist_this_recv = NULL; // gather from all procs
        double *mass_dist_this_recv = NULL;
        int *atomtype_dist_this_recv = NULL;
        int nparticles_face_i_this_all = 0;

        for(int i=0; i<nfaceids; ++i)
        {
            MPI_Sum_Scalar(nparticles_face_this[i], nparticles_face_i_this_all, world);

            if(nparticles_face_i_this_all > 0)
            {
                MPI_Allgather_Vector(&radius_dist_face_local_this  [i][0], nparticles_face_this[i], radius_dist_this_recv, world);
                MPI_Allgather_Vector(&mass_dist_face_local_this    [i][0], nparticles_face_this[i], mass_dist_this_recv, world);
                MPI_Allgather_Vector(&atomtype_dist_face_local_this[i][0], nparticles_face_this[i], atomtype_dist_this_recv, world);

                for(int p = 0; p < nparticles_face_i_this_all; ++p)
                {
                    double rad  = radius_dist_this_recv[p]/cg_; // fg particle radius
                    double mass = mass_dist_this_recv[p]/cg3_;  // fg particle mass
                    int type    = atomtype_dist_this_recv[p];
                    distributions_face_[i][ConstantParticleTemplateSphere(rad, mass, type)] += cg3_; // number of fg particles!
                }

                delete [] radius_dist_this_recv; radius_dist_this_recv = NULL;
                delete [] mass_dist_this_recv; mass_dist_this_recv = NULL;
                delete [] atomtype_dist_this_recv; atomtype_dist_this_recv = NULL;
                nparticles_face_i_this_all = 0;
            }
        }
    }

    // sum per face values from all processes
    MPI_Sum_Vector(&mass_face_this[0],nfaceids,world);
    MPI_Sum_Vector(&nparticles_face_this[0],nfaceids,world);
    MPI_Sum_Vector(&average_vx_face_out_this[0],nfaceids,world);
    MPI_Sum_Vector(&average_vy_face_out_this[0],nfaceids,world);
    MPI_Sum_Vector(&average_vz_face_out_this[0],nfaceids,world);

    mass_ += mass_this;
    nparticles_ += nparticles_this; // total number of cg particles that crossed any mesh face since last inquiry
    average_vx_out_ += average_vx_out_this;
    average_vy_out_ += average_vy_out_this;
    average_vz_out_ += average_vz_out_this;
    property_sum_ += property_this;

    // global, same values across all processors
    for(int i=0; i<nfaceids; ++i)
    {
        mass_face_[i] += mass_face_this[i];
        nparticles_face_[i] += nparticles_face_this[i];
        average_vx_face_out_[i] += average_vx_face_out_this[i];
        average_vy_face_out_[i] += average_vy_face_out_this[i];
        average_vz_face_out_[i] += average_vz_face_out_this[i];
    }
}

/* ----------------------------------------------------------------------
   perform particle deletion of marked particles
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct
------------------------------------------------------------------------- */

void FixMassflowMeshFace::pre_exchange()
{
    if(delete_atoms_)
    {
        double mass_deleted_this_ = 0.;
        int nparticles_deleted_this_ = 0.;
        int *atom_map_array = atom->get_map_array();

        // delete particles

        while (atom_tags_delete_.size() > 0)
        {
            int iPart = atom->map(atom_tags_delete_[0]);

            mass_deleted_this_ += atom->rmass[iPart];
            nparticles_deleted_this_++;

            atom->avec->copy(atom->nlocal-1,iPart,1);

            // manipulate atom map array
            // need to do this since atom map is needed for deletion
            atom_map_array[atom->tag[atom->nlocal-1]] = iPart;

            atom->nlocal--;

            atom_tags_delete_.erase(atom_tags_delete_.begin());
        }

        atom_tags_delete_.clear();

        MPI_Sum_Scalar(mass_deleted_this_,world);
        MPI_Sum_Scalar(nparticles_deleted_this_,world);

        mass_deleted_ += mass_deleted_this_;
        nparticles_deleted_ += nparticles_deleted_this_;

        // tags and maps
        if(nparticles_deleted_this_)
        {
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

void FixMassflowMeshFace::write_restart(FILE *fp)
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

void FixMassflowMeshFace::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  mass_ = list[n++];
  t_count_ = list[n++];
  mass_last_ = list[n++];
  nparticles_last_ = static_cast<int>(list[n++]);
  mass_deleted_ = list[n++];
  nparticles_deleted_ = static_cast<int>(list[n++]);
}

/* ----------------------------------------------------------------------
   output # of faces
------------------------------------------------------------------------- */

double FixMassflowMeshFace::compute_scalar()
{
    return faceid2index_.size();
}

/* ----------------------------------------------------------------------
   always returns fg number of particles (also when coarse-grained!)
------------------------------------------------------------------------- */

double FixMassflowMeshFace::compute_vector(int index)
{
    if(reset_t_count_)
    {
        delta_t_ = t_count_;
        t_count_ = 0.;
        reset_t_count_ = false;
    }

    switch(index)
    {
    case 0:
        return mass_;
    case 1:
        return cg3_*static_cast<double>(nparticles_); // total # of particles passed through mesh since sim start
    case 2:
        return delta_t_ == 0. ? 0. : (mass_-mass_last_)/delta_t_;
    case 3:
        return delta_t_ == 0. ? 0. : cg3_*static_cast<double>(nparticles_-nparticles_last_)/delta_t_;
    case 4:
        return mass_deleted_;
    case 5:
        return cg3_*static_cast<double>(nparticles_deleted_);
    case 6:
        return cg3_*static_cast<double>(nparticles_-nparticles_last_);
    case 7:
#ifdef FAVRE_AVERAGED
        if (mass_ > 0.)
            return average_vx_out_/mass_;
#else
        if (nparticles_ > nparticles_last_)
            return average_vx_out_/(nparticles_-nparticles_last_);
#endif
        return 0.;
    case 8:
#ifdef FAVRE_AVERAGED
        if (mass_ > 0.)
            return average_vy_out_/mass_;
#else
        if (nparticles_ > nparticles_last_)
            return average_vy_out_/(nparticles_-nparticles_last_);
#endif
        return 0.;
    case 9:
#ifdef FAVRE_AVERAGED
        if (mass_ > 0.)
            return average_vz_out_/mass_;
#else
        if (nparticles_ > nparticles_last_)
            return average_vz_out_/(nparticles_-nparticles_last_);
#endif
        return 0.;
    case 10:
        if(fix_property_)
            return property_sum_;
    default:
        return 0.;
    }
}


double FixMassflowMeshFace::compute_array(int i, int j)
{
  if(reset_t_count_)
  {
      delta_t_ = t_count_;
      t_count_ = 0.;
      reset_t_count_ = false;
  }

  switch(j)
  {
  case 0:
      return mass_face_[i];
  case 1:
      return cg3_*static_cast<double>(nparticles_face_[i]);
  case 2:
      return delta_t_ == 0. ? 0. : (mass_face_[i]-mass_face_last_[i])/delta_t_;
  case 3:
      return delta_t_ == 0. ? 0. : cg3_*static_cast<double>(nparticles_face_[i]-nparticles_face_last_[i])/delta_t_;
  case 4:
      return mass_deleted_;
  case 5:
      return cg3_*static_cast<double>(nparticles_deleted_);
  case 6:
      return cg3_*static_cast<double>(nparticles_face_[i]-nparticles_face_last_[i]);
  case 7:
#ifdef FAVRE_AVERAGED
      if (mass_face_[i] > 0.)
          return average_vx_face_out_[i]/mass_face_[i];
#else
      if (nparticles_face_[i]-nparticles_face_last_[i])
          return average_vx_face_out_[i]/(nparticles_face_[i]-nparticles_face_last_[i]);
#endif
      return 0.;
  case 8:
#ifdef FAVRE_AVERAGED
      if (mass_face_[i] > 0.)
          return average_vy_face_out_[i]/mass_face_[i];
#else
      if (nparticles_face_[i]-nparticles_face_last_[i])
          return average_vy_face_out_[i]/(nparticles_face_[i]-nparticles_face_last_[i]);
#endif
      return 0.;
  case 9:
#ifdef FAVRE_AVERAGED
      if (mass_face_[i] > 0.)
          return average_vz_face_out_[i]/mass_face_[i];
#else
      if (nparticles_face_[i]-nparticles_face_last_[i])
          return average_vz_face_out_[i]/(nparticles_face_[i]-nparticles_face_last_[i]);
#endif
      return 0.;
  default:
      return 0.;
  }
}


double FixMassflowMeshFace::compute_array_by_id(int face_id, int j)
{
  if(faceid2index_.find(face_id) == faceid2index_.end())
      error->fix_error(FLERR, this, "Invalid face id!");

  return compute_array(faceid2index_[face_id], j);
}


const DiscreteParticleDistribution& FixMassflowMeshFace::get_distribution_by_id(int face_id)
{
  if(faceid2index_.find(face_id) == faceid2index_.end())
    error->fix_error(FLERR, this, "Invalid face id!");

  return distributions_face_[faceid2index_[face_id]];
}


bool FixMassflowMeshFace::confirm_classification(int ibody, int iPart, int tri, int current_side)
{
  int opposite_side = current_side == OUTSIDE ? INSIDE : OUTSIDE;

  if((ibody > -1) ? ((*ms_counter_)(ibody) == opposite_side) : (fix_counter_->get_vector_atom_int(iPart) == opposite_side))
  {
    double pref_alt[3];
    double delta[3];
    bool convex = true;
    TriMesh *mesh = fix_mesh_->triMesh();

    for(int i = 0; i < 3; ++i)
    {
      mesh->node(tri,i,pref_alt);
      double dot_tri = 0.;
      vectorSubtract3D(pref_alt,pref_,delta);
      dot_tri = vectorDot3D(delta,nvec_);

      if(insideOut_) dot_tri =-dot_tri;

      if(dot_tri > 0.)
      {
        convex = false;
        break;
      }
    }

    if(convex ? (current_side == OUTSIDE) : (current_side == INSIDE))
    {
      if (ibody > -1) (*ms_counter_)(ibody) = current_side;
      else            fix_counter_->set_vector_atom_int(iPart, current_side);

      return true;
    }
  }

  return false;
}
