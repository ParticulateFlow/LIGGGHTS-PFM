/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Copyright 2009-2012 JKU Linz
   Copyright 2012-2014 DCS Computing GmbH, Linz
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
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fix_insert_stream.h"
#include "fix_mesh_surface.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "vector_liggghts.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_property_atom_tracer_stream.h"
#include "fix_particledistribution_discrete.h"
#include "fix_multisphere.h"
#include "multisphere.h"
#include "fix_template_sphere.h"
#include "particleToInsert.h"
#include "tri_mesh_planar.h"

enum{FACE_NONE,FACE_MESH,FACE_CIRCLE};

using namespace LAMMPS_NS;
using namespace FixConst;

/*NL*/ #define LMP_DEBUGMODE_FIXINSERT_STREAM false
/*NL*/ #define LMP_DEBUG_OUT_FIXINSERT_STREAM screen

#define FIX_INSERT_NTRY_SUBBOX 500
#define FIX_INSERT_STREAM_TINY 1e-14

/* ---------------------------------------------------------------------- */

FixInsertStream::FixInsertStream(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"insertion_face") == 0)
    {
      //NP should be possible to define either mesh or face
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      int f_i = modify->find_fix(arg[iarg+1]);
      if (f_i == -1) error->fix_error(FLERR,this,"Could not find fix mesh/surface id you provided");
      if (strncmp(modify->fix[f_i]->style,"mesh",4))
        error->fix_error(FLERR,this,"The fix belonging to the id you provided is not of type mesh");
      ins_face = (static_cast<FixMeshSurface*>(modify->fix[f_i]))->triMesh();
      ins_face->useAsInsertionMesh(false);
      face_style = FACE_MESH;
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"extrude_length") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      extrude_length = atof(arg[iarg+1]);
      if(extrude_length < 0. ) error->fix_error(FLERR,this,"invalid extrude_length");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"duration") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      duration = atoi(arg[iarg+1]);
      if(duration < 1 ) error->fix_error(FLERR,this,"'duration' can not be < 1");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"parallel") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      if(strcmp("yes",arg[iarg+1]) == 0)
        parallel = true;
      else if(strcmp("no",arg[iarg+1]) == 0)
        parallel = false;
      else error->fix_error(FLERR,this,"expecting 'yes' or 'no' for 'parallel'");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"ntry_mc") == 0) {
      if (iarg+2 > narg) error->fix_error(FLERR,this,"not enough arguments");
      ntry_mc = atoi(arg[iarg+1]);
      if(ntry_mc < 1000) error->fix_error(FLERR,this,"ntry_mc must be > 1000");
      iarg += 2;
      hasargs = true;
    } else error->fix_error(FLERR,this,"unknown keyword or wrong keyword order");
  }

  fix_release = NULL;
  i_am_integrator = false;

  tracer = NULL;
  ntracer = 0;

  ins_fraction = 0.;
  do_ins_fraction_calc = true;

  //NP execute end of step
  nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixInsertStream::~FixInsertStream()
{
    if(tracer) delete []tracer;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::post_create()
{
  FixInsert::post_create();

  // only register property if I am the first fix/insert/stream in the simulation
  //NP 14 values: original position to integrate (3), start step (1),
  //NP            release step (1), integration velocity (3), vel(3), omega (3)

  if(modify->n_fixes_style(style) == 1)
  {
        const char* fixarg[22];
        fixarg[0]="release_fix_insert_stream";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="release_fix_insert_stream";
        fixarg[4]="vector"; //NP 1 scalar per particle to be registered
        fixarg[5]="yes";    //NP restart yes
        fixarg[6]="yes";    //NP communicate ghost no
        fixarg[7]="no";    //NP communicate rev yes
        fixarg[8]="0.";     // x
        fixarg[9]="0.";     // y
        fixarg[10]="0.";    // z
        fixarg[11]="0.";    // insertion step
        fixarg[12]="0.";    // release_step
        fixarg[13]="0.";    // v_integrate_x
        fixarg[14]="0.";    // v_integrate_y
        fixarg[15]="0.";    // v_integrate_z
        fixarg[16]="0.";    // v_insert_x
        fixarg[17]="0.";    // v_insert_y
        fixarg[18]="0.";    // v_insert_z
        fixarg[19]="0.";    // omega_x
        fixarg[20]="0.";    // omega_y
        fixarg[21]="0.";    // omega_z
        modify->add_fix_property_atom(22,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::pre_delete(bool unfixflag)
{
    // delete if I am the last fix of this style to be deleted
    if(unfixflag && modify->n_fixes_style(style) == 1)
        modify->delete_fix("release_fix_insert_stream");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init_defaults()
{
    face_style = FACE_NONE;
    extrude_length = 0.;

    extrude_length_min = extrude_length_max = 0.;

    duration = 0;

    parallel = false;

    ntry_mc = 100000;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::register_tracer_callback(FixPropertyAtomTracerStream* tr)
{
    // just return if I already have this callback
    for(int i = 0; i < ntracer; i++)
        if(tracer[i] == tr) return;

    FixPropertyAtomTracerStream** tracer_new = new FixPropertyAtomTracerStream*[ntracer+1];

    for(int i = 0; i < ntracer; i++)
        tracer_new[i] = tracer[i];

    tracer_new[ntracer] = tr;
    ntracer++;
    delete []tracer;
    tracer = tracer_new;
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixInsertStream::calc_insertion_properties()
{
    double dt,dot,extrude_vec[3],t1[3],t2[3];

    //NP do NOT calc insertion fraction here since
    //NP data copy done at init() only

    // error check on insertion face
    if(face_style == FACE_NONE)
        error->fix_error(FLERR,this,"must define an insertion face");

    // check properties of insertion face
    if(face_style == FACE_MESH)
    {
        // check if face planar
        if(!ins_face->isPlanar())
            error->fix_error(FLERR,this,"command requires a planar face for insertion");

        //NP check for /insertion mesh
        if(all_in_flag)
        {
            if(!dynamic_cast<TriMeshPlanar*>(ins_face))
                error->fix_error(FLERR,this,"using all_in yes requires you to use a fix mesh/surface/planar");
        }

        // get normal vector of face 0
        ins_face->surfaceNorm(0,normalvec);

        //NP if (screen) printVec3D(screen,"normalvec",normalvec);

        // flip normal vector so dot product with v_insert is > 0
        dot = vectorDot3D(v_insert,normalvec);
        if(dot < 0) vectorFlip3D(normalvec);

        // calc v normal
        dot = vectorDot3D(v_insert,normalvec);
        vectorCopy3D(normalvec,v_normal);
        vectorScalarMult3D(v_normal,dot);

        // error check on v normal
        if(vectorMag3D(v_normal) < 1.e-3)
          error->fix_error(FLERR,this,"insertion velocity projected on face normal is < 1e-3");

        // get reference point on face
        ins_face->node(0,0,p_ref);
    }
    else error->fix_error(FLERR,this,"FixInsertStream::calc_insertion_properties(): Implementation missing");

    // error check on insertion velocity
    if(vectorMag3D(v_insert) < 1e-5)
        error->fix_error(FLERR,this,"insertion velocity too low");

    // further error-checks
    if(insert_every == -1 && extrude_length == 0.)
      error->fix_error(FLERR,this,"must define either 'insert_every' or 'extrude_length'");
    if(insert_every > -1 && extrude_length > 0.)
      error->fix_error(FLERR,this,"must not provide both 'insert_every' and 'extrude_length'");
    if(extrude_length > 0. && duration > 0)
      error->fix_error(FLERR,this,"must not provide both 'extrude_length' and 'duration'");

    dt = update->dt;

    // if extrude_length given, calculate insert_every
    if(insert_every == -1)
    {
        // no duration allowed here (checked before)

        if(extrude_length < 3.*max_r_bound() && (all_in_flag || check_ol_flag))
            error->fix_error(FLERR,this,"'extrude_length' is too small");
        // add TINY for resolving round-off
        insert_every = static_cast<int>((extrude_length+FIX_INSERT_STREAM_TINY)/(dt*vectorMag3D(v_normal)));
        /*NL*///if (screen) fprintf(screen,"insert_every %d, extrude_length %f, dt %f, vectorMag3D(v_normal) %f\n",insert_every,extrude_length,dt,vectorMag3D(v_normal));
        if(insert_every == 0)
          error->fix_error(FLERR,this,"insertion velocity too high or extrude_length too low");
    }
    // if insert_every given, calculate extrude_length
    // take into account duration can be != insert_every
    else
    {
        if(insert_every < 1) error->fix_error(FLERR,this,"'insert_every' must be > 0");

        // duration = insert_every by default (if already > 0, defined directly)
        if(duration == 0) duration = insert_every;
        else if (duration > insert_every) error->fix_error(FLERR,this,"'duration' > 'insert_every' not allowed");

        extrude_length = static_cast<double>(duration) * dt * vectorMag3D(v_normal);
        /*NL*/ //if (screen) fprintf(screen,"extrude_length %f, max_r_bound() %f, duration %d\n",extrude_length,max_r_bound(),duration);
        if(extrude_length < 3.*max_r_bound())
          error->fix_error(FLERR,this,"'insert_every' or 'vel' is too small, or radius of inserted particles too large");
    }

    /*NL*///if (screen) fprintf(screen,"insert_every %d, duration %d, v_normal %f, extrude_length %f dt %e\n",insert_every,duration,vectorMag3D(v_normal),extrude_length,dt);

    // ninsert - if ninsert not defined directly, calculate it
    if(ninsert == 0 && ninsert_exists)
    {
        if(massinsert > 0.) ninsert = static_cast<int>((massinsert+FIX_INSERT_STREAM_TINY) / fix_distribution->mass_expect());
        else error->fix_error(FLERR,this,"must define either 'nparticles' or 'mass'");
    }

    // flow rate
    if(nflowrate == 0.)
    {
        if(massflowrate == 0.) error->fix_error(FLERR,this,"must define either 'massrate' or 'particlerate'");
        nflowrate = massflowrate / fix_distribution->mass_expect();
    }
    else massflowrate = nflowrate * fix_distribution->mass_expect();

    // ninsert_per and massinsert
    ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);
    if(ninsert_exists) massinsert = static_cast<double>(ninsert) * fix_distribution->mass_expect();

    // calculate bounding box of extruded face
    if(face_style == FACE_MESH)
    {
        // get bounding box for face
        ins_face->getGlobalBoundingBox().getBoxBounds(ins_vol_xmin,ins_vol_xmax);

        // get bounding box for extruded face - store in t1,t2
        vectorScalarMult3D(normalvec,-extrude_length,extrude_vec);
        vectorAdd3D(ins_vol_xmin,extrude_vec,t1);
        vectorAdd3D(ins_vol_xmax,extrude_vec,t2);

        // take min and max
        vectorComponentMin3D(ins_vol_xmin,t1,ins_vol_xmin);
        vectorComponentMax3D(ins_vol_xmax,t2,ins_vol_xmax);

        /*NL*///if (screen) printVec3D(screen,"ins_vol_xmin",ins_vol_xmin);
        /*NL*///if (screen) printVec3D(screen,"ins_vol_xmin",ins_vol_xmax);
    }
    else error->fix_error(FLERR,this,"Missing implementation in calc_insertion_properties()");

    extrude_length_min = 0.;
    extrude_length_max = extrude_length;
}

/* ---------------------------------------------------------------------- */

int FixInsertStream::setmask()
{
    int mask = FixInsert::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init()
{
    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::init() start\n");

    FixInsert::init();

    if(fix_multisphere && v_randomSetting != RANDOM_CONSTANT)
        error->fix_error(FLERR,this,"Currently only fix insert/stream with multisphere particles only supports constant velocity");

    fix_release = static_cast<FixPropertyAtom*>(modify->find_fix_property("release_fix_insert_stream","property/atom","vector",5,0,style));
    if(!fix_release) error->fix_error(FLERR,this,"Internal error if fix insert/stream");

    i_am_integrator = modify->i_am_first_of_style(this);

    //NP currently only one fix rigid allowed
    /*NP
    if(fix_multisphere) fix_multisphere->set_v_integrate(v_normal);
    if(!i_am_integrator && fix_multisphere)
        error->fix_error(FLERR,this,"Currently only one fix insert/stream is allowed with multisphere particles");
    */

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::init() end\n");

    // error check on insertion face
    if(face_style == FACE_NONE)
        error->fix_error(FLERR,this,"must define an insertion face");

    //NP disallow movement because shallow mesh copy would not know about it
    if(ins_face->isMoving() || ins_face->isScaling())
        error->fix_error(FLERR,this,"cannot translate, rotate, scale mesh which is used for particle insertion");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::setup_pre_exchange()
{

}

/* ---------------------------------------------------------------------- */

double FixInsertStream::insertion_fraction()
{
    /*NL*/ if(0 == comm->me && screen && ins_face->isMoving()) fprintf(screen,"Ins face MOVING\n");
    // have to re-calculate insertion fraction for my subbox
    // in case subdomains of simulation box are changing
    //NP ATTENTION domain->box_change_domain is false even if box is changing size, i.e.
    //NP subboxes are changing size
    if(domain->box_change || do_ins_fraction_calc || ins_face->isMoving())
        calc_ins_fraction();

    return ins_fraction;
}

/* ----------------------------------------------------------------------
   calculate insertion fraction for my subbox
   has to be called at initialization and before every insertion in case
   box is changing
------------------------------------------------------------------------- */

void FixInsertStream::calc_ins_fraction()
{
    //NP similar to Region::volume_mc

    do_ins_fraction_calc = false;

    double pos[3], boxedgevec[3], dot;
    int n_in_local = 0, n_test = ntry_mc;

    for(int i = 0; i < n_test; i++)
    {
        generate_random_global(pos);

        if(domain->is_in_subdomain(pos))
           n_in_local++;
    }

    ins_fraction = static_cast<double>(n_in_local)/static_cast<double>(n_test);

    // also calculate min and max extrusion
    // this can speed up insertion if extrusion volume extends across multiple procs

    if(parallel)
    {
        extrude_length_min = extrude_length;
        extrude_length_max = 0.;

        for(int ix = 0; ix < 2; ix++)
            for(int iy = 0; iy < 2; iy++)
                for(int iz = 0; iz < 2; iz++)
                {
                    vectorConstruct3D
                    (
                        boxedgevec,
                        (ix == 0 ? domain->sublo[0] : domain->subhi[0]) - p_ref[0],
                        (iy == 0 ? domain->sublo[1] : domain->subhi[1]) - p_ref[1],
                        (iz == 0 ? domain->sublo[2] : domain->subhi[2]) - p_ref[2]
                    );

                    dot = -vectorDot3D(boxedgevec,normalvec);
                    /*NL*/ //if (screen) fprintf(screen,"dot %f\n",dot);
                    if(dot > 0. && dot < extrude_length)
                    {
                        extrude_length_max = MathExtraLiggghts::max(extrude_length_max,dot);
                        extrude_length_min = MathExtraLiggghts::min(extrude_length_min,dot);
                    }
                    else if(dot < 0.)
                        extrude_length_min = 0.;
                    else if(dot >= extrude_length)
                        extrude_length_max = extrude_length;
                }
        if(extrude_length_min == extrude_length)
            extrude_length_min = 0.;
        if(extrude_length_max == 0.)
            extrude_length_max = extrude_length;
    }

    /*NL*/// if (screen) fprintf(screen,"proc %d: fraction %f extrude_length_min %f extrude_length_max %f\n",
    /*NL*///         comm->me,ins_fraction,extrude_length_min,extrude_length_max);

    double ins_fraction_all;
    MPI_Sum_Scalar(ins_fraction,ins_fraction_all,world);
    if(ins_fraction_all < 0.9 || ins_fraction_all > 1.1)
        error->fix_error(FLERR,this,"insertion volume could not be distributed properly in parallel. "
                                     "Bad decomposition or insertion face extrusion is too small or outside domain");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::pre_insert()
{
    if((!domain->is_in_domain(ins_vol_xmin) || !domain->is_in_domain(ins_vol_xmax)) && comm->me == 0)
      error->warning(FLERR,"Fix insert/stream: Extruded insertion face extends outside domain, may not insert all particles correctly");

    //NP should do error check as in FixInsertPack::calc_ninsert_this() here
}

/* ---------------------------------------------------------------------- */

inline int FixInsertStream::is_nearby(int i)
{
    double pos_rel[3], pos_projected[3], t[3];
    double **x = atom->x;

    /*NL*///if(screen && atom->tag[i] == 230) fprintf(screen,"pref = %f %f %f\n",p_ref[0],p_ref[1],p_ref[2]);

    vectorSubtract3D(x[i],p_ref,pos_rel);
    double dist_normal = vectorDot3D(pos_rel,normalvec);

    /*NL*///if(screen && atom->tag[i] == 230) fprintf(screen,"pos_rel %f %f %f\n",pos_rel[0],pos_rel[1],pos_rel[2]);
    /*NL*///if(screen && atom->tag[i] == 230) fprintf(screen,"dist_normal %f  extrude_length %f maxrad %f\n",dist_normal,extrude_length,maxrad);

    // on wrong side of extrusion
    if(dist_normal > maxrad) return 0;

    /*NL*///if (screen) fprintf(screen,"1\n");

    // on right side of extrusion, but too far away
    // 3*maxrad b/c extrude_length+rad is max extrusion for overlapcheck yes
    if(dist_normal < -(extrude_length + 3.*maxrad)) return 0;

    /*NL*///if (screen) fprintf(screen,"2\n");

    // on right side of extrusion, within extrude_length
    // check if projection is on face or not

    vectorScalarMult3D(normalvec,dist_normal,t);
    vectorAdd3D(x[i],t,pos_projected);

    //TODO also should check if projected point is NEAR surface

    return ins_face->isOnSurface(pos_projected);
}

BoundingBox FixInsertStream::getBoundingBox() const {
  BoundingBox bb = ins_face->getGlobalBoundingBox();

  const double cut = 3.*maxrad;
  const double delta = -(extrude_length + 2.*cut);
  bb.extrude(delta, normalvec);
  bb.shrinkToSubbox(domain->sublo, domain->subhi);

  // extend to include ghost particles
  const double extend = cut + extend_cut_ghost();
  bb.extendByDelta(extend);

  return bb;
}

/* ----------------------------------------------------------------------
   generate random positions on insertion face
   extrude by random length in negative face normal direction
     currently only implemented for all_in_flag = 0
     since would be tedious to check/implement otherwise
------------------------------------------------------------------------- */

inline void FixInsertStream::generate_random(double *pos, double rad)
{
    double r, ext[3];

    // generate random position on the mesh
    //NP position is _not_ restricted to my subbox
    //NP otherwise parallelization would not work if extrusion volume is
    //NP distributed across multiple processors over height
    if(all_in_flag)
        ins_face->generateRandomOwnedGhostWithin(pos,rad);
        //NP ins_face->generateRandomSubboxWithin(pos,rad);
    else
        ins_face->generateRandomOwnedGhost(pos);
        //NP ins_face->generateRandomSubbox(pos);

    // extrude the position
    //NP min pos: max_rbound,
    //NP max pos: extrude_length - max_rbound or extrude_length
    //NP   extrude_length - max_rbound for strict non-overlapping
    //NP   if overlap is checked for, can go a bit beyond so
    //NP   stream is more continuous

    if(check_ol_flag)
        r = -1.*(random->uniform()*(extrude_length_max         ) + rad + extrude_length_min);
    else
        r = -1.*(random->uniform()*(extrude_length_max - 2.*rad) + rad + extrude_length_min);

    vectorScalarMult3D(normalvec,r,ext);
    vectorAdd3D(pos,ext,pos);
}

/* ----------------------------------------------------------------------
   generate random positions on shallow copy insertion face
   extrude by random length in negative face normal direction
     currently only implemented for all_in_flag = 0
     since would be tedious to check/implement otherwise
------------------------------------------------------------------------- */

inline void FixInsertStream::generate_random_global(double *pos)
{
    double r, ext[3];

    // generate random position on the mesh
    //NP position is not restricted to my subbox
    ins_face->generateRandomOwnedGhost(pos);

    /*NL*/// if (screen) printVec3D(screen,"pos bef",pos);

    // extrude the position
    r = -1.*(random->uniform()*extrude_length);
    vectorScalarMult3D(normalvec,r,ext);
    vectorAdd3D(pos,ext,pos);

    /*NL*/// if (screen) printVec3D(screen,"pos aft",pos);
    /*NL*/ //error->all(FLERR,"end");
}

/* ----------------------------------------------------------------------
   generate random positions within extruded face
   perform overlap check via xnear if requested
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertStream::x_v_omega(int ninsert_this_local,int &ninserted_this_local, int &ninserted_spheres_this_local, double &mass_inserted_this_local)
{
    ninserted_this_local = ninserted_spheres_this_local = 0;
    mass_inserted_this_local = 0.;

    int nins;
    double pos[3];
    ParticleToInsert *pti;

    double omega_tmp[] = {0.,0.,0.};

    int ntry = 0;
    int maxtry = ninsert_this_local * maxattempt;

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::x_v_omega() start, maxtry %d\n",maxtry);

    /*NL*/ //if(screen && maxtry > 0) fprintf(screen,"proc %d: ninsert_this_local, maxtry %d sublo/hi %f %f, extrude_length_max %f extrude_length_min %f\n",
    /*NL*/ //              comm->me,ninsert_this_local,maxtry,domain->sublo[2],domain->subhi[2],extrude_length_max,extrude_length_min);

    // no overlap check
    // insert with v_normal, no omega
    if(!check_ol_flag)
    {
        for(int itotal = 0; itotal < ninsert_this_local; itotal++)
        {
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rad_to_insert = pti->r_bound_ins;

            do
            {
                generate_random(pos,rad_to_insert);
                ntry++;
            }
            while(ntry < maxtry && (!domain->is_in_subdomain(pos)));

            //NP only insert if could successfully randomize position
            if(ntry < maxtry)
            {
                // could randomize vel, omega, quat here

                if(quat_random_)
                    MathExtraLiggghts::random_unit_quat(random,quat_insert);

                nins = pti->set_x_v_omega(pos,v_normal,omega_tmp,quat_insert);

                ninserted_spheres_this_local += nins;
                mass_inserted_this_local += pti->mass_ins;
                ninserted_this_local++;
            }
        }
    }
    // overlap check
    // account for maxattempt
    // pti checks against xnear and adds self contributions
    else
    {
        while(ntry < maxtry && ninserted_this_local < ninsert_this_local)
        {
            pti = fix_distribution->pti_list[ninserted_this_local];
            double rad_to_insert = pti->r_bound_ins;

            nins = 0;
            while(nins == 0 && ntry < maxtry)
            {
                do
                {
                    generate_random(pos,rad_to_insert);
                    ntry++;

                }
                while(ntry < maxtry && ((!domain->is_in_subdomain(pos)) || (domain->dist_subbox_borders(pos) < rad_to_insert)));

                //NP only insert if could successfully randomize position
                if(ntry < maxtry)
                {
                    // could randomize vel, omega, quat here

                    if(quat_random_)
                        MathExtraLiggghts::random_unit_quat(random,quat_insert);
                    /*NL*/ //if (screen) fprintf(screen,"quat %f %f %f %f\n",quat_insert[0],quat_insert[1],quat_insert[2],quat_insert[3]);

                    /*NL*///if (screen) fprintf(screen,"domain->dist_subbox_borders(pos) %f rad_to_insert %f\n",domain->dist_subbox_borders(pos),rad_to_insert);
                    nins = pti->check_near_set_x_v_omega(pos,v_normal,omega_tmp,quat_insert,neighList);
                }
            }

            if(nins > 0)
            {
                ninserted_spheres_this_local += nins;
                mass_inserted_this_local += pti->mass_ins;
                ninserted_this_local++;
            }
        }
    }

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::x_v_omega() end\n");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::finalize_insertion(int ninserted_spheres_this_local)
{
    // nins particles have been inserted on this proc, set initial position, insertion step and release step according to pos
    //NP assume last ninserted_spheres_this_local have been inserted by this command

    int n_steps = -1;
    int step = update->ntimestep;
    int ilo = atom->nlocal - ninserted_spheres_this_local;
    int ihi = atom->nlocal;

    double pos_rel[3], dist_normal;
    double **x = atom->x;
    double dt = update->dt;

    double **release_data = fix_release->array_atom;

    MultisphereParallel *multisphere = NULL;
    if(fix_multisphere) multisphere = &fix_multisphere->data();

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::finalize_insertion() start, ilo %d, ihi %d, nlocal %d v_normal %f %f %f\n",ilo,ihi,atom->nlocal,v_normal[0],v_normal[1],v_normal[2]);

    for(int i = ilo; i < ihi; i++)
    {
        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::finalize_insertion() i %d, 0\n",i);

        //NP calc_n_steps sets insertion step and velocity
        if(multisphere)
            n_steps = fix_multisphere->calc_n_steps(i,p_ref,normalvec,v_normal);

        if(!multisphere || n_steps == -1)
        {
            vectorSubtract3D(p_ref,x[i],pos_rel);
            dist_normal = vectorDot3D(pos_rel,normalvec);
            n_steps = static_cast<int>((dist_normal+FIX_INSERT_STREAM_TINY)/(vectorMag3D(v_normal)*dt));
        }

        // first 3 values is original position to integrate
        vectorCopy3D(x[i],release_data[i]);

        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::finalize_insertion() i %d, 1\n",i);

        // 4th value is insertion step
        release_data[i][3] = static_cast<double>(step);

        // 5th value is step to release
        release_data[i][4] = static_cast<double>(step + n_steps);

        // 6-8th value is integration velocity
        vectorCopy3D(v_normal,&release_data[i][5]);

        // set inital conditions
        // randomize vel, omega, quat here
        double v_toInsert[3],omega_toInsert[3];

        vectorCopy3D(v_insert,v_toInsert);
        vectorCopy3D(omega_insert,omega_toInsert);

        //NP uniform gaussian work only for single spheres
        //NP for multisphere, would have to do this outside here
        //NP since release step of single particles is different
        //NP from release step for body

        // could ramdonize vel, omega, quat here
        generate_random_velocity(v_toInsert);

        // 9-11th value is velocity, 12-14 is omega
        vectorCopy3D(v_toInsert,&release_data[i][8]);
        vectorCopy3D(omega_toInsert,&release_data[i][11]);

        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::finalize_insertion() i %d, 2\n",i);
    }

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::finalize_insertion() end\n");

    //NP callback for stream tracers
    if(ntracer)
        for(int i = 0; i < ntracer; i++)
            tracer[i]->mark_tracers(ilo,ihi);
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::end_of_step()
{
    int r_step, i_step;

    int step = update->ntimestep;
    int nlocal = atom->nlocal;
    double **release_data = fix_release->array_atom;
    double time_elapsed, dist_elapsed[3], v_integrate[3], *v_toInsert, *omega_toInsert;
    double dt = update->dt;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int *mask = atom->mask;

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::end_of_step() start\n");

    // only one fix handles the integration
    if(!i_am_integrator) return;

    for(int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if(release_data[i][3] == 0.) continue;

            i_step = static_cast<int>(release_data[i][3]+FIX_INSERT_STREAM_TINY);
            r_step = static_cast<int>(release_data[i][4]+FIX_INSERT_STREAM_TINY);
            vectorCopy3D(&release_data[i][5],v_integrate);

            if(step > r_step)
            {
                continue;
            }
            else if (r_step == step)
            {
                //NP dont do this for multisphere, skip to next i in for loop
                if(fix_multisphere && fix_multisphere->belongs_to(i) >= 0)
                {
                    //NP v, omage stored in all particles
                    v_toInsert = &release_data[i][8];
                    omega_toInsert = &release_data[i][11];
                    fix_multisphere->release(i,v_toInsert,omega_toInsert);
                    continue;
                }

                // integrate with constant vel and set v,omega

                time_elapsed = (step - i_step) * dt;

                // particle moves with v_integrate
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                double *x_ins = release_data[i];

                // set x,v,omega
                // zero out force, torque

                vectorAdd3D(x_ins,dist_elapsed,x[i]);

                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);

                v_toInsert = &release_data[i][8];
                omega_toInsert = &release_data[i][11];

                vectorCopy3D(v_toInsert,v[i]);
                vectorCopy3D(omega_toInsert,omega[i]);

                /*NL*///if(fix_rm) error->one(FLERR,"must set params for frm");
            }
            // step < r_step, only true for inserted particles
            //   b/c r_step is 0 for all other particles
            // integrate with constant vel
            else
            {
                //NP do for multisphere particles as well
                //NP but is overridden by set_xv()

                time_elapsed = (step - i_step) * dt;

                // particle moves with v_integrate
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                double *x_ins = release_data[i];

                // set x,v,omega
                vectorAdd3D(x_ins,dist_elapsed,x[i]);
                vectorCopy3D(v_integrate,v[i]);
                vectorZeroize3D(omega[i]);

                // zero out force, torque
                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);
            }
        }
    }

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStream::end_of_step() end\n");
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::reset_timestep(bigint newstep,bigint oldstep)
{
    //NP have to reset release data in case time step
    //NP is changed
    reset_releasedata(newstep,oldstep);
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::reset_releasedata(bigint newstep,bigint oldstep)
{
  //NP need to reset releasedata in case of restart, since
  //NP reset_timestep might have been called

  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **release_data = fix_release->array_atom;

  for(int i = 0; i < nlocal; i++)
  {
        /*NL*/ //if (screen) fprintf(screen,"release_data[i][4]-release_data[i][3] %f\n",release_data[i][4]-release_data[i][3]);

        // first 3 values is original position to integrate
        vectorCopy3D(x[i],release_data[i]);

        // 4th value is insertion step
        release_data[i][3] -= static_cast<double>(oldstep-newstep);

        // 5th value is step to release
        release_data[i][4] -= static_cast<double>(oldstep-newstep);

        // 6-8th value is integration velocity
        vectorCopy3D(v_normal,&release_data[i][5]);
  }
}
