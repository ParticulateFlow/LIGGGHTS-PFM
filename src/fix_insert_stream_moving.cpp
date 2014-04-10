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
#include "fix_insert_stream_moving.h"
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
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "tri_mesh_planar.h"

enum{FACE_NONE,FACE_MESH,FACE_CIRCLE};

using namespace LAMMPS_NS;
using namespace FixConst;

/*NL*/ #define LMP_DEBUGMODE_FIXINSERT_STREAM false
/*NL*/ #define LMP_DEBUG_OUT_FIXINSERT_STREAM screen

#define FIX_INSERT_NTRY_SUBBOX 500
#define FIX_INSERT_STREAM_TINY 1e-14

/* ---------------------------------------------------------------------- */

FixInsertStreamMoving::FixInsertStreamMoving(LAMMPS *lmp, int narg, char **arg) :
  FixInsertStream(lmp, narg, arg)
{
    ins_face_planar = 0;
}

/* ---------------------------------------------------------------------- */

FixInsertStreamMoving::~FixInsertStreamMoving()
{
}

/* ---------------------------------------------------------------------- */

void FixInsertStreamMoving::post_create()
{
  FixInsert::post_create();

  // only register property if I am the first fix/insert/stream in the simulation
  //NP 16 values: triangle ID (1), bary coords (3), extrude length (1)
  //NP           start step (1), release step (1), integration velocity (3), vel(3), omega(3)

  if(modify->n_fixes_style(style) == 1)
  {
        const char* fixarg[24];
        fixarg[0]="release_fix_insert_stream";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="release_fix_insert_stream";
        fixarg[4]="vector"; //NP 1 scalar per particle to be registered
        fixarg[5]="yes";    //NP restart yes
        fixarg[6]="yes";    //NP communicate ghost no
        fixarg[7]="no";    //NP communicate rev yes
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fixarg[11]="0.";
        fixarg[12]="0.";
        fixarg[13]="0.";
        fixarg[14]="0.";
        fixarg[15]="0.";
        fixarg[16]="0.";
        fixarg[17]="0.";
        fixarg[18]="0.";
        fixarg[19]="0.";
        fixarg[20]="0.";
        fixarg[21]="0.";
        fixarg[22]="0.";
        fixarg[23]="0.";
        modify->add_fix_property_atom(24,const_cast<char**>(fixarg),style);
  }
}

/* ---------------------------------------------------------------------- */

void FixInsertStreamMoving::init()
{
    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::init() start\n");

    FixInsert::init();

    fix_release = static_cast<FixPropertyAtom*>(modify->find_fix_property("release_fix_insert_stream","property/atom","vector",10,0,style));
    if(!fix_release) error->fix_error(FLERR,this,"Internal error if fix insert/stream");

    i_am_integrator = modify->i_am_first_of_style(this);

    //NP currently fix rigid disallowed
    if(fix_multisphere)
        error->fix_error(FLERR,this,"cannot use with multi-sphere");

    //NP disallow scaling because it would undermine
    if(ins_face->isScaling())
        error->fix_error(FLERR,this,"cannot scale mesh which is used for particle insertion");

    //NP have to use a planar face to be able to search coords
    ins_face_planar = dynamic_cast<TriMeshPlanar*>(ins_face);
    if(!ins_face_planar)
        error->fix_error(FLERR,this,"requires you to use a fix mesh/surface/planar");
}

/* ---------------------------------------------------------------------- */

void FixInsertStreamMoving::finalize_insertion(int ninserted_spheres_this_local)
{
    // nins particles have been inserted on this proc, set initial position, insertion step and release step according to pos

    int n_steps = -1;
    int step = update->ntimestep;
    int ilo = atom->nlocal - ninserted_spheres_this_local;
    int ihi = atom->nlocal;

    int tri_id;
    double bary[3], dist_normal;
    double **x = atom->x;
    double dt = update->dt;

    double **release_data = fix_release->array_atom;

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::finalize_insertion() start, ilo %d, ihi %d, nlocal %d\n",ilo,ihi,atom->nlocal);

    for(int i = ilo; i < ihi; i++)
    {
        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::finalize_insertion() i %d, 0\n",i);

        // calculate relevant information
        ins_face_planar->locatePosition(x[i],tri_id,bary,dist_normal);

        /*NL*/ //double xtest[3];
        /*NL*/ //ins_face_planar->constructPositionFromBary(tri_id,bary,xtest);
        /*NL*/ //fprintf(screen,"dist_normal %f\n",dist_normal);
        /*NL*/ //printVec3D(screen,"x[i]",x[i]);
        /*NL*/ //printVec3D(screen,"bary",bary);
        /*NL*/ //printVec3D(screen,"xtest",xtest);

        n_steps = static_cast<int>((dist_normal+FIX_INSERT_STREAM_TINY)/(vectorMag3D(v_normal)*dt));

        // 0 = tri ID
        release_data[i][0] = static_cast<double>(tri_id);

        // 1,2,3 = bary coos
        vectorCopy3D(bary,&(release_data[i][1]));

        // 4 = normal dist
        release_data[i][4] = dist_normal;

        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::finalize_insertion() i %d, 1\n",i);

        // 5 = insertion step
        release_data[i][5] = static_cast<double>(step);

        // 6 = step to release
        release_data[i][6] = static_cast<double>(step + n_steps);

        // 7,8,9 = integration velocity
        vectorCopy3D(v_normal,&release_data[i][7]);

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
        if(v_randomSetting==1)
        {
            v_toInsert[0] = v_insert[0] + v_insertFluct[0] * 2.0 * (random->uniform()-0.50);
            v_toInsert[1] = v_insert[1] + v_insertFluct[1] * 2.0 * (random->uniform()-0.50);
            v_toInsert[2] = v_insert[2] + v_insertFluct[2] * 2.0 * (random->uniform()-0.50);
        }
        else if(v_randomSetting==2)
        {
            v_toInsert[0] = v_insert[0] + v_insertFluct[0] * random->gaussian();
            v_toInsert[1] = v_insert[1] + v_insertFluct[1] * random->gaussian();
            v_toInsert[2] = v_insert[2] + v_insertFluct[2] * random->gaussian();
        }

        // 10 11 12 is velocity, 13 14 15 is omega
        vectorCopy3D(v_toInsert,&release_data[i][10]);
        vectorCopy3D(omega_toInsert,&release_data[i][13]);

        /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::finalize_insertion() i %d, 2\n",i);
    }

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::finalize_insertion() end\n");
}

/* ---------------------------------------------------------------------- */

void FixInsertStreamMoving::end_of_step()
{
    int r_step, i_step, tri_id;

    int step = update->ntimestep;
    int nlocal = atom->nlocal;
    double **release_data = fix_release->array_atom;
    double time_elapsed, dist_elapsed[3], v_integrate[3], *v_toInsert, *omega_toInsert;
    double x_ins[3], bary[3], dist_extrude0[3], dist_extrude[3];
    double dt = update->dt;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int *mask = atom->mask;

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::end_of_step() start\n");

    // only one fix handles the integration
    if(!i_am_integrator) return;

    for(int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            if(release_data[i][5] == 0.) continue;

            i_step = static_cast<int>(release_data[i][5]+FIX_INSERT_STREAM_TINY);
            r_step = static_cast<int>(release_data[i][6]+FIX_INSERT_STREAM_TINY);
            vectorCopy3D(&release_data[i][7],v_integrate);

            if(step > r_step) continue;
            else if (r_step == step)
            {
                // integrate with constant vel and set v,omega

                time_elapsed = (step - i_step) * dt;

                // construct position
                tri_id = static_cast<int>(release_data[i][0]);
                vectorCopy3D(&(release_data[i][1]),bary);
                if(!ins_face_planar->constructPositionFromBary(tri_id,bary,x_ins))
                    error->one(FLERR,"Internal error");

                // particle moves with v_integrate
                vectorScalarMult3D(normalvec,-release_data[i][4],dist_extrude0);
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                vectorAdd3D(dist_extrude0,dist_elapsed,dist_extrude);

                // set x,v,omega
                // zero out force, torque

                vectorAdd3D(x_ins,dist_extrude,x[i]);

                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);

                v_toInsert = &release_data[i][10];
                omega_toInsert = &release_data[i][13];

                vectorCopy3D(v_toInsert,v[i]);
                vectorCopy3D(omega_toInsert,omega[i]);
            }
            // step < r_step, only true for inserted particles
            //   b/c r_step is 0 for all other particles
            // integrate with constant vel
            else
            {
                time_elapsed = (step - i_step) * dt;

                // construct position
                tri_id = static_cast<int>(release_data[i][0]);
                vectorCopy3D(&(release_data[i][1]),bary);
                if(!ins_face_planar->constructPositionFromBary(tri_id,bary,x_ins))
                    error->one(FLERR,"Internal error");

                // particle moves with v_integrate
                vectorScalarMult3D(normalvec,-release_data[i][4],dist_extrude0);
                vectorScalarMult3D(v_integrate,time_elapsed,dist_elapsed);
                vectorAdd3D(dist_extrude0,dist_elapsed,dist_extrude);

                // set x,v,omega
                vectorAdd3D(x_ins,dist_extrude,x[i]);
                vectorCopy3D(v_integrate,v[i]);
                vectorZeroize3D(omega[i]);

                // zero out force, torque
                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);
            }
        }
    }

    /*NL*/ if(LMP_DEBUGMODE_FIXINSERT_STREAM) fprintf(LMP_DEBUG_OUT_FIXINSERT_STREAM,"FixInsertStreamMoving::end_of_step() end\n");
}

/* ---------------------------------------------------------------------- */

void FixInsertStreamMoving::reset_releasedata(bigint newstep,bigint oldstep)
{
  //NP need to reset releasedata in case of restart, since
  //NP reset_timeste might have been called

  int nlocal = atom->nlocal;
  double **x = atom->x;
  double **release_data = fix_release->array_atom;
  double bary[3], dist_normal;
  int tri_id;

  for(int i = 0; i < nlocal; i++)
  {
        // calculate relevant information
        ins_face_planar->locatePosition(x[i],tri_id,bary,dist_normal);

        // 0 = tri ID
        release_data[i][0] = static_cast<double>(tri_id);

        // 1,2,3 = bary coos
        vectorCopy3D(bary,&(release_data[i][1]));

        // 4 = normal dist
        release_data[i][4] = dist_normal;

        // 5 = insertion step
        release_data[i][5] -= static_cast<double>(oldstep-newstep);

        // 6 = step to release
        release_data[i][6] -= static_cast<double>(oldstep-newstep);

        // 7,8,9 = integration velocity
        vectorCopy3D(v_normal,&release_data[i][7]);
  }
}
