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

#include "particleToInsert_multisphere.h"
#include "math.h"
#include "error.h"
#include "vector_liggghts.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "domain.h"
#include "comm.h"
#include "fix.h"
#include "fix_multisphere.h"
#include "string.h"
#include "math_extra_liggghts.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ParticleToInsertMultisphere::ParticleToInsertMultisphere(LAMMPS* lmp,int ns) : ParticleToInsert(lmp,ns)
{
    memory->create(displace,nspheres,3,"displace");

    for(int i = 0; i < nspheres; i++)
       vectorZeroize3D(displace[i]);
}

/* ---------------------------------------------------------------------- */

ParticleToInsertMultisphere::~ParticleToInsertMultisphere()
{
    memory->destroy(displace);
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertMultisphere::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{
    double disp_glob[3];

    vectorCopy3D(x,xcm_ins);
    vectorCopy4D(quat,quat_ins);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    //if(!isUnitQuat4D(quat_ins)) error->warning(FLERR,"quaternion rotation untested in ParticleToInsertMultisphere");

    /*NP test quaternion rotation
    double test[] = {1.,2.,3.},testr[3];
    double quatt[] = {1.,0.,0.,0.};
    MathExtraLiggghts::vec_quat_rotate(test,quatt,testr);
    printVec3D(screen,"original vector",test);
    printVec3D(screen,"rotated vector",testr);
    */

    MathExtraLiggghts::vec_quat_rotate(ex_space,quat);
    MathExtraLiggghts::vec_quat_rotate(ey_space,quat);
    MathExtraLiggghts::vec_quat_rotate(ez_space,quat);

    for(int j = 0; j < nspheres; j++)
    {
        MathExtraLiggghts::local_coosys_to_cartesian(disp_glob,displace[j],ex_space,ey_space,ez_space);
        vectorAdd3D(x,disp_glob,x_ins[j]);
    }

    /*NL*///printVec3D(screen,"xcm",xcm_ins);

    return nspheres;
}

/* ---------------------------------------------------------------------- */

//NP checks against xnear list
//NP returns # inerted spheres
//NP if >1, increase nbodies

int ParticleToInsertMultisphere::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{

    // check every sphere against all others in xnear
    // if no overlap add to xnear
    double del[3],disp_glob[3], rsq, radsum;
    double ex_space_try[3], ey_space_try[3], ez_space_try[3];

    // rotate if needed


    // calculate x_ins for this quaternion
    // do this in a "try" step since we do not know if we will succeed

    //if(!isUnitQuat4D(quat_ins)) error->one(FLERR,"quaternion rotation untested in ParticleToInsertMultisphere");
    MathExtraLiggghts::vec_quat_rotate(ex_space,quat,ex_space_try);
    MathExtraLiggghts::vec_quat_rotate(ey_space,quat,ey_space_try);
    MathExtraLiggghts::vec_quat_rotate(ez_space,quat,ez_space_try);
    for(int j = 0; j < nspheres; j++)
    {
        MathExtraLiggghts::local_coosys_to_cartesian(disp_glob,displace[j],ex_space_try,ey_space_try,ez_space_try);
        vectorAdd3D(x,disp_glob,x_ins[j]);
    }

    /*NL*/// int kk = 0;

    for(int i = 0; i < nnear; i++)
    {
        for(int j = 0; j < nspheres; j++)
        {
           vectorSubtract3D(x_ins[j],xnear[i],del);
           rsq = vectorMag3DSquared(del);
           radsum = radius_ins[j] + xnear[i][3];

           /*NL*/// printVec3D(screen,"x_ins[j]",x_ins[j]);
           /*NL*/// fprintf(screen,"rsq %f\n",rsq);
           /*NL*/// fprintf(screen,"radsum %f\n",radsum);
           /*NL*/// kk++; if (kk > 100) error->all(FLERR,"end");

           // no success in overlap
           if (rsq <= radsum*radsum) return 0;
        }
    }

    // no overlap with any other - success

    vectorCopy3D(x,xcm_ins);
    vectorCopy4D(quat,quat_ins);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // set axis to the value we succeeded for, x_ins already up to date
    vectorCopy3D(ex_space_try,ex_space);
    vectorCopy3D(ey_space_try,ey_space);
    vectorCopy3D(ez_space_try,ez_space);

    // add to xnear
    for(int j = 0; j < nspheres; j++)
    {
        vectorCopy3D(x_ins[j],xnear[nnear]);
        xnear[nnear][3] = radius_ins[j];
        nnear++;
    }

    return nspheres;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsertMultisphere::insert()
{
    int inserted = 0;
    int nfix = modify->nfix;
    Fix **fix = modify->fix;

    // perform the actual particle insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups
    //NP subdomain check on xcm is ok for single particle insertion
    //NP    since comm->exchange called right after
    //NP subdomain check on xcm is necessary for rigid body insertion
    //NP   so all particles of one body are inserted on the same proc
    //NP   for this reasom, this makes it possible for fix rigid to init all atoms correctly

    for(int i = 0; i < nspheres; i++)
    {
        /*NL*///fprintf(screen,"proc %d tyring to insert particle at position %f %f %f\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
        //NP do not need subdomain check any longer since have processor-local lists anyway
        //if (domain->is_in_extended_subdomain(xcm_ins))
        //{
            inserted++;
            atom->avec->create_atom(atom_type,x_ins[i]);
            /*NL*///fprintf(screen,"proc %d inserting particle at position %f %f %f, mass %e\n",comm->me,x_ins[i][0],x_ins[i][1],x_ins[i][2],mass_ins);
            int m = atom->nlocal - 1;
            atom->mask[m] = 1 | groupbit;
            atom->radius[m] = radius_ins[i];
            atom->density[m] = density_ins;

            //NP for interaction, the total mass of the template is seen
            atom->rmass[m] = mass_ins; //NP 4.*M_PI/3.*radius_ins[i]*radius_ins[i]*radius_ins[i]*density_ins;

            //NP v, omega set 0 here, calculated from rigid
            //NP body constraint afterwards
            vectorZeroize3D(atom->v[m]);
            vectorZeroize3D(atom->omega[m]);
            vectorZeroize3D(atom->f[m]);
            vectorZeroize3D(atom->torque[m]);

            for (int j = 0; j < nfix; j++)
            {
               if (fix[j]->create_attribute) fix[j]->set_arrays(m);
            }
        //}
    }

    // now rigid body insertion
    //NP as mentioned above the proc owning xcm has inserted all particles
    //NP called only on proc owning xcm
    int init_flag = 0;
    int nlocal = atom->nlocal;

    if(modify->n_fixes_style("multisphere") != 1)
        error->one(FLERR,"Multi-sphere particle inserted: You have to use exactly one fix multisphere");

    FixMultisphere *fix_multisphere = static_cast<FixMultisphere*>(modify->find_fix_style("multisphere",0));

    /*NL*///fprintf(screen,"inserting body at xcm position %f %f %f\n",xcm_ins[0],xcm_ins[1],xcm_ins[2]);
    fix_multisphere->data().add_body(nspheres,xcm_ins,xcm_to_xbound,r_bound_ins, v_ins, omega_ins, mass_ins,
                                density_ins,atom_type,type_ms,inertia,ex_space,ey_space,ez_space,displace);

    // set displace correctly, set body to -2
    //NP -2 means belongs to body, but does not know yet to which
    //NP must assume last nspheres on this proc belong to this body
    int i = 0;
    for(int isphere = nlocal-nspheres; isphere < nlocal; isphere++)
        fix_multisphere->set_body_displace(isphere,displace[i++],-2);

    return inserted;
}

/* ---------------------------------------------------------------------- */

void ParticleToInsertMultisphere::scale_pti(double r_scale)
{
    error->one(FLERR,"scale_pti not implemented in ParticleToInsertMultisphere");
}

/* ---------------------------------------------------------------------- */

void ParticleToInsertMultisphere::random_rotate(double rn1,double rn2, double rn3)
{
    /*NL*/ //fprintf(screen,"will do rotate\n");
    /*NL*/ //for(int i=0;i<nspheres;i++) fprintf(screen,"particle %d: x %e %e %e\n",i,x_ins[i][0],x_ins[i][1],x_ins[i][2]);

    //NP only do something if I am a multisphere particle
    if(nspheres==1)return;

    double *vert_before_rot;
    double vert_after_rot[3];

    //NP 3 random angles
    double phix=rn1*2.*M_PI;
    double phiy=rn2*2.*M_PI;
    double phiz=rn3*2.*M_PI;

    double cos_phix = cos(phix);
    double cos_phiy = cos(phiy);
    double cos_phiz = cos(phiz);
    double sin_phix = sin(phix);
    double sin_phiy = sin(phiy);
    double sin_phiz = sin(phiz);

    //NP rotate ex_space, ey_space, ez_space around the 3 angles
    for(int i=0;i<3;i++)
    {
        if     (i==0) vert_before_rot=ex_space;
        else if(i==1) vert_before_rot=ey_space;
        else if(i==2) vert_before_rot=ez_space;

        vert_after_rot[0] = vert_before_rot[0]*cos_phiy*cos_phiz+vert_before_rot[1]*(cos_phiz*sin_phix*sin_phiy-cos_phix*sin_phiz)+vert_before_rot[2]*(cos_phix*cos_phiz*sin_phiy+sin_phix*sin_phiz);
        vert_after_rot[1] = vert_before_rot[0]*cos_phiy*sin_phiz+vert_before_rot[2]*(-cos_phiz*sin_phix+cos_phix*sin_phiy*sin_phiz)+vert_before_rot[1]*(cos_phix*cos_phiz+sin_phix*sin_phiy*sin_phiz);
        vert_after_rot[2] = vert_before_rot[2]*cos_phix*cos_phiy+vert_before_rot[1]*cos_phiy*sin_phix-vert_before_rot[0]*sin_phiy;

        /*NL*/ //vert_after_rot[0] = vert_before_rot[0]*cos_phiy*cos_phiy+vert_before_rot[1]*(cos_phiy*sin_phix*sin_phiy-cos_phix*sin_phiy)+vert_before_rot[2]*(cos_phix*cos_phiy*sin_phiy+sin_phix*sin_phiy);
        /*NL*/ //vert_after_rot[1] = vert_before_rot[0]*cos_phiy*sin_phiy+vert_before_rot[2]*(-cos_phiy*sin_phix+cos_phix*sin_phiy*sin_phiy)+vert_before_rot[1]*(cos_phix*cos_phiy+sin_phix*sin_phiy*sin_phiy);
        /*NL*/ //vert_after_rot[2] = vert_before_rot[2]*cos_phix*cos_phiy+vert_before_rot[1]*cos_phiy*sin_phix-vert_before_rot[0]*sin_phiy;

        if     (i==0) for(int j=0;j<3;j++) ex_space[j]=vert_after_rot[j];
        else if(i==1) for(int j=0;j<3;j++) ey_space[j]=vert_after_rot[j];
        else if(i==2) for(int j=0;j<3;j++) ez_space[j]=vert_after_rot[j];
    }

    for(int i=0;i<nspheres;i++)
    {
        x_ins[i][0] = xcm_ins[0] + ex_space[0]*displace[i][0] +   ey_space[0]*displace[i][1] +   ez_space[0]*displace[i][2];
        x_ins[i][1] = xcm_ins[1] + ex_space[1]*displace[i][0] +   ey_space[1]*displace[i][1] +   ez_space[1]*displace[i][2];
        x_ins[i][2] = xcm_ins[2] + ex_space[2]*displace[i][0] +   ey_space[2]*displace[i][1] +   ez_space[2]*displace[i][2];
    }

    /*NL*/ //fprintf(screen,"did rotate\n");
    /*NL*/ //for(int i=0;i<nspheres;i++) fprintf(screen,"particle %d: x %e %e %e\n",i,x_ins[i][0],x_ins[i][1],x_ins[i][2]);
}
