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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix_multisphere.h"
#include "domain_wedge.h"
#include "math_extra.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "respa.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include "force.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_template_multisphere.h"
#include "neighbor.h"
#include "fix_gravity.h"
#include "vector_liggghts.h"
#include "mpi_liggghts.h"
#include "atom_vec.h"
#include "math_extra_liggghts.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define TOLERANCE 1.0e-6
#define EPSILON 1.0e-7
#define MAXJACOBI 50
#define DELTA_GROW 10000

/*NL*/ #define LMP_DEBUGMODE_RIGID_MS false //(update->ntimestep>845 && comm->me==0)
/*NL*/ #define LMP_DEBUGMODE_RIGID_MS_OUT screen //(update->ntimestep>845 && comm->me==0)

enum {LOOP_LOCAL,LOOP_ALL};

/* ---------------------------------------------------------------------- */

FixMultisphere::FixMultisphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  do_modify_body_forces_torques_(false),
  multisphere_(*(new MultisphereParallel(lmp))),
  fix_corner_ghost_(0),
  fix_delflag_(0),
  fix_existflag_(0),
  fix_gravity_(0),
  fw_comm_flag_(MS_COMM_UNDEFINED),
  rev_comm_flag_(MS_COMM_UNDEFINED),
  body_(NULL),
  displace_(NULL),
  ntypes_(0),
  Vclump_(0)
{

  //NP check if style is molecular
  if(atom->molecular == 1)
    error->fix_error(FLERR,this,"Must NOT use a hybrid sphere/molecular atom style with fix multisphere (use sphere only)");

  atom->molecule_flag = 1;
  grow_arrays(atom->nmax);

  //NP make an exclusion by molecule id;only particles in the fix rigid group are included
  char **modarg;
  modarg = new char*[3];
  modarg[2] = new char[50];
  modarg[0] = (char*) "exclude";
  modarg[1] = (char*) "molecule";
  strcpy(modarg[2],arg[1]);
  neighbor->modify_params(3,modarg);
  delete [] modarg[2];
  delete []modarg;

  //NP per atom creation and restart
  restart_global = 1;
  restart_peratom = 1;
  restart_pbc = 1;
  atom->add_callback(0);
  atom->add_callback(1);

  // fix handles properties that need to be initialized at particle creation
  create_attribute = 1;

  //NP this fix can force reneighboring
  force_reneighbor = 1;
  next_reneighbor = -1;

  //NP modified C.K.
  // is now local data, not global
  local_flag = 1;

  size_local_rows = 0;    //NP init with 0 particles
  size_local_cols = 12;           // 0 = vector, N = columns in local array
  local_freq = 1;

  size_peratom_cols = 0;

  vector_flag = 1;
  size_vector = 0; // no bodies present at creation

  global_freq = 1;
  extarray = 0;

  //NP max # comm for case f, torque, flag
  comm_forward = 7;

  //NP max # comm for case pos, vel, omega, flag
  comm_reverse = 10;
}

/* ---------------------------------------------------------------------- */

FixMultisphere::~FixMultisphere()
{
    atom->delete_callback(id,0);
    atom->delete_callback(id,1);

    delete &multisphere_;

    memory->destroy(displace_);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::post_create()
{
    //NP register corner ghost
    if(!fix_corner_ghost_)
    {
        const char* fixarg[9];
        fixarg[0]="cornerghost";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="cornerghost";
        fixarg[4]="scalar";
        fixarg[5]="no";    // restart
        fixarg[6]="no";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_corner_ghost_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    //NP register deletion flag
    if(!fix_delflag_)
    {
        const char* fixarg[9];
        fixarg[0]="delflag";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="delflag";
        fixarg[4]="scalar";
        fixarg[5]="yes";     // restart
        fixarg[6]="no";      // communicate ghost
        fixarg[7]="yes";     // communicate rev
        fixarg[8]="0.";
        fix_delflag_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }
    //NP register deletion flag
    if(!fix_existflag_)
    {
        const char* fixarg[9];
        fixarg[0]="existflag";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="existflag";
        fixarg[4]="scalar";
        fixarg[5]="no";     // restart
        fixarg[6]="no";      // communicate ghost
        fixarg[7]="yes";     // communicate rev
        fixarg[8]="1.";
        fix_existflag_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
    }

    //NP in case of restart: see comment in FixMultisphere::restart
    if(modify->have_restart_data(this))
    {
        evflag = 0;
        set_xv(LOOP_LOCAL);
    }
}

/* ---------------------------------------------------------------------- */

int FixMultisphere::setmask()
{
    int mask = 0;
    mask |= INITIAL_INTEGRATE;
    mask |= PRE_EXCHANGE;
    mask |= PRE_NEIGHBOR;
    mask |= FINAL_INTEGRATE;
    return mask;
}

/* ---------------------------------------------------------------------- */

double FixMultisphere::max_r_bound() const
{
    return multisphere_.max_r_bound();
}

/* ---------------------------------------------------------------------- */

double FixMultisphere::extend_cut_ghost() const
{
    //NP this is to extend ghost region
    //NP one rbound is enough since just need to ensure proc has
    //NP ghost atoms for all owned bodies
    return /*NP*2.*/max_r_bound();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::add_body_finalize()
{
    //NP need to do tasks that are done in FixRigid::setup()
    //NP this function called via particle templates

    //NP invoke set_v()
    //NP since not set for particles belonging to newly inserted bodies
    //NP but do NOT loop ghosts at this stage, new particles do not have ghosts at this point


    /*NL*/ if(LMP_DEBUGMODE_RIGID_MS) fprintf(LMP_DEBUGMODE_RIGID_MS_OUT,"FixMultisphere::add_body_finalize doing setup\n");

    multisphere_.id_extend_body_extend(body_);
    multisphere_.generate_map();
    multisphere_.reset_forces(true);
    set_xv(LOOP_LOCAL); //NP only loop local particles, ghosts do not exist yet for these new particles
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::init() //NP modified C.K.
{
  // lots of error checks and warnings

  if(0 == atom->map_style)
      error->fix_error(FLERR,this,"requires an 'atom_modify map' command to allocate an atom map");

  if(!atom->rmass_flag || !atom->omega_flag)
    error->fix_error(FLERR,this,"need per-atom mass and omega");

  if(domain->dimension != 3)
    error->fix_error(FLERR,this,"works with 3D simulations only");

  if(modify->n_fixes_style("heat/gran") > 0)
    error->fix_error(FLERR,this,"is not compatible with heat transfer simulations");

  if(domain->triclinic || dynamic_cast<DomainWedge*>(domain))
    error->fix_error(FLERR,this,"does not work with triclinic or wedge box");

  if (strstr(update->integrate_style,"respa"))
    error->fix_error(FLERR,this,"does not work with respa");
    //step_respa = ((Respa *) update->integrate)->step;

  //NP this is because reverse comm would maybe erase force and torque values??
  if(force->newton) error->fix_error(FLERR,this,"requires newton 'off'");

  //NP check if a fix gravity is registered
  if(modify->n_fixes_style("gravity") > 1)
    error->fix_error(FLERR,this,"only one fix gravity supported");
  fix_gravity_ = static_cast<FixGravity*>(modify->find_fix_style("gravity",0));

  // warn if more than one rigid fix
  if(modify->n_fixes_style("rigid") + modify->n_fixes_style("multisphere") > 1)
    error->warning(FLERR,"More than one fix rigid / fix multisphere");

  fix_remove_.clear();

  // timestep info

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dtq = 0.5 * update->dt;


}

/* ---------------------------------------------------------------------- */

void FixMultisphere::add_remove_callback(FixRemove *ptr)
{
    fix_remove_.push_back(ptr);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup(int vflag)
{
  // calc MS comm properties
  ntypes_ = modify->n_fixes_style("particletemplate/multisphere");

  if(Vclump_) delete []Vclump_;
  Vclump_ = new double [ntypes_+1];

  for(int ifix = 0; ifix < ntypes_; ifix++)
  {
      FixTemplateMultisphere *ftm =  static_cast<FixTemplateMultisphere*>(modify->find_fix_style("particletemplate/multisphere",ifix));
      int itype = ftm->type();
      Vclump_[itype] = ftm->volexpect();
  }

  //NP guesstimate virial as 2x the set_v contribution
  int i,n;
  int nlocal = atom->nlocal;

  // virial setup before call to set_v

  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (vflag_global)
    for (n = 0; n < 6; n++) virial[n] *= 2.0;
  if (vflag_atom) {
    for (i = 0; i < nlocal; i++)
      for (n = 0; n < 6; n++)
        vatom[i][n] *= 2.0;
  }

  //NP execute communication routine
  calc_force();

  /*NL*/ if(LMP_DEBUGMODE_RIGID_MS && screen) fprintf(screen,"FixMultisphere::setup finished on proc %d\n",comm->me);
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup_pre_exchange()
{
    pre_exchange();
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::setup_pre_neighbor()
{
    /*NL*/ //if (screen) fprintf(screen,"SUPN A step %d proc %d: nbody %d\n",update->ntimestep,comm->me,n_body());
    pre_neighbor();
    /*NL*/ //if (screen) fprintf(screen,"SUPN B step %d proc %d: nbody %d\n",update->ntimestep,comm->me,n_body());
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::set_arrays(int i)
{
    /*NL*/ //if (screen) fprintf(screen,"set arrays called\n");
    body_[i] = -1;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixMultisphere::copy_arrays(int i, int j,int delflag)
{
    body_[j] = body_[i];
    displace_[j][0] = displace_[i][0];
    displace_[j][1] = displace_[i][1];
    displace_[j][2] = displace_[i][2];
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::initial_integrate(int vflag)
{
  double dtfm;
  int timestep = update->ntimestep;
  double **xcm = multisphere_.xcm_.begin();
  double **vcm = multisphere_.vcm_.begin();
  double **fcm = multisphere_.fcm_.begin();
  double **torquecm = multisphere_.torquecm_.begin();
  double **ex_space = multisphere_.ex_space_.begin();
  double **ey_space = multisphere_.ey_space_.begin();
  double **ez_space = multisphere_.ez_space_.begin();
  double **angmom = multisphere_.angmom_.begin();
  double **omega = multisphere_.omega_.begin();
  double **quat = multisphere_.quat_.begin();
  double **inertia = multisphere_.inertia_.begin();
  double *masstotal = multisphere_.masstotal_.begin();
  int *start_step = multisphere_.start_step_.begin();
  double **v_integrate = multisphere_.v_integrate_.begin();
  bool **fflag = multisphere_.fflag_.begin();
  bool **tflag = multisphere_.tflag_.begin();
  int nbody = multisphere_.n_body();

  if(strstr(style,"nointegration"))
    return;

  /*NL*/ //if(screen) fprintf(screen,"nbody_all %d\n",n_body_all());
  /*NL*/ //if(screen && map(7833) >= 0) fprintf(screen,"proc %d has body %d at step %d\n",comm->me,7833,update->ntimestep);
  bool has_stream = modify->n_fixes_style("insert/stream") > 0;

  for (int ibody = 0; ibody < nbody; ibody++) {

    if(has_stream && timestep < start_step[ibody])
    {
        vectorCopy3D(v_integrate[ibody],vcm[ibody]);

        // update xcm by full step
        xcm[ibody][0] += dtv * vcm[ibody][0];
        xcm[ibody][1] += dtv * vcm[ibody][1];
        xcm[ibody][2] += dtv * vcm[ibody][2];
        /*NL*/ //if (screen) printVec3D(screen,"vcm[ibody] FixMultisphere::initial_integrate",vcm[ibody]);
        /*NL*/ //if (screen) printVec3D(screen,"xcm[ibody] FixMultisphere::initial_integrate",xcm[ibody]);
        /*NL*/ //error->one(FLERR,"end");
        /*NL*/ //if (screen) fprintf(screen,"this should not be\n");
        continue;
    }

    /*NL*/ //if(LMP_DEBUGMODE_RIGID) fprintf(LMP_DEBUG_OUTP_RIGID,"FixRigid::initial_integrate for body %d\n",ibody);

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];

    if(fflag[ibody][0]) vcm[ibody][0] += dtfm * fcm[ibody][0];
    if(fflag[ibody][1]) vcm[ibody][1] += dtfm * fcm[ibody][1];
    if(fflag[ibody][2]) vcm[ibody][2] += dtfm * fcm[ibody][2];

    // update xcm by full step

    xcm[ibody][0] += dtv * vcm[ibody][0];
    xcm[ibody][1] += dtv * vcm[ibody][1];
    xcm[ibody][2] += dtv * vcm[ibody][2];

    /*NL*/ //if (screen) {
    /*NL*/ //fprintf(screen,"masstotal[ibody] %f\n",masstotal[ibody]);
    /*NL*/ //printVec3D(screen,"xcm[ibody]",xcm[ibody]);
    /*NL*/ //printVec3D(screen,"vcm[ibody]",vcm[ibody]);
    /*NL*/ //printVec3D(screen,"fcm[ibody]",fcm[ibody]);
    /*NL*/ //}

    // update angular momentum by 1/2 step

    if(tflag[ibody][0]) angmom[ibody][0] += dtf * torquecm[ibody][0];
    if(tflag[ibody][1]) angmom[ibody][1] += dtf * torquecm[ibody][1];
    if(tflag[ibody][2]) angmom[ibody][2] += dtf * torquecm[ibody][2];

    // compute omega at 1/2 step from angmom at 1/2 step and current q
    // update quaternion a full step via Richardson iteration
    // returns new normalized quaternion, also updated omega at 1/2 step
    // update ex,ey,ez to reflect new quaternion

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
    MathExtra::richardson(quat[ibody],angmom[ibody],omega[ibody],
                          inertia[ibody],dtq);
    MathExtra::q_to_exyz(quat[ibody],
                         ex_space[ibody],ey_space[ibody],ez_space[ibody]);

    /*NL*/ //bool eval =  true; //13500 < update->ntimestep && 14000 > update->ntimestep;

    /*NL*/ //if(screen && tag(ibody) == 88 && eval) {
    /*NL*/ //     fprintf(screen,"step " BIGINT_FORMAT " proc %d, xcm %f %f %f vcm %f %f %f omegacm %f %f %f\n",
    /*NL*/ //                                      update->ntimestep,comm->me,
    /*NL*/ //                                     xcm[ibody][0],xcm[ibody][1],xcm[ibody][2],
    /*NL*/ //                                      vcm[ibody][0],vcm[ibody][1],vcm[ibody][2],
    /*NL*/ //                                      data().omega_(ibody)[0],data().omega_(ibody)[1],data().omega_(ibody)[2]);
    /*NL*/ //    fprintf(screen,"step " BIGINT_FORMAT " proc %d, fcm %f %f %f torquecm %f %f %f\n",
    /*NL*/ //                                     update->ntimestep,comm->me,
    /*NL*/ //                                      fcm[ibody][0],fcm[ibody][1],fcm[ibody][2],
    /*NL*/ //                                      torquecm[ibody][0],torquecm[ibody][1],torquecm[ibody][2]);
    /*NL*/ //}
  }

  // virial setup before call to set_xv

  if (vflag) v_setup(vflag);
  else evflag = 0;

  // set coords/orient and velocity/rotation of atoms in rigid bodies
  // from quarternion and omega

  set_xv();

  //NP reverse comm of x,v,omega
  rev_comm_flag_ = MS_COMM_REV_X_V_OMEGA;
  reverse_comm();

  //NP forward comm of x, v, omega done via verlet
}

/* ---------------------------------------------------------------------- */

void FixMultisphere::final_integrate()
{
  double dtfm;
  int timestep = update->ntimestep;
  double **vcm = multisphere_.vcm_.begin();
  double **fcm = multisphere_.fcm_.begin();
  double **torquecm = multisphere_.torquecm_.begin();
  double **ex_space = multisphere_.ex_space_.begin();
  double **ey_space = multisphere_.ey_space_.begin();
  double **ez_space = multisphere_.ez_space_.begin();
  double **angmom = multisphere_.angmom_.begin();
  double **omega = multisphere_.omega_.begin();
  double **inertia = multisphere_.inertia_.begin();
  double *masstotal = multisphere_.masstotal_.begin();
  int *start_step = multisphere_.start_step_.begin();
  bool **fflag = multisphere_.fflag_.begin();
  bool **tflag = multisphere_.tflag_.begin();
  int nbody = multisphere_.n_body();

  //NP sum over atoms to get force and torque on rigid body

  //NP differences to FixRigid::final_integrate() version:
  //NP + most importantly, do not do the allreduce - each proc owns its own bodies
  //NP + do a forward comm and a reverse comm instead
  //NP + loop over owned and ghost atoms
  //NP + check if a fix gravity is registered. If yes, add contribution to rigid body

  // calculate forces and torques on body

  calc_force();

  if(strstr(style,"nointegration"))
    return;

  bool has_stream = modify->n_fixes_style("insert/stream") > 0;

  // resume integration
  for (int ibody = 0; ibody < nbody; ibody++)
  {
    if(has_stream && timestep < start_step[ibody]) continue;

    /*NL*/ //bool eval =  13500 < update->ntimestep && 14000 > update->ntimestep;

    /*NL*/ //if(screen && tag(ibody) == 60 && eval) {
    /*NL*/ //     fprintf(screen,"step " BIGINT_FORMAT " final integrate proc %d, xcm %f %f %f vcm %f %f %f omegacm %f %f %f\n",
    /*NL*/ //                                      update->ntimestep,comm->me,
    /*NL*/ //                                     xcm[ibody][0],xcm[ibody][1],xcm[ibody][2],
    /*NL*/ //                                      vcm[ibody][0],vcm[ibody][1],vcm[ibody][2],
    /*NL*/ //                                      data().omega_(ibody)[0],data().omega_(ibody)[1],data().omega_(ibody)[2]);}

    // update vcm by 1/2 step

    dtfm = dtf / masstotal[ibody];
    if(fflag[ibody][0]) vcm[ibody][0] += dtfm * fcm[ibody][0];
    if(fflag[ibody][1]) vcm[ibody][1] += dtfm * fcm[ibody][1];
    if(fflag[ibody][2]) vcm[ibody][2] += dtfm * fcm[ibody][2];

    // update angular momentum by 1/2 step

    if(tflag[ibody][0]) angmom[ibody][0] += dtf * torquecm[ibody][0];
    if(tflag[ibody][1]) angmom[ibody][1] += dtf * torquecm[ibody][1];
    if(tflag[ibody][2]) angmom[ibody][2] += dtf * torquecm[ibody][2];

    MathExtra::angmom_to_omega(angmom[ibody],ex_space[ibody],ey_space[ibody],
                               ez_space[ibody],inertia[ibody],omega[ibody]);
  }

  //NP set velocity/rotation of atoms in rigid bodies
  //NP virial is already setup from initial_integrate
  set_v();

  //NP reverse comm of v,omega
  rev_comm_flag_ = MS_COMM_REV_V_OMEGA;
  reverse_comm();

  //NP forward comm of v,omega
  fw_comm_flag_ = MS_COMM_FW_V_OMEGA;
  forward_comm();
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

void FixMultisphere::calc_force()
{
  int ibody;
  tagint *image = atom->image;
  double **x = atom->x;
  double **f_atom = atom->f;
  double **torque_atom = atom->torque;
  double f_one[3],torque_one[3];
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  double **xcm = multisphere_.xcm_.begin();
  double *masstotal = multisphere_.masstotal_.begin();
  double **fcm = multisphere_.fcm_.begin();
  double **dragforce_cm = multisphere_.dragforce_cm_.begin();
  double **torquecm = multisphere_.torquecm_.begin();
  int nbody = multisphere_.n_body();

  //NP forward comm of forces and torques
  //NP so have current forces and torques stored in ghosts
  //NP loop not only to nlocal, but also over ghosts
  fw_comm_flag_ = MS_COMM_FW_F_TORQUE;
  forward_comm();

  double unwrap[3],dx,dy,dz;

  // set force and torque to 0
  // do not reset external torques
  multisphere_.reset_forces(false);

  //NP correct forces if necessary
  if(do_modify_body_forces_torques_)
        modify_body_forces_torques();

  // calculate forces and torques of bodies
  for (int i = 0; i < nlocal+nghost; i++)
  {
    //NP skip if atom not in rigid body
    if(body_[i] < 0) continue;

    //NP body ID stored in atom is global
    //NP need to know where stored in my data
    ibody = map(body_[i]);

    //NP skip if body not owned by this proc
    if (ibody < 0) continue;

    //NP skip if periodic ghost of owned particle
    //NP since would double-count forces in this case
    if(!domain->is_owned_or_first_ghost(i))
        continue;

    /*NL*///if (screen) fprintf(screen,"A step " BIGINT_FORMAT ",body tag %d, atom tag %d: force %f %f %f\n",update->ntimestep,tag(ibody),atom->tag[i],fcm[ibody][0],fcm[ibody][1],fcm[ibody][2]);

    vectorCopy3D(f_atom[i],f_one);
    vectorCopy3D(torque_atom[i],torque_one);

    fcm[ibody][0] += f_one[0];
    fcm[ibody][1] += f_one[1];
    fcm[ibody][2] += f_one[2];

    domain->unmap(x[i],image[i],unwrap);
    dx = unwrap[0] - xcm[ibody][0];
    dy = unwrap[1] - xcm[ibody][1];
    dz = unwrap[2] - xcm[ibody][2];

    //NP modified C.K.
    //NP make sure
    //NP this is important for ghost atoms across periodic boundaries
    if(i >= nlocal)
        domain->minimum_image(dx,dy,dz);

    //NP torque due to atom force and torque
    torquecm[ibody][0] += dy*f_one[2] - dz*f_one[1] + torque_one[0];
    torquecm[ibody][1] += dz*f_one[0] - dx*f_one[2] + torque_one[1];
    torquecm[ibody][2] += dx*f_one[1] - dy*f_one[0] + torque_one[2];

    /*NL*/ //bool eval = 3100 < update->ntimestep && 3175 > update->ntimestep;

    /*NL*/ //if(screen && (tag(ibody) == 88) && eval) {
    /*NL*/ //     fprintf(screen,"step " BIGINT_FORMAT " proc %d atom tag %d, ghost %s,pos %f %f %f f %f %f %f torque_one %f %f %f, dx %f dy %f dz %f\n",
    /*NL*/ //                                      update->ntimestep,comm->me,atag[i],i>=nlocal?"yes":"no",
    /*NL*/ //                                      x[i][0],x[i][1],x[i][2],
    /*NL*/ //                                      f_one[0],f_one[1],f_one[2],
    /*NL*/ //                                      torque_one[0],torque_one[1],torque_one[2],dx,dy,dz);
    /*NL*/ //}

    /*NL*///if (screen) fprintf(screen,"B step %d,body tag %d, atom tag %d: force %f %f %f\n",update->ntimestep,tag(ibody),atom->tag[i],fcm[ibody][0],fcm[ibody][1],fcm[ibody][2]);
  }

  /*NL*/ //bool eval = 18000 < update->ntimestep && 20000 > update->ntimestep;
  /*NL*/ //ibody = map(238);
  /*NL*/ //if(screen && ibody >= 0 && eval) {
  /*NL*/ //     fprintf(screen,"step " BIGINT_FORMAT " fcm on body after summation (proc %d) %f %f %f torquecm %f %f %f\n",
  /*NL*/ //                                      update->ntimestep,comm->me,
  /*NL*/ //                                      fcm[ibody][0],fcm[ibody][1],fcm[ibody][2],
  /*NL*/ //                                      torquecm[ibody][0],torquecm[ibody][1],torquecm[ibody][2]);
  /*NL*/ //}

  // add external forces on bodies, such as gravity, dragforce

  if(fix_gravity_)
  {
      double grav[3];
      fix_gravity_->get_gravity(grav);
      for (ibody = 0; ibody < nbody; ibody++)
      {
            fcm[ibody][0] += masstotal[ibody]*grav[0];
            fcm[ibody][1] += masstotal[ibody]*grav[1];
            fcm[ibody][2] += masstotal[ibody]*grav[2];
            /*NL*///if (screen) fprintf(screen,"C step %d,body tag %d: force %f %f %f\n",update->ntimestep,tag(ibody),fcm[ibody][0],fcm[ibody][1],fcm[ibody][2]);
      }
  }

  for (ibody = 0; ibody < nbody; ibody++)
  {
      vectorAdd3D(fcm[ibody],dragforce_cm[ibody],fcm[ibody]);
      /*NL*///if (screen) fprintf(screen,"body tag %d: forcecm %f %f %f, torquecm %f %f %f\n",tag(ibody),fcm[ibody][0],fcm[ibody][1],fcm[ibody][2],torquecm[ibody][0],torquecm[ibody][1],torquecm[ibody][2]);
  }
}

/* ----------------------------------------------------------------------
   set space-frame coords and velocity of each atom in each rigid body
   set orientation and rotation of extended particles
   x = Q displace + Xcm, mapped back to periodic box
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

//NP adapted to loop also over ghosts
void FixMultisphere::set_xv()
{
    set_xv(LOOP_ALL);
}

/*-----------------------------------------------------------------------*/

void FixMultisphere::set_xv(int ghostflag)
{
  int ibody;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double vr[6];

  tagint *image = atom->image;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega_one = atom->omega;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  double **xcm = multisphere_.xcm_.begin();
  double **vcm = multisphere_.vcm_.begin();
  double **ex_space = multisphere_.ex_space_.begin();
  double **ey_space = multisphere_.ey_space_.begin();
  double **ez_space = multisphere_.ez_space_.begin();
  double **omega = multisphere_.omega_.begin();

  int nloop = 0;

  if(ghostflag == LOOP_ALL) nloop = nlocal+nghost;
  else if(ghostflag == LOOP_LOCAL) nloop = nlocal;
  else error->all(FLERR,"Illegal call to FixMultisphere::set_v");

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // set x and v of each atom

  for (int i = 0; i < nloop; i++) {

    /*NL*/ //if(screen && atom->tag[i] == 23497) fprintf(screen,"atom tag 23497 has body %d\n",body_[i]);

    if (body_[i] < 0) continue;
    ibody = map(body_[i]);

    if (ibody < 0) continue;

    xbox = (image[i] & IMGMASK) - IMGMAX;
    ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
    zbox = (image[i] >> IMG2BITS) - IMGMAX;

    /*NL*/// if (screen) fprintf(screen,"boxes %d %d %d image %d\n",xbox,ybox,zbox,image[i]);

    // save old positions and velocities for virial

    if (evflag) {
      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    // x = displacement from center-of-mass, based on body orientation
    // v = vcm + omega around center-of-mass

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],displace_[i],x[i]);

    v[i][0] = omega[ibody][1]*x[i][2] - omega[ibody][2]*x[i][1] + vcm[ibody][0];
    v[i][1] = omega[ibody][2]*x[i][0] - omega[ibody][0]*x[i][2] + vcm[ibody][1];
    v[i][2] = omega[ibody][0]*x[i][1] - omega[ibody][1]*x[i][0] + vcm[ibody][2];

    // add center of mass to displacement
    // map back into periodic box via xbox,ybox,zbox
    // for triclinic, would have to add in box tilt factors as well

    x[i][0] += xcm[ibody][0] - xbox*xprd;
    x[i][1] += xcm[ibody][1] - ybox*yprd;
    x[i][2] += xcm[ibody][2] - zbox*zprd;

    omega_one[i][0] = omega[ibody][0];
    omega_one[i][1] = omega[ibody][1];
    omega_one[i][2] = omega[ibody][2];

    /*NL*/ //bool eval = true; //50000 < update->ntimestep && 55000 > update->ntimestep;

    /*NL*/ //if(screen && tag(ibody) == 1 && eval) {
    /*NL*/ //     fprintf(screen,"step " BIGINT_FORMAT " proc %d atom tag %d ghost %s, boxes %d %d %d \n",
    /*NL*/ //                                      update->ntimestep,comm->me,atag[i],i>=nlocal?"yes":"no",
    /*NL*/ //                                      xbox,ybox,zbox);
    /*NL*/ //}


    /*NL*///if(screen && !domain->is_in_domain(x[i]))
    /*NL*///{
    /*NL*///    fprintf(screen,"particle %d, tag %d\n",i,atom->tag[i]);
    /*NL*///    printVec3D(screen,"omega_one[i]",omega_one[i]);
    /*NL*///    fprintf(screen,"part 1 %f\n",ex_space[ibody][0]*displace[i][0] +      ey_space[ibody][0]*displace[i][1] +      ez_space[ibody][0]*displace[i][2]);
    /*NL*///    fprintf(screen,"part 2 %f\n",xcm[ibody][0] - xbox*xprd);
    /*NL*///    fprintf(screen,"part 2 xbox*xprd %f\n",xbox*xprd);
    /*NL*///    fprintf(screen,"image %d\n",image[i]);
    /*NL*///    fprintf(screen,"xbox %d\n",xbox);
    /*NL*///    fprintf(screen,"xprd %f\n",xprd);
    /*NL*///    printVec3D(screen,"x",x[i]);
    /*NL*///    printVec3D(screen,"displace",displace[i]);
    /*NL*///    printVec3D(screen,"displace",displace[i]);
    /*NL*///    printVec3D(screen,"xcm",xcm[ibody]);
    /*NL*///    printVec3D(screen,"ex_space",ex_space[ibody]);
    /*NL*///    printVec3D(screen,"ey_space",ey_space[ibody]);
    /*NL*///    printVec3D(screen,"ez_space",ez_space[ibody]);
    /*NL*///    error->one(FLERR,"not in domain");
    /*NL*///}

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c final_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }
}

/* ----------------------------------------------------------------------
   set space-frame velocity of each atom in a rigid body
   set omega and angmom of extended particles
   v = Vcm + (W cross (x - Xcm))
------------------------------------------------------------------------- */

//NP adapted to loop also over ghosts
void FixMultisphere::set_v()
{
    set_v(LOOP_ALL);
}

/*-----------------------------------------------------------------------*/

void FixMultisphere::set_v(int ghostflag)
{
  int ibody;
  int xbox,ybox,zbox;
  double x0,x1,x2,v0,v1,v2,fc0,fc1,fc2,massone;
  double delta[3],vr[6];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **omega_one = atom->omega;
  int *type = atom->type;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  double **vcm = multisphere_.vcm_.begin();
  double **omega = multisphere_.omega_.begin();
  double **ex_space = multisphere_.ex_space_.begin();
  double **ey_space = multisphere_.ey_space_.begin();
  double **ez_space = multisphere_.ez_space_.begin();

  int nloop = 0;

  if(ghostflag == LOOP_ALL) nloop = nlocal+nghost;
  else if(ghostflag == LOOP_LOCAL) nloop = nlocal;
  else error->all(FLERR,"Illegal call to FixMultisphere::set_v");

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  // set v of each atom

  for (int i = 0; i < nloop; i++) {
    if (body_[i] < 0) continue;
    ibody = map(body_[i]);
    if (ibody < 0) continue;

    MathExtra::matvec(ex_space[ibody],ey_space[ibody],ez_space[ibody],displace_[i],delta);

    // save old velocities for virial

    if (evflag) {
      v0 = v[i][0];
      v1 = v[i][1];
      v2 = v[i][2];
    }

    v[i][0] = omega[ibody][1]*delta[2] - omega[ibody][2]*delta[1] + vcm[ibody][0];
    v[i][1] = omega[ibody][2]*delta[0] - omega[ibody][0]*delta[2] + vcm[ibody][1];
    v[i][2] = omega[ibody][0]*delta[1] - omega[ibody][1]*delta[0] + vcm[ibody][2];

    omega_one[i][0] = omega[ibody][0];
    omega_one[i][1] = omega[ibody][1];
    omega_one[i][2] = omega[ibody][2];

    // virial = unwrapped coords dotted into body constraint force
    // body constraint force = implied force due to v change minus f external
    // assume f does not include forces internal to body
    // 1/2 factor b/c initial_integrate contributes other half
    // assume per-atom contribution is due to constraint force on that atom

    if (evflag) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fc0 = massone*(v[i][0] - v0)/dtf - f[i][0];
      fc1 = massone*(v[i][1] - v1)/dtf - f[i][1];
      fc2 = massone*(v[i][2] - v2)/dtf - f[i][2];

      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;

      x0 = x[i][0] + xbox*xprd;
      x1 = x[i][1] + ybox*yprd;
      x2 = x[i][2] + zbox*zprd;

      vr[0] = 0.5*x0*fc0;
      vr[1] = 0.5*x1*fc1;
      vr[2] = 0.5*x2*fc2;
      vr[3] = 0.5*x0*fc1;
      vr[4] = 0.5*x0*fc2;
      vr[5] = 0.5*x1*fc2;

      v_tally(1,&i,1.0,vr);
    }
  }
}

/* ----------------------------------------------------------------------
   delete atoms belonging to deleted bodies
------------------------------------------------------------------------- */

void FixMultisphere::pre_exchange()
{
    AtomVec *avec = atom->avec;

    // reset last trigger for re-neigh
    next_reneighbor = -1;

    //NP have to do this here, b/c deletion at other place in code would
    //NP mess up paralleliztion

    double *delflag = fix_delflag_->vector_atom;
    /*NL*/ //double *existflag = fix_existflag_->vector_atom;
    int i = 0;

    /*NL*/ //double sum = vectorSumN(existflag,atom->nlocal);
    /*NL*/ //if (screen) fprintf(screen,"sum pre_exchange %f nlocal %d\n",sum,atom->nlocal);

    while(i < atom->nlocal)
    {
        //NP important to use round() here
        //NP never check if(double_variable) !!!!
        /*NL*/ //if (screen) fprintf(screen,"delfag particle %d: %f\n",atom->tag[i],delflag[i]);
        if(round(delflag[i]) == 1.)
        {
            /*NL*/ //if (screen) fprintf(screen,"step " BIGINT_FORMAT " proc %d deleting particle tag %d\n",update->ntimestep,comm->me,atom->tag[i]);
            avec->copy(atom->nlocal-1,i,1);
            atom->nlocal--;
        }
        else i++;
    }
}

/* ----------------------------------------------------------------------
   communicate body

   remap xcm of each rigid body back into periodic simulation box
   done during pre_neighbor so will be after call to pbc()
     and after fix_deform::pre_exchange() may have flipped box
   use domain->remap() in case xcm is far away from box
     due to 1st definition of rigid body or due to box flip
   if don't do this, then atoms of a body which drifts far away
     from a triclinic box will be remapped back into box
     with huge displacements when the box tilt changes via set_x()

   exchange bodies with stencil procs

   check for lost bodies and remove them
     also mark atoms belonging to lost bodies for deletion

   communicate displace, image
------------------------------------------------------------------------- */

void FixMultisphere::pre_neighbor()
{
    //NP reset corner ghost flag

    int nall = atom->nlocal + atom->nghost;
    double *corner_ghost = fix_corner_ghost_->vector_atom;
    vectorZeroizeN(corner_ghost,nall);

    /*NL*/ //if (screen) fprintf(screen,"step %d on proc %d: x-sublo %f x-subhi %f\n",update->ntimestep,comm->me,domain->sublo[0],domain->subhi[0]);
    /*NL*/ //if (screen) fprintf(screen,"step %d proc %d: nbody %d, nall %d\n",update->ntimestep,comm->me,n_body(),atom->nlocal+atom->nghost);

    //NP communicate displace, image, body to ghosts
    //NP since needed in set_xv(), set_v()
    //NP ok to communicate only once after re-neighboring
    //NP since can change only at processor exchange

    //NP communicate body first, since needed for check_lost_atoms()
    fw_comm_flag_ = MS_COMM_FW_BODY;
    forward_comm();

    //NP callback to fix remove
    //NP need this here before exchange()
    //NP since then the removal list generated at
    //NP FixRemove::pre_exchange() is still up-to-date
    for(size_t irem = 0; irem < fix_remove_.size(); irem++)
        (fix_remove_[irem])->delete_bodies();

    //NP re-map bodies, then exchange with stencil procs
    //NP need to do that here because body and atom data
    //NP must be sync'd before checking for lost atoms

    //NP   need  fw comm of image since not communicated on
    //NP   reneighboring steps in Comm::borders()

    //NP   need reverse comm of atom image flag since
    //NP   changed on ghosts by remap_bodies
    //NP fw comm image and displace
    //NP note that image has changed due to remap_bodies() call
    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();
    multisphere_.remap_bodies(body_);
    rev_comm_flag_ = MS_COMM_REV_IMAGE;
    reverse_comm();
    multisphere_.exchange();

    //NP after all bodies are communicated, need to find out how many we have
    multisphere_.calc_nbody_all();

    //NP need to re-set map as well
    //NP as exchange has happened so bodies deleted and added
    /*NL*/ //if (screen) fprintf(screen,"generating map after exchange");
    multisphere_.generate_map();

    //NP check for lost atoms (e.g. lost during atom communication via verlet,
    //NP fixed boundaries etc)
    //NP if atom is lost that belongs to rigid body, delete whole body
    //NP and mark other atoms belonging to this rigid body for deletion

    // set deletion flag
    // if any deleted atoms, do re-neigh in 100 steps at latest to remove
    // remainder particles
    double   *delflag =   fix_delflag_->vector_atom;
    double *existflag = fix_existflag_->vector_atom;
    vectorZeroizeN(delflag,atom->nlocal+atom->nghost);
    vectorZeroizeN(existflag,atom->nlocal+atom->nghost);

    if(multisphere_.check_lost_atoms(body_,delflag,existflag))
        next_reneighbor = update->ntimestep + 100;

    /*NL*/ //double sum = vectorSumN(existflag,atom->nlocal);
    /*NL*/ //if (screen) fprintf(screen,"sum pre %f nlocal %d\n",sum,atom->nlocal);

    //NP need to send deletion flag from ghosts to owners
    fix_delflag_->do_reverse_comm();
    fix_existflag_->do_reverse_comm();

    /*NL*/ //sum = vectorSumN(existflag,atom->nlocal);
    /*NL*/ //if (screen) fprintf(screen,"sum post %f nlocal %d\n",sum,atom->nlocal);

    //NP fw comm image and displace
    //NP note that image has changed due to remap_bodies() call
    fw_comm_flag_ = MS_COMM_FW_IMAGE_DISPLACE;
    forward_comm();

    // merge delflag and existflag

    int nlocal = atom->nlocal;
    delflag =   fix_delflag_->vector_atom;
    existflag = fix_existflag_->vector_atom;
    for(int i = 0; i < nlocal; i++)
    {
            delflag[i] = (round(existflag[i]) == 0) ? 1. : delflag[i];
    }
}

/* ----------------------------------------------------------------------
   count # of degrees-of-freedom removed by fix_rigid for atoms in igroup
------------------------------------------------------------------------- */

int FixMultisphere::dof(int igroup)
{
    int n = 0;

    if(0 == comm->me)
        error->warning(FLERR,"Energy calculated for multisphere particles is currently not correct");
    //error->all(FLERR,"For multisphere particles, please use compute ke and instead of the thermo keyword ke");

    return n;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixMultisphere::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*3 * sizeof(double);
  bytes += maxvatom*6 * sizeof(double);

  // add Multisphere memory usage

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixMultisphere::grow_arrays(int nmax)
{
    /*NL*/ //if (screen) fprintf(screen,"called with nmax %d\n",nmax);
    body_ = memory->grow(body_,nmax,"rigid:body_");
    memory->grow(displace_,nmax,3,"rigid:displace");
    atom->molecule = body_;
}

/* ----------------------------------------------------------------------
   extract values
------------------------------------------------------------------------- */
/*
void * FixMultisphere::extract(char *name, int &len1, int &len2)
{
    return multisphere_.extract(name,len1,len2);
}*/

/* ----------------------------------------------------------------------
   return attributes of a rigid body
   12 values per body
   xcm = 1,2,3; vcm = 4,5,6; fcm = 7,8,9; torque = 10,11,12
------------------------------------------------------------------------- */

double** FixMultisphere::get_dump_ref(int &nb, int &nprop, char* prop)
{
    error->one(FLERR,"TODO");
  return NULL;
}
