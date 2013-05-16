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

#define DELTA 10000

#include "multisphere.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "vector_liggghts.h"

/* ----------------------------------------------------------------------
   constructor / destructor
------------------------------------------------------------------------- */

Multisphere::Multisphere(LAMMPS *lmp) :
  Pointers(lmp),
  customValues_(*(new CustomValueTracker(lmp))),

  nbody_(0),
  nbody_all_(0),
  mapTagMax_(0), //NP this makes tags start at 1!!!
  mapArray_(0),

  id_ (*customValues_.addElementProperty< ScalarContainer<int> >("id","comm_exchange_borders"/*ID does never change*/,"frame_invariant","restart_yes")),

  xcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("xcm","comm_exchange_borders","frame_invariant", "restart_yes")),
  vcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("vcm","comm_exchange_borders","frame_invariant", "restart_yes")),
  fcm_          (*customValues_.addElementProperty< VectorContainer<double,3> >("fcm","comm_none","frame_invariant", "restart_no")),
  torquecm_     (*customValues_.addElementProperty< VectorContainer<double,3> >("torque","comm_none","frame_invariant", "restart_no")),
  dragforce_cm_ (*customValues_.addElementProperty< VectorContainer<double,3> >("dragforce","comm_none","frame_invariant", "restart_no")),

  angmom_ (*customValues_.addElementProperty< VectorContainer<double,3> >("angmom","comm_exchange_borders","frame_invariant", "restart_yes")),
  omega_  (*customValues_.addElementProperty< VectorContainer<double,3> >("omega","comm_exchange_borders","frame_invariant", "restart_yes")),
  quat_   (*customValues_.addElementProperty< VectorContainer<double,4> >("quat","comm_exchange_borders","frame_invariant", "restart_yes")),

  density_   (*customValues_.addElementProperty< ScalarContainer<double> >("density","comm_exchange_borders","frame_invariant","restart_yes")),
  masstotal_ (*customValues_.addElementProperty< ScalarContainer<double> >("masstotal","comm_exchange_borders","frame_invariant","restart_yes")),
  inertia_   (*customValues_.addElementProperty< VectorContainer<double,3> >("inertia","comm_exchange_borders","frame_invariant", "restart_yes")),
  ex_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ex_space","comm_exchange_borders","frame_invariant", "restart_yes")),
  ey_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ey_space","comm_exchange_borders","frame_invariant", "restart_yes")),
  ez_space_  (*customValues_.addElementProperty< VectorContainer<double,3> >("ez_space","comm_exchange_borders","frame_invariant", "restart_yes")),

  nrigid_    (*customValues_.addElementProperty< ScalarContainer<int> >("nrigid","comm_exchange_borders","frame_invariant", "restart_yes")),

  imagebody_ (*customValues_.addElementProperty< ScalarContainer<int> >("imagebody","comm_exchange_borders","frame_invariant", "restart_yes")),
  remapflag_ (*customValues_.addElementProperty< VectorContainer<int,4> >("remapflag","comm_none","frame_invariant", "restart_no")),

  fflag_ (*customValues_.addElementProperty< VectorContainer<bool,3> >("fflag","comm_exchange_borders","frame_invariant", "restart_yes")),
  tflag_ (*customValues_.addElementProperty< VectorContainer<bool,3> >("tflag","comm_exchange_borders","frame_invariant", "restart_yes")),

  start_step_ (*customValues_.addElementProperty< ScalarContainer<int> >("start_step","comm_exchange_borders","frame_invariant", "restart_yes")),
  v_integrate_(*customValues_.addElementProperty< VectorContainer<double,3> >("v_integrate","comm_exchange_borders","frame_invariant", "restart_yes")),

  r_bound_       (*customValues_.addElementProperty< ScalarContainer<double> >("r_bound","comm_exchange_borders","frame_invariant", "restart_yes")),
  xcm_to_xbound_ (*customValues_.addElementProperty< VectorContainer<double,3> >("xcm_to_xbound","comm_exchange_borders","frame_invariant", "restart_yes"))
{

}

Multisphere::~Multisphere()
{
    delete &customValues_;

    // deallocate map memory if exists
    if(mapArray_) clear_map();
}

/* ----------------------------------------------------------------------
   add a new body
------------------------------------------------------------------------- */

void Multisphere::add_body(int nspheres, double *xcm_ins, double *xcm_to_xbound_ins,
               double r_bound_ins, double *v_ins, double *omega_ins, double mass_ins,double dens_ins,
               double *inertia_ins, double *ex_space_ins, double *ey_space_ins, double *ez_space_ins,
               double **displace_ins, int start_step_ins,double *v_integrate_ins)
{
    int n = nbody_;

    // set initialize ID for element
    // ID starts from 0
    //NP -1 means no ID given, will be done by id_extend()
    id_.add(-1);

    bool flags[3] = {true,true,true};
    double zerovec[3] = {0.,0.,0.};

    xcm_.add(xcm_ins);
    vcm_.add(v_ins);
    fcm_.addZero();
    torquecm_.addZero();
    dragforce_cm_.addZero();

    angmom_.addZero();
    omega_.add(omega_ins);
    quat_.addZero();

    density_.add(dens_ins);
    masstotal_.add(mass_ins);
    inertia_.add(inertia_ins);
    ex_space_.add(ex_space_ins);
    ey_space_.add(ey_space_ins);
    ez_space_.add(ez_space_ins);

    nrigid_.add(nspheres);
    imagebody_.add((512 << 20) | (512 << 10) | 512);
    remapflag_.addZero();

    fflag_.add(flags);
    tflag_.add(flags);

    start_step_.add(start_step_ins);
    if(v_integrate_ins)
        v_integrate_.add(v_integrate_ins);
    else
        v_integrate_.add(zerovec);

    r_bound_.add(r_bound_ins);
    xcm_to_xbound_.add(xcm_to_xbound_ins);

    // calculate q and ang momentum

    MathExtra::exyz_to_q
    (
        ex_space_(n),ey_space_(n),ez_space_(n),
        quat_(n)
    );

    MathExtraLiggghts::angmom_from_omega
    (
        omega_(n),
        ex_space_(n),ey_space_(n),ez_space_(n),
        inertia_(n),
        angmom_(n)
    );

    /*NL*/// printVec3D(screen,"omega_(n)",omega_(n));
    /*NL*/// printVec3D(screen,"angmom_(n)",angmom_(n));

    // increase local body counter
    nbody_++;

    //NP need to set body, displace for new particles

    //NP need to set x,v,omega for new particles
    //NP based on xcm, vcm, omega

    //NP this is done via FixMultisphere::add_body_finalize()
}

/* ----------------------------------------------------------------------
   remap bodies
------------------------------------------------------------------------- */

void Multisphere::remap_bodies(int *body)
{

  int original,oldimage,newimage;

  // adjust body

  for (int ibody = 0; ibody < nbody_; ibody++) {
    original = imagebody_(ibody);
    domain->remap(xcm_(ibody),imagebody_(ibody));

    if (original == imagebody_(ibody)) remapflag_(ibody)[3] = 0;
    else {
      oldimage = original & 1023;
      newimage = imagebody_(ibody) & 1023;
      remapflag_(ibody)[0] = newimage - oldimage;
      oldimage = (original >> 10) & 1023;
      newimage = (imagebody_(ibody) >> 10) & 1023;
      remapflag_(ibody)[1] = newimage - oldimage;
      oldimage = original >> 20;
      newimage = imagebody_(ibody) >> 20;
      remapflag_(ibody)[2] = newimage - oldimage;
      remapflag_(ibody)[3] = 1;
    }
  }

  //NP adjust image flags of any atom in a rigid body whose xcm was remapped
  int *atomimage = atom->image;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;

  int ibody,idim,otherdims;

  for (int i = 0; i < nlocal+nghost; i++)
  {
    if(body[i] < 0) continue;
    ibody = map(body[i]);

    if (ibody < 0) continue;
    if (remapflag_(ibody)[3] == 0) continue;

    if (remapflag_(ibody)[0]) {
      idim = atomimage[i] & 1023;
      otherdims = atomimage[i] ^ idim;
      idim -= remapflag_(ibody)[0];
      idim &= 1023;
      atomimage[i] = otherdims | idim;
    }
    if (remapflag_(ibody)[1]) {
      idim = (atomimage[i] >> 10) & 1023;
      otherdims = atomimage[i] ^ (idim << 10);
      idim -= remapflag_(ibody)[1];
      idim &= 1023;
      atomimage[i] = otherdims | (idim << 10);
    }
    if (remapflag_(ibody)[2]) {
      idim = atomimage[i] >> 20;
      otherdims = atomimage[i] ^ (idim << 20);
      idim -= remapflag_(ibody)[2];
      idim &= 1023;
      atomimage[i] = otherdims | (idim << 20);
    }
  }
}

/* ----------------------------------------------------------------------
   add unique ids to any body with id = -1
   new ids are grouped by proc and start after max current tag
   called after creating new atoms
------------------------------------------------------------------------- */

void Multisphere::id_extend_body_extend(int *body)
{
  int idmax, idmax_all;
  int nlocal = atom->nlocal;

  // calc total # of bodies
  calc_nbody_all();

  // return if no bodies are present
  if(nbody_all_ == 0)
    return;

  // idmax = max id for all bodies

  idmax = id_.max();
  MPI_Max_Scalar(idmax,idmax_all,this->world);

  // mapTagMax_ cannot get smaller - so ensure IDs are given only once
  //NP initialially mapTagMax_=0 --> so tags start at 1

  mapTagMax_ = MathExtraLiggghts::max(mapTagMax_,idmax_all);

  /*NL*/ //fprintf(this->screen,"id_extend proc %d: idmax %d, idmax_all %d mapTagMax_ %d nbody_ %d\n",comm->me,idmax,idmax_all,mapTagMax_,nbody_);

  // noid = # of bodies I own with no id (id = -1)
  // noid_sum = # of total bodies on procs <= me with no tag
  // also check if # of atoms with no body id is consistent
  // nobody = # of atoms newly inserted with no body associated (body = -2)
  // nobody_check = # of atoms that should be newly inserted

  int noid = 0;
  int nobody = 0, nobody_first = 0, nobody_check = 0;
  for (int i = 0; i < nbody_; i++)
  {
    if (id_(i) == -1)
    {
        noid++;
        nobody_check += nrigid_(i);
    }
  }

  for(int i = 0; i < nlocal; i++)
  {
    if(body[i] == -2)
    {
       if(nobody == 0)
         nobody_first = i;
       nobody++;
    }
  }

  if(nobody != nobody_check)
    error->one(FLERR,"Internal error: # of atoms with no associated body inconsistent");

  int noid_sum;
  MPI_Scan(&noid,&noid_sum,1,MPI_INT,MPI_SUM,world);

  // itag = 1st new tag that my untagged bodies should use
  // give atoms body ID as well
  //NP assume all new particles are at the end of the list
  //NP body === -2 means means belongs to body, but does not know yet to which
  //NP important to differentiate between -1 (does not belong to any body)
  //NP and -2 - in case spheres and multi-spheres are mixed

  int itag = mapTagMax_ + noid_sum - noid + 1;
  for (int ibody = 0; ibody < nbody_; ibody++)
  {
    if (id_(ibody) == -1)
    {
        id_(ibody) = itag;
        /*NL*/ //fprintf(screen,"proc %d setting id of local body %d to %d\n",comm->me,ibody,id_(ibody));

        if(nobody_first == nlocal-1)
            error->one(FLERR,"Internal error: atom body id inconsistent");

        for(int iatom = nobody_first; iatom < nobody_first+nrigid_(ibody); iatom++)
        {
            if(body[iatom] != -2)
                error->one(FLERR,"Internal error: atom body id inconsistent");
            body[iatom] = itag;
            /*NL*/ //if((itag==39||itag==40)) fprintf(screen,"atom tag %d belongs to body tag %d\n",atom->tag[iatom],itag);
        }
        nobody_first += nrigid_(ibody);

        // search for next particle with no body associated
        while(nobody_first < nlocal-1 && body[nobody_first] != -2)
            nobody_first++;

        itag++;
    }
  }
}

/* ----------------------------------------------------------------------
   clear and generate a global map for global-local lookup
------------------------------------------------------------------------- */

void Multisphere::clear_map()
{
    // deallocate old memory
    memory->destroy(mapArray_);
    mapArray_ = NULL;
}

void Multisphere::generate_map()
{
    int idmax, idmax_all;

    // deallocate old memory if exists
    if(mapArray_) clear_map();

    if(nbody_all_ == 0)
        return;

    // get max ID of all proc
    idmax = id_.max();
    MPI_Max_Scalar(idmax,idmax_all,world);
    mapTagMax_ = MathExtraLiggghts::max(mapTagMax_,idmax_all);

    /*NL*/ //fprintf(this->screen,"generate_map proc %d: idmax %d, mapTagMax_ %d nbody_ %d\n",comm->me,idmax,mapTagMax_,nbody_);

    // alocate and initialize new array
    // IDs start at 1, have to go up to (inclusive) mapTagMax_
    //NP mapTagMax_ can be -1 if only newly inserted bodies are present
    //NP in this case, nothing happens here
    memory->create(mapArray_,mapTagMax_+1,"Multisphere:mapArray_");
    for(int i = 0; i < mapTagMax_+1; i++)
        mapArray_[i] = -1;

    // build map
    for (int i = nbody_-1; i >= 0; i--)
    {
        /*NL*/// fprintf(this->screen,"id_(i) %d\n",id_(i));
        mapArray_[id_(i)] = i;
    }
}

/* ----------------------------------------------------------------------
   check for lost atoms
------------------------------------------------------------------------- */

bool Multisphere::check_lost_atoms(int *body, double *atom_delflag)
{
    int body_tag,ibody,i;
    int nall = atom->nlocal + atom->nghost;
    int deleted = 0;

    int *nrigid_current = new int[nbody_];
    int *delflag = new int[nbody_];
    vectorZeroizeN(nrigid_current,nbody_);
    vectorZeroizeN(delflag,nbody_);

    /*NL*///fprintf(screen,"CHECKING lost atoms for %d bodies on proc %d\n",nbody_,comm->me);

    //NP check if nrigid is consistent for each body
    //NP   with atoms body data
    //NP if inconstistent, delete body and particles
    //NP must be called on reneighboring step
    for(i = 0; i < nall; i++)
    {
        body_tag = body[i];
        /*NL*///fprintf(screen,"proc %d step %d particle id %d, local %d (nlocal %d nghost %d) body_tag %d map(body_tag) %d is_periodic_ghost_of_owned %d\n",
        /*NL*///                comm->me,update->ntimestep,atom->tag[i],i,atom->nlocal,atom->nghost,body_tag,map(body_tag),domain->is_periodic_ghost_of_owned(i));

        //NP if particle belongs to a rigid body and the body is owned on this proc
        //NP specific case if ghost in periodic system with procgrid[idim] == 1
        //NP have to exclude double counts for this case
        if(body_tag >= 0 && map(body_tag) >= 0 && 0 == domain->is_periodic_ghost_of_owned(i))
            nrigid_current[map(body_tag)]++;
    }

    //NP mark bodies for deletion
    //NP do not actually delete them yet b/c need map function later
    for(ibody = 0; ibody < nbody_; ibody++)
    {
        if(nrigid_current[ibody] > nrigid_(ibody))
        {
            /*NL*/fprintf(screen,"proc %d FLAGGING removal of body tag %d ibody %d on step "BIGINT_FORMAT", nrigid_current %d nrigid %d\n",
            /*NL*/                comm->me,tag(ibody),ibody,update->ntimestep,nrigid_current[ibody],nrigid_(ibody));
            error->one(FLERR,"Internal error in multisphere method");
        }
        if(nrigid_current[ibody] != nrigid_(ibody))
        {
            /*NL*///fprintf(screen,"proc %d FLAGGING removal of body tag %d ibody %d on step "BIGINT_FORMAT", nrigid_current %d nrigid %d\n",
            /*NL*///                comm->me,tag(ibody),ibody,update->ntimestep,nrigid_current[ibody],nrigid_(ibody));
            delflag[ibody] = 1;
            /*NL*///error->one(FLERR,"end");
        }
    }

    //NP mark atoms belonging to bodies to be deleted for deletion as well
    for(i = 0; i < nall; i++)
    {
        body_tag = body[i];

        if(body_tag < 0) continue;

        ibody = map(body_tag);

        //NP remove if particle belongs to body and body is owned on this proc
        //NP also re-set body flag
        if(ibody >= 0 && delflag[ibody])
        {
           /*NL*///fprintf(screen,"step %d: proc %d deleting atom %d, belongs to body %d\n",update->ntimestep,comm->me,atom->tag[i],body_tag);
           atom_delflag[i] = 1.;
           body[i] = -1;
           deleted = 1;
        }
    }

    //NP delete bodies
    ibody = 0;
    while(ibody < nbody_)
    {
        if(delflag[ibody] == 1)
        {
            /*NL*///fprintf(screen,"DELETING body tag %d ibody %d on step %d\n",tag(ibody),ibody,update->ntimestep);
            delflag[ibody] = delflag[nbody_-1];
            remove_body(ibody);
        }
        else ibody++;
    }

    calc_nbody_all();

    MPI_Max_Scalar(deleted,world);

    delete []nrigid_current;
    delete []delflag;

    if(deleted == 1) return true;
    return false;
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   if don't recognize name, return NULL

   important: returns GLOBAL lengths (tag_max)
              because this is required by CfdDatacouplingMPI
------------------------------------------------------------------------- */

void *Multisphere::extract(char *name, int &len1, int &len2)
{
    // scalars

    len1 = len2 = 1;

    if (strcmp(name,"nbody") == 0) return (void *) &nbody_;
    if (strcmp(name,"nbody_all") == 0) return (void *) &nbody_all_;

    // per-body properties

    len1 = mapTagMax_;
    ContainerBase *cb = customValues_.getElementPropertyBase(name);

    if(NULL == cb)
    {
        len1 = len2 = -1;
        return NULL;
        /*
        fprintf(screen,"ERROR: Property %s not found, which is required for multi-sphere\n",name);
        error->all(FLERR,"Property required for multi-sphere not found");
        */
    }

    len2 = cb->lenVec();
    if(cb->nVec() != 1)
       error->all(FLERR,"Internal error, cannot use multi-vector containers");
    return cb->begin_slow_dirty();
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   only return if data of type double**
   if don't recognize name or not double ** data, return NULL
------------------------------------------------------------------------- */

double** Multisphere::extract_double_vector(char *name)
{
    VectorContainer<double,3> *cb = customValues_.getElementProperty<VectorContainer<double,3> >(name);
    if(!cb) return 0;
    return cb->begin();
}

/* ----------------------------------------------------------------------
   return a pointer to a named internal variable
   only return if data of type double*
   if don't recognize name or not double * data, return NULL
------------------------------------------------------------------------- */

double* Multisphere::extract_double_scalar(char *name)
{
    ScalarContainer<double> *cb = customValues_.getElementProperty<ScalarContainer<double> >(name);
    if(!cb) return 0;
    return cb->begin();
}
