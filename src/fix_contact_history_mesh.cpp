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

#include "mpi.h"
#include "string.h"
#include "stdio.h"
#include "fix_contact_history_mesh.h"
#include "atom.h"
#include "fix_mesh_surface.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactHistoryMesh::FixContactHistoryMesh(LAMMPS *lmp, int narg, char **arg) :
  FixContactHistory(lmp, narg, arg),
  ipage1_(0),
  dpage1_(0),
  ipage2_(0),
  dpage2_(0),
  keeppage_(0),
  keepflag_(0),
  mesh_(0),
  fix_neighlist_mesh_(0),
  fix_nneighs_(0),
  build_neighlist_(true)
{
  //NP init performed in base class

  // parse args
  Fix *f = modify->find_fix_id(arg[iarg_++]);
  if(!f || strncmp(f->style,"mesh/surface",12) )
    error->fix_error(FLERR,this,"wrong ID for fix mesh/surface");
  mesh_ = (static_cast<FixMeshSurface*>(f))->triMesh();
  fix_neighlist_mesh_ = (static_cast<FixMeshSurface*>(f))->meshNeighlist();

  swap_ = new double[dnum_];

  // initial allocation of delflag
  keepflag_ = (bool **) memory->srealloc(keepflag_,atom->nmax*sizeof(bool *),
                                      "contact_history:keepflag");
}

/* ---------------------------------------------------------------------- */

FixContactHistoryMesh::~FixContactHistoryMesh()
{
  // delete locally stored arrays

  if(ipage1_) delete [] ipage1_;
  if(dpage1_) delete [] dpage1_;
  if(ipage2_) delete [] ipage2_;
  if(dpage2_) delete [] dpage2_;

  if(keeppage_) {
    for(int i = 0; i < comm->nthreads; i++) {
      delete keeppage_[i];
      keeppage_[i] = NULL;
    }
    delete [] keeppage_;
    keeppage_ = NULL;
  }

  //NP have been free'd above
  ipage_ = 0;
  dpage_ = 0;

  delete [] swap_;

  if(keepflag_) memory->sfree(keepflag_);
}

/* ---------------------------------------------------------------------- */

int FixContactHistoryMesh::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= MIN_PRE_FORCE;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_EXCHANGE;
  mask |= MIN_PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::init()
{
  FixContactHistory::init();

  char *fix_nneighs_name = new char[strlen(mesh_->mesh_id())+1+14];
  sprintf(fix_nneighs_name,"n_neighs_mesh_%s",mesh_->mesh_id());

  fix_nneighs_ = static_cast<FixPropertyAtom*>(modify->find_fix_property(fix_nneighs_name,"property/atom","scalar",0,0,this->style));

  delete [] fix_nneighs_name;
}

/* ----------------------------------------------------------------------
  create pages if first time or if neighbor pgsize/oneatom has changed
  note that latter could cause shear history info to be discarded
------------------------------------------------------------------------- */

void FixContactHistoryMesh::allocate_pages()
{
  if ((ipage_ == NULL) || (pgsize_ != neighbor->pgsize) || (oneatom_ != neighbor->oneatom)) {

	//NP evaluates to false at first time allocation
    bool use_first = ipage_ == ipage2_;

    delete [] ipage1_;
    delete [] dpage1_;
    delete [] ipage2_;
    delete [] dpage2_;

    if(keeppage_) {
      for(int i = 0; i < comm->nthreads; i++) {
        delete keeppage_[i];
        keeppage_[i] = NULL;
      }
      delete [] keeppage_;
      keeppage_ = NULL;
    }

    pgsize_ = neighbor->pgsize;
    oneatom_ = neighbor->oneatom;
    int nmypage = comm->nthreads;
    ipage1_ = new MyPage<int>[nmypage];
    dpage1_ = new MyPage<double>[nmypage];
    ipage2_ = new MyPage<int>[nmypage];
    dpage2_ = new MyPage<double>[nmypage];
    keeppage_ = new MyPage<bool>*[nmypage];
    for (int i = 0; i < nmypage; i++) {
      ipage1_[i].init(oneatom_,pgsize_);
      dpage1_[i].init(oneatom_,pgsize_);
      ipage2_[i].init(oneatom_,pgsize_);
      dpage2_[i].init(oneatom_,pgsize_);
    }

#if defined(_OPENMP)
    #pragma omp parallel
    {
      const int tid = omp_get_thread_num();
      // make sure page is allocated in memory near core
      keeppage_[tid] = new MyPage<bool>();
      keeppage_[tid]->init(oneatom_,pgsize_);
    }
#else
    for (int i = 0; i < nmypage; i++) {
      keeppage_[i] = new MyPage<bool>();
      keeppage_[i]->init(oneatom_,pgsize_);
    }
#endif

    if(use_first)
    {
        ipage_ = ipage1_;
        dpage_ = dpage1_;
    }
    else
    {
        ipage_ = ipage2_;
        dpage_ = dpage2_;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::setup_pre_exchange()
{
    pre_exchange();
}
/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_setup_pre_exchange()
{
  pre_exchange();
}

/* ----------------------------------------------------------------------
   sort contacthistory so need
------------------------------------------------------------------------- */

void FixContactHistoryMesh::pre_exchange()
{
    //NP recent_restart is set in case is called out of setup_pre_exchange
    if(!recent_restart)
        sort_contacts();

   // set maxtouch = max # of partners of any owned atom
   // bump up comm->maxexchange_fix if necessary

   int nlocal = atom->nlocal;

   maxtouch_ = 0;
   for (int i = 0; i < nlocal; i++) maxtouch_ = MAX(maxtouch_,npartner_[i]);
   comm->maxexchange_fix = MAX(comm->maxexchange_fix,(dnum_+1)*maxtouch_+1);
}

/* ----------------------------------------------------------------------
   need to execute here since neighlist is refreshed in setup_pre_force(int foo)
   if this is not done, data will get out-of-sync

   need not care about overwriting hist values since always have pointer
   to most current data
------------------------------------------------------------------------- */

void FixContactHistoryMesh::setup_pre_force(int dummy)
{
   pre_neighbor();
   pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::pre_neighbor()
{
    build_neighlist_ = true;
    /*NL*/ //fprintf(screen,"***building neigh list at step "BIGINT_FORMAT"\n",update->ntimestep);
}

/* ----------------------------------------------------------------------
   swap 2 pages data structures where contacthistory is stored
   do this here because exchange and sort can rag the data structure
------------------------------------------------------------------------- */

//NP this is called after FixNeighlistMesh::pre_force()
//NP because of the order of creation in FixWallGran::post_create()
//NP this is important so building new neigh list before refreshing contact hist
//NP contact hist is allocated for all neighs, no need to know # neighs before alloc

void FixContactHistoryMesh::pre_force(int dummy)
{
    if(!build_neighlist_)
        return;
    build_neighlist_ = false;

    //NP remove old contacts which are still in the history list
    //NP but not in the new mesh neigh list
    //NP this might be e.g. caused by particles moving through PBCs
    //NP in this case, contact "jumps" from one element to the other
    //NP do this here since have up-to-date mesh neighlist here
    //NP which was created in FixNeighlistMesh::pre_force()

    cleanUpContactJumps();

    //NP prev refers to time-steps before current rebuild
    //NP next refers to time-steps after current rebuild

    int nneighs_next;
    int nlocal = atom->nlocal;
    int *partner_prev;
    double *contacthistory_prev;

    MyPage<int>    *ipage_next = (ipage_ == ipage1_) ? ipage2_ : ipage1_;
    MyPage<double> *dpage_next = (dpage_ == dpage1_) ? dpage2_ : dpage1_;

    ipage_next->reset();
    dpage_next->reset();

    for (int i = 0; i < nlocal; i++)
    {
        nneighs_next = fix_nneighs_->get_vector_atom_int(i);

        //NP store pointer to previous data
        partner_prev = partner_[i];
        contacthistory_prev = contacthistory_[i];

        //NP get new storage for partner at next
        //NP get new storage for contact history at next
        partner_[i] = ipage_next->get(nneighs_next);
        vectorInitializeN(partner_[i],nneighs_next,-1);
        contacthistory_[i] = dpage_next->get(nneighs_next*dnum_);
        vectorZeroizeN(contacthistory_[i],nneighs_next*dnum_);

        /*NL*/ //fprintf(screen,"nneighs_next %d npartner_[i] %d\n",nneighs_next,npartner_[i]);

        if(npartner_[i] > nneighs_next)
        {
            /*NL*/ //fprintf(screen,"tag %d nneighs_next %d npartner_[i] %d\n",atom->tag[i],nneighs_next,npartner_[i]);
            error->one(FLERR,"internal error");
        }

        //NP copy from current to next
        //NP need to loop from 0..npartner_[i]-1 only
        //NP since list has been sorted at this point so that active contacts
        //NP are at the beginning of the list

        for(int ipartner = 0; ipartner < npartner_[i]; ipartner++)
        {
            if(partner_prev[ipartner] < 0)
                error->one(FLERR,"internal error");

            partner_[i][ipartner] = partner_prev[ipartner];
            vectorCopyN(&(contacthistory_prev[ipartner*dnum_]),&(contacthistory_[i][ipartner*dnum_]),dnum_);
            /*NL*/ //fprintf(screen,"new data %e %e %e\n",contacthistory_[i][ipartner*dnum_],contacthistory_[i][ipartner*dnum_+1],contacthistory_[i][ipartner*dnum_+2]);
        }
   }

   //NP important not to re-set pages here because unpack_exchange() might want to add to it

    //NP switch current and next
    ipage_ = ipage_next;
    dpage_ = dpage_next;
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::sort_contacts()
{
    int nlocal = atom->nlocal;
    int nneighs, first_empty, last_filled;

    /*NL*/ //fprintf(screen,"starting sort\n");

    for(int i = 0; i < nlocal; i++)
    {
        nneighs = fix_nneighs_->get_vector_atom_int(i);

        if(0 == nneighs)
            continue;

        do
        {
            first_empty = last_filled = -1;

            for(int j = 0; j < nneighs; j++)
            {
                if(-1 == first_empty && -1 == partner_[i][j])
                    first_empty = j;
                if(partner_[i][j] >= 0)
                    last_filled = j;
            }

            if(first_empty > -1 && last_filled > -1 && first_empty < last_filled)
                swap(i,first_empty,last_filled,true);
        }
        while(first_empty > -1 && last_filled > -1 && first_empty < last_filled);
    }

    /*NL*/ //fprintf(screen,"ending sort\n");
}

/* ----------------------------------------------------------------------
     mark all contacts for deletion
------------------------------------------------------------------------- */

void FixContactHistoryMesh::markAllContacts()
{
    int nlocal = atom->nlocal;
    keeppage_[0]->reset(true);

    for(int i = 0; i < nlocal; i++)
    {
      const int nneighs = fix_nneighs_->get_vector_atom_int(i);
      keepflag_[i] = keeppage_[0]->get(nneighs);
    }
}

/* ----------------------------------------------------------------------
     mark all contacts for deletion
------------------------------------------------------------------------- */

void FixContactHistoryMesh::resetDeletionPage(int tid)
{
    // keep pages are initalized with 0 (= false)
    keeppage_[tid]->reset(true);
}

/* ----------------------------------------------------------------------
     mark all contacts for deletion
------------------------------------------------------------------------- */

void FixContactHistoryMesh::markForDeletion(int tid, int ifrom, int ito)
{
    for(int i = ifrom; i < ito; i++)
    {
      const int nneighs = fix_nneighs_->get_vector_atom_int(i);
      keepflag_[i] = keeppage_[tid]->get(nneighs);
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::cleanUpContacts()
{
    cleanUpContacts(0, atom->nlocal);
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::cleanUpContacts(int ifrom, int ito)
{
    /*NL*/ //if(10024 == update->ntimestep) {
    /*NL*/ //       int iDeb = atom->map(DEBUG_P_TAG);
    /*NL*/ //       int nn = static_cast<int>(round(fix_nneighs_->vector_atom[iDeb]));
    /*NL*/ //       for(int kk = 0; kk < nn; kk++)
    /*NL*/ //          fprintf(screen,"       neigh PRE-pre-clean %d: tri id %d\n",kk,partner_[iDeb][kk]); }

    //NP delete contacts that no longer exist (i.e. delflag set)

    for(int i = ifrom; i < ito; i++)
    {
        const int nneighs = fix_nneighs_->get_vector_atom_int(i);

        for(int j = 0; j < nneighs; j++)
        {
            // delete values
            //NP do not swap with last in array, since paged data structure
            //NP doesn't support this
            if(!keepflag_[i][j])
            {
                if(partner_[i][j] > -1)
                {
                    npartner_[i]--;
                    /*NL*/ //if(DEBUG_P_TAG == atom->tag[i]) {fprintf(screen,"step "BIGINT_FORMAT" deleting contact: tri id %d, npartner %d, nneighs %d\n",update->ntimestep,partner_[i][j],npartner_[i],nneighs);
                    /*NL*/ //       for(int kk = 0; kk < nneighs; kk++)
                    /*NL*/ //         fprintf(screen,"    neigh %d: tri id %d\n",kk,partner_[i][kk]); }
                }
                partner_[i][j] = -1;
                vectorZeroizeN(&(contacthistory_[i][j*dnum_]),dnum_);
            }
        }
        /*NL*/ //if(DEBUG_P_TAG == atom->tag[i]) {
        /*NL*/ //       for(int kk = 0; kk < nneighs; kk++)
        /*NL*/ //          fprintf(screen,"    neigh post-clean %d: tri id %d\n",kk,partner_[i][kk]); }
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::cleanUpContactJumps()
{
    int nlocal = atom->nlocal;
    int iTri;

    //NP delete contacts that no longer exist (i.e. not found in mesh neigh list)
    //NP important to loop here from 0..npartner-1 because when this fct is called
    //NP nneighs is already the new #of neighbors

    for(int i = 0; i < nlocal; i++)
    {
        int ipartner = 0;
        while (ipartner < npartner_[i])
        {
            //NP contacts should be sorted at this point
            if(partner_[i][ipartner] < 0)
                error->one(FLERR,"internal error");

            iTri = mesh_->map(partner_[i][ipartner]);
            /*NL*/ //if(7998 == atom->tag[i])  fprintf(screen,"tri id %d found in old list %s\n",partner_[i][ipartner],fix_neighlist_mesh_->contactInList(iTri,i)?"true":"false");

            //NP delete also non-owned triangles in list
            //NP can happen when > 1 proc in periodic dimension
            if(iTri == -1 || (iTri > -1 && !fix_neighlist_mesh_->contactInList(iTri,i)))
            {
                /*NL*/ //if(7998 == atom->tag[i]) fprintf(screen,"step "BIGINT_FORMAT" deleting contact jump: tri id %d, npartner %d\n",update->ntimestep,partner_[i][ipartner],npartner_[i]);

                //NP delete non-existing contact
                //NP need to swap here so that there are no "holes" in the contact list
                //NP necessary since this is assumed in pre_force and old number of neighbors
                //NP is not known
                //NP so need to keep list without holes
                //NP do not swap delflag since not allocated at this point in parallel
                //NP since unpack_exchange() does not allocate delflag
                partner_[i][ipartner] = -1;
                vectorZeroizeN(&(contacthistory_[i][ipartner*dnum_]),dnum_);
                swap(i,ipartner,npartner_[i]-1,false);
                npartner_[i]--;

                /*NL*/ //error->one(FLERR,"catch");
            }
            else
                ipartner++;
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::reset_history()
{
    int nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
    {
        const int nneighs = fix_nneighs_->get_vector_atom_int(i);

        // zeroize values
        for(int j = 0; j < nneighs; j++)
            vectorZeroizeN(&(contacthistory_[i][j*dnum_]),dnum_);
    }
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_setup_pre_force(int dummy)
{
  if (*computeflag_) pre_force(0);
  *computeflag_ = 0;
}

/* ---------------------------------------------------------------------- */

void FixContactHistoryMesh::min_pre_force(int dummy)
{
  pre_force(0);
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistoryMesh::grow_arrays(int nmax)
{
  FixContactHistory::grow_arrays(nmax);
  keepflag_ = (bool **) memory->srealloc(keepflag_,nmax*sizeof(bool *),
                                      "contact_history:keepflag");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistoryMesh::copy_arrays(int i, int j, int delflag)
{
  // just copy pointers for partner and shearpartner
  // b/c can't overwrite chunk allocation inside ipage,dpage
  // incoming atoms in unpack_exchange just grab new chunks
  // so are orphaning chunks for migrating atoms
  // OK, b/c will reset ipage,dpage on next reneighboring

  FixContactHistory::copy_arrays(i,j,delflag);
  keepflag_[j] = keepflag_[i];
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistoryMesh::unpack_exchange(int nlocal, double *buf)
{
  // allocate new chunks from ipage,dpage for incoming values
  //NP do not allocate new chunk for delflag, since value is not needed
  //NP between time-steps

  int m = 0;
  const int nneighs = fix_nneighs_->get_vector_atom_int(nlocal);

  /*NL*/ //fprintf(screen,"unpacking (id %s), nlocal %d, nneighs %d npartner %d \n",id,nlocal,nneighs,static_cast<int> (buf[m]));

  //NP allocate nneigh storage positions instead of npartner in base class

  npartner_[nlocal] = ubuf(buf[m++]).i;
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(nneighs);
  contacthistory_[nlocal] = dpage_->get(dnum_*nneighs);

  //NP unpack values for contacts
  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = ubuf(buf[m++]).i;
    /*NL*/ //fprintf(screen,"atom %d: unpacking partner %d\n",atom->tag[nlocal],partner_[nlocal][n]);
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = buf[m++];
      /*NL*/ //fprintf(screen,"atom %d: unpacking hist %e\n",atom->tag[nlocal],contacthistory_[nlocal][n*dnum_+d]);
    }
  }

  //NP set initial values for neighs
  for (int n = npartner_[nlocal]; n < nneighs; n++) {
    partner_[nlocal][n] = -1;
    for (int d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = 0.;
    }
  }

  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactHistoryMesh::unpack_restart(int nlocal, int nth)
{
  // ipage = NULL if being called from granular pair style init()

  if (ipage_ == NULL) allocate_pages();

  // skip to Nth set of extra values

  double **extra = atom->extra;

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  // allocate new chunks from ipage,dpage for incoming values
  //NP allocate nneigh storage positions instead of npartner in base class

  int d;

  npartner_[nlocal] = ubuf(extra[nlocal][m++]).i;
  maxtouch_ = MAX(maxtouch_,npartner_[nlocal]);
  partner_[nlocal] = ipage_->get(npartner_[nlocal]);
  contacthistory_[nlocal] = dpage_->get(npartner_[nlocal]*dnum_);

  /*NL*/ //fprintf(screen,"npartner_[nlocal] %d\n",npartner_[nlocal]);

  //NP unpack values for contacts
  for (int n = 0; n < npartner_[nlocal]; n++) {
    partner_[nlocal][n] = ubuf(extra[nlocal][m++]).i;
    for (d = 0; d < dnum_; d++) {
      contacthistory_[nlocal][n*dnum_+d] = extra[nlocal][m++];
      /*NL*/ //fprintf(screen,"unpacking %e\n",contacthistory_[nlocal][n*dnum_+d]);
    }
  }
}

/* ----------------------------------------------------------------------
   pack state of Fix into one write
------------------------------------------------------------------------- */

void FixContactHistoryMesh::write_restart(FILE *fp)
{
    FixContactHistory::write_restart(fp);

    //NP sort contacts so only contacts, not neighs are written
    //NP pack_*() implementations in base class operate in npartner_
    //NP not nneighs
    sort_contacts();
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixContactHistoryMesh::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(int *);
  bytes += nmax * sizeof(double *);

  int nmypage = comm->nthreads;
  for (int i = 0; i < nmypage; i++) {
    bytes += ipage1_[i].size();
    bytes += dpage1_[i].size();
    bytes += ipage2_[i].size();
    bytes += dpage2_[i].size();
    bytes += keeppage_[i]->size();
  }

  return bytes;
}
