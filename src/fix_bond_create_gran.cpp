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
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "fix_bond_create_gran.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define MAX(A,B) (((A) > (B)) ? (A) : (B))

/* ---------------------------------------------------------------------- */

FixBondCreateGran::FixBondCreateGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all(FLERR,"Illegal fix bond/create command");

  MPI_Comm_rank(world,&me);

  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal fix bond/create command");

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extvector = 0;

  iatomtype = atoi(arg[4]);
  jatomtype = atoi(arg[5]);
  double cutoff = atof(arg[6]);
  btype = atoi(arg[7]);
  newperts = atoi(arg[8]);

  if (iatomtype < 1 || iatomtype > atom->ntypes ||
      jatomtype < 1 || jatomtype > atom->ntypes)
    error->all(FLERR,"Invalid atom type in fix bond/create command");
  if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/create command");
  if (btype < 1 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in fix bond/create command");

  cutsq = cutoff*cutoff;

  // optional keywords

  imaxbond = 0;
  inewtype = iatomtype;
  jmaxbond = 0;
  jnewtype = jatomtype;
  fraction = 1.0;
  int seed = 12345;

  int iarg = 9;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"iparam") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      imaxbond = atoi(arg[iarg+1]);
      inewtype = atoi(arg[iarg+2]);
      if (imaxbond < 0) error->all(FLERR,"Illegal fix bond/create command");
      if (inewtype < 1 || inewtype > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix bond/create command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"jparam") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      jmaxbond = atoi(arg[iarg+1]);
      jnewtype = atoi(arg[iarg+2]);
      if (jmaxbond < 0) error->all(FLERR,"Illegal fix bond/create command");
      if (jnewtype < 1 || jnewtype > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix bond/create command");
      iarg += 3;
    } else if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix bond/create command");
      fraction = atof(arg[iarg+1]);
      seed = atoi(arg[iarg+2]);
      if (fraction < 0.0 || fraction > 1.0)
        error->all(FLERR,"Illegal fix bond/create command");
      if (seed <= 0) error->all(FLERR,"Illegal fix bond/create command");
      iarg += 3;
    } else error->all(FLERR,"Illegal fix bond/create command");
  }

  // error check

  if (atom->molecular == 0)
    error->all(FLERR,"Cannot use fix bond/create with non-molecular systems");
  if (iatomtype == jatomtype &&
      ((imaxbond != jmaxbond) || (inewtype != jnewtype)))
    error->all(FLERR,"Inconsistent iparam/jparam values in fix bond/create command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // perform initial allocation of atom-based arrays
  // register with Atom class
  // bondcount values will be initialized in setup()

  bondcount = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  countflag = 0;

  // set comm sizes needed by this fix

  comm_forward = newperts + 2; //NP modified C.K. 2;
  comm_reverse = newperts + 1; //NP modified C.K. 2;

  // allocate arrays local to this fix

  nmax = 0;
  npartner = NULL;
  partner = NULL;
  probability = NULL;
  //NP modified C.K. distsq = NULL;

  // zero out stats

  createcount = 0;
  createcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondCreateGran::~FixBondCreateGran()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete random;

  // delete locally stored arrays

  memory->sfree(bondcount);
  memory->sfree(npartner);
  memory->destroy(partner); //NP modified C.K.
  memory->sfree(probability); //NP modified C.K.
  //NP modified C.K. memory->sfree(partner);
  //NP modified C.K. memory->sfree(distsq);

  //NP do _not_  delete this fix here - should stay active
  //NP modify->delete_fix("propagate_bonds_gran");
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::post_create()
{
    //NP register a fix to propagate granular bonds across processors
    char* fixarg[4];

    fixarg[0]="propagate_bonds_gran";
    fixarg[1]="all";
    fixarg[2]="bond/propagate/gran";
    modify->add_fix(3,fixarg);
}


/* ---------------------------------------------------------------------- */

int FixBondCreateGran::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixBondCreateGran::modify_param(int narg, char **arg)
{
    if(narg != 2) error->all(FLERR,"Illegal fix_modify command");
    if(strcmp(arg[0],"every") == 0) nevery = atoi(arg[1]);
    else error->all(FLERR,"Illegal fix_modify command");
    return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::init()
{
  if(force->pair == NULL) error->all(FLERR,"Fix bond/create cutoff is longer than pairwise cutoff");

  if(!(force->bond_match("gran")))
     error->all(FLERR,"Fix bond/create can only be used together with dedicated 'granular' bond styles");

  // check cutoff for iatomtype,jatomtype - cutneighsq is used here
  double cutsq_limit = sqrt(force->pair->cutsq[iatomtype][jatomtype]) + neighbor->skin;
  cutsq_limit *= cutsq_limit;
  if (force->pair == NULL || cutsq > cutsq_limit)
    error->all(FLERR,"Fix bond/create cutoff is longer than pairwise cutoff");
/*
  // require special bonds = 0,1,1

  int flag = 0;
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) flag = 1;
  if (force->special_coul[1] != 0.0 || force->special_coul[2] != 1.0 ||
      force->special_coul[3] != 1.0) flag = 1;
  if (flag) error->all(FLERR,"Fix bond/create requires special_bonds = 0,1,1");

  // warn if angles, dihedrals, impropers are being used

  if (force->angle || force->dihedral || force->improper) {
    if (me == 0)
      error->warning(FLERR,"Created bonds will not create angles, "
                     "dihedrals, or impropers");
  }
*/
  // need a half neighbor list, built when ever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::setup(int vflag)
{
  int i,j,m;

  // compute initial bondcount if this is first run
  // can't do this earlier, like in constructor or init, b/c need ghost info

  if (countflag) return;
  countflag = 1;

  // count bonds stored with each bond I own
  // if newton bond is not set, just increment count on atom I
  // if newton bond is set, also increment count on atom J even if ghost
  // bondcount is long enough to tally ghost atom counts

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

  for (i = 0; i < nall; i++) bondcount[i] = 0;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] == btype) {
        bondcount[i]++;
        if (newton_bond) {
          m = atom->map(bond_atom[i][j]);
          if (m < 0)
            error->one(FLERR,"Could not count initial bonds in fix bond/create");
          bondcount[m]++;
        }
      }
    }

  // if newton_bond is set, need to sum bondcount

  commflag = 0;
  if (newton_bond) comm->reverse_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::post_integrate()
{
  int i,j,k,m,ii,jj,inum,jnum,itype,jtype,n1,n3,possible;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,min,max;
  int *ilist,*jlist,*numneigh,**firstneigh,*slist;
  int flag;

  if (nevery == 0 || update->ntimestep % nevery ) return;

  // need updated ghost atom positions

  comm->forward_comm();

  // forward comm of bondcount, so ghosts have it

  commflag = 0;
  comm->forward_comm_fix(this);

  // resize bond partner list and initialize it
  // probability array overlays distsq array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(partner,nmax,newperts,"fix_bond_create:partner"); //NP modified C.K.
    npartner = (int*) memory->srealloc(npartner,nmax*sizeof(int),"fix_bond_create:npartner"); //NP modified C.K.
    probability = (double *)memory->srealloc(probability,nmax*sizeof(double),"fix_bond_create:probability"); //NP modified C.K.
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    for(int j = 0; j< newperts; j++) partner[i][j] = 0; //NP modified C.K.
    npartner[i] = 0;  //NP modified C.K.
    probability[i] = 1.; //NP modified C.K.
  }

  // loop over neighbors of my atoms
  // each atom sets one closest eligible partner atom ID to bond with

  double **x = atom->x;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  flag = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];

      possible = 0;
      if (itype == iatomtype && jtype == jatomtype) {
        if ((imaxbond == 0 || bondcount[i] < imaxbond) &&
            (jmaxbond == 0 || bondcount[j] < jmaxbond))
          possible = 1;
      } else if (itype == jatomtype && jtype == iatomtype) {
        if ((jmaxbond == 0 || bondcount[i] < jmaxbond) &&
            (imaxbond == 0 || bondcount[j] < imaxbond))
          possible = 1;
      }
      if (!possible) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;

      if(already_bonded(i,j)) continue;

      if(npartner[i]==newperts || npartner[j]==newperts)
      {
          flag = 1;
          continue;
      }

      //NP modified C.K.
      partner[i][npartner[i]] = tag[j];
      partner[j][npartner[j]] = tag[i];
      npartner[i]++;
      npartner[j]++;

    }
  }

  if(flag) error->warning(FLERR,"Could not generate all possible bonds");

  // reverse comm of distsq and partner
  // not needed if newton_pair off since I,J pair was seen by both procs

  commflag = 1;
  if (force->newton_pair) comm->reverse_comm_fix(this);

  // each atom now knows its partners
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  if (fraction < 1.0) {
    for (i = 0; i < nlocal; i++)
      if (npartner[i]) probability[i] = random->uniform(); //NP modified C.K.
  }

  commflag = 1;
  comm->forward_comm_fix(this);

  // create bonds for atoms I own
  // if other atom is owned by another proc, it should create same bond
  // if both atoms list each other as bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  int **special = atom->special;
  int newton_bond = force->newton_bond;
  int n_bondhist = atom->n_bondhist;
  double ***bond_hist = atom->bond_hist;

  int ncreate = 0;
  for (i = 0; i < nlocal; i++) {
    if (npartner[i] == 0) continue; //NP modified C.K.

    for(k = 0; k < npartner[i]; k++)
    {
        j = atom->map(partner[i][k]);

        int found = 0;
        for(int jp = 0; jp < npartner[j]; jp++)
            if(partner[j][jp] == tag[i]) found = 1;
        if (!found) error->all(FLERR,"Internal fix bond/create error");

        // apply probability constraint
        // MIN,MAX insures values are added in same order on different procs

        if (fraction < 1.0) {
          min = MIN(probability[i],probability[j]);
          max = MAX(probability[i],probability[j]);
          if (0.5*(min+max) >= fraction) continue;
        }

        /*NL*///fprintf(screen,"creating bond btw atoms %d and %d (i has now %d bonds) at step %d\n",i,j,num_bond[i]+1,update->ntimestep);

        // if newton_bond is set, only store with I or J
        // if not newton_bond, store bond with both I and J

        if (!newton_bond || tag[i] < tag[j]) {
          if (num_bond[i] == atom->bond_per_atom)
            error->one(FLERR,"New bond exceeded bonds per atom in fix bond/create");
          bond_type[i][num_bond[i]] = btype;
          bond_atom[i][num_bond[i]] = tag[j];

          //reset history
          for (int ih = 0; ih < n_bondhist; ih++) {
              bond_hist[i][num_bond[i]][ih] = 0.;
          }
          num_bond[i]++;
        }
        // increment bondcount, convert atom to new type if limit reached

        bondcount[i]++;
        if (type[i] == iatomtype) {
          if (bondcount[i] == imaxbond) type[i] = inewtype;
        } else {
          if (bondcount[i] == jmaxbond) type[i] = jnewtype;
        }

        // count the created bond once

        if (tag[i] < tag[j]) ncreate++;
    }
  }

  // tally stats

  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  createcounttotal += createcount;
  atom->nbonds += createcount;

  /*NL*/ if(createcount && comm->me == 0) fprintf(screen,"Created %d bonds at timestep "BIGINT_FORMAT"\n",createcount,update->ntimestep);

  // trigger reneighboring if any bonds were formed

  if (createcount) next_reneighbor = update->ntimestep;
}

inline bool FixBondCreateGran::already_bonded(int i,int j)
{
    for(int k = 0; k < atom->num_bond[i]; k++)
      if(atom->bond_atom[i][k] == atom->tag[j]) return true;
    return false;
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::post_integrate_respa(int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

int FixBondCreateGran::pack_comm(int n, int *list, double *buf,
                             int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;

  if (commflag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<int>(bondcount[j]);
    }
    return 1;

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<int>(npartner[j]);
      for(int k=0; k<newperts; k++) //NP communicate all slots, also the empty ones
        buf[m++] = static_cast<int>(partner[j][k]);
      buf[m++] = probability[j];
    }
    return newperts + 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 0) {
    for (i = first; i < last; i++)
      bondcount[i] = static_cast<int> (buf[m++]);

  } else {
    for (i = first; i < last; i++) {
      npartner[i] = static_cast<int> (buf[m++]);
      for(int k=0; k<newperts; k++)
        partner[i][k] = static_cast<int>(buf[m++]);
      probability[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixBondCreateGran::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (commflag == 0) {
    for (i = first; i < last; i++)
      buf[m++] = bondcount[i];
    return 1;

  } else {
    for (i = first; i < last; i++) {
      buf[m++] = static_cast<int>(npartner[i]);
      for(int k=0; k<newperts; k++)
        buf[m++] = static_cast<int>(partner[i][k]);
    }
    return newperts + 2;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondCreateGran::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  int flag = 0;
  int nnew;

  if (commflag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      bondcount[j] += static_cast<int> (buf[m++]);
    }
  } else {
    //NP add new bonds coming from other proc
    for (i = 0; i < n; i++) {
      j = list[i];
      nnew = static_cast<int> (buf[m++]);
      if(nnew+npartner[j] > newperts)
      {
          flag = 1;
          nnew = newperts - npartner[j];
      }
      for(int k = npartner[j]; k < npartner[j]+newperts; k++)
      {
          if(k >= npartner[j]+nnew) m++; //NP do not do anything if
          else partner[j][k] = static_cast<int> (buf[m++]);
      }
      npartner[i] += nnew;
    }
  }
  if(flag) error->warning(FLERR,"Could not generate all possible bonds");
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixBondCreateGran::grow_arrays(int nmax)
{
  bondcount = (int *)
    memory->srealloc(bondcount,nmax*sizeof(int),"bond/create:bondcount");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixBondCreateGran::copy_arrays(int i, int j)
{
  bondcount[j] = bondcount[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixBondCreateGran::pack_exchange(int i, double *buf)
{
  buf[0] = bondcount[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixBondCreateGran::unpack_exchange(int nlocal, double *buf)
{
  bondcount[nlocal] = static_cast<int> (buf[0]);
  return 1;
}

/* ---------------------------------------------------------------------- */

double FixBondCreateGran::compute_vector(int n)
{
  if (n == 1) return (double) createcount;
  return (double) createcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixBondCreateGran::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*2 * sizeof(int);
  bytes += newperts*nmax * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
