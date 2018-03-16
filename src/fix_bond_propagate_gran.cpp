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
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include "fix_bond_propagate_gran.h"
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

using namespace LAMMPS_NS;
using namespace FixConst;

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixBondPropagateGran::FixBondPropagateGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    restart_global = 1;
}

/* ---------------------------------------------------------------------- */

FixBondPropagateGran::~FixBondPropagateGran()
{
}

/* ---------------------------------------------------------------------- */

int FixBondPropagateGran::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */
//NP  copy partner info from neighbor lists to atom arrays
//NP  so can be exchanged with atoms
//NP  IMPORTANT: newtonflag = 1 assumed for all history values
void FixBondPropagateGran::pre_exchange()
{
  int i1,i2,ip,n;

  int **bondlist = neighbor->bondlist;
  double **bondhistlist = neighbor->bondhistlist;
  int nbondlist = neighbor->nbondlist;

  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  double ***bond_hist = atom->bond_hist;
  int n_bondhist = atom->n_bondhist;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;

  int newton_bond = force->newton_bond;

  //NP task 1
  //NP propagate bond contact history
  if(!n_bondhist) return;

  /*NL*/ //if(screen) fprintf(screen,"step " BIGINT_FORMAT ": FixBondPropagateGran::pre_exchange() (id %s)\n",update->ntimestep,id);
  /*NL*/ //if(screen) fprintf(screen,"nbondlist_propagate %d\n",nbondlist);

  for (n = 0; n < nbondlist; n++) {

    if(bondlist[n][3]) continue; //do not copy broken bonds

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    /*NL*///error->one(FLERR,"bond op...");

    if (newton_bond || i1 < nlocal)
    {
        ip = -1;
        for(int k = 0; k < num_bond[i1]; k++)
            if(bond_atom[i1][k] == tag[i2])
            {
                ip = k;
                break;
            }

        if(ip == -1)
        {
            /*NL*/ if(screen) fprintf(screen,"step " BIGINT_FORMAT ": tags %d %d \n",
            /*NL*/         update->ntimestep,atom->tag[i1],atom->tag[i2]);
            error->one(FLERR,"Failed to operate on granular bond history during copy i1");
        }

        for(int k = 0; k < n_bondhist; k++)
           bond_hist[i1][ip][k] = bondhistlist[n][k];
    }

    if (newton_bond || i2 < nlocal)
    {
        ip = -1;
        for(int k = 0; k < num_bond[i2]; k++)
            if(bond_atom[i2][k] == tag[i1])
            {
                ip = k;
                break;
            }

        if(ip == -1)
        {
            /*NL*/ if(screen) fprintf(screen,"step " BIGINT_FORMAT ": tags %d %d \n",
            /*NL*/         update->ntimestep,atom->tag[i1],atom->tag[i2]);
            error->one(FLERR,"Failed to operate on granular bond history during copy i2");
        }

        /*NL*///error->one(FLERR,"bond op");
        for(int k = 0; k < n_bondhist; k++)
           bond_hist[i2][ip][k] = -bondhistlist[n][k];
    }
  }

  //NP task 2
  //NP remove broken bonds
  //NP should be done equally on all processors

  for (n = 0; n < nbondlist; n++) {
    if(!bondlist[n][3]) continue; // continue if not broken

    i1 = bondlist[n][0];
    i2 = bondlist[n][1];

    //NP if the bond is broken, we remove it from
    //NP both atom data

    // delete bond from atom I if I stores it
    // atom J will also do this

    if (newton_bond || i1 < nlocal)
    {
        ip = -1;
        for(int k = 0; k < num_bond[i1]; k++)
            if(bond_atom[i1][k] == tag[i2])
            {
                ip = k;
                break;
            }

        if(ip == -1)
        {
            /*NL*/ if(screen) fprintf(screen,"step " BIGINT_FORMAT ": tags %d %d \n",
            /*NL*/         update->ntimestep,atom->tag[i1],atom->tag[i2]);
            error->one(FLERR,"Failed to operate on granular bond history during deletion1");
        }

        remove_bond(i1,ip);
    }

    if (newton_bond || i2 < nlocal)
    {
        ip = -1;
        for(int k = 0; k < num_bond[i2]; k++)
            if(bond_atom[i2][k] == tag[i1])
            {
                ip = k;
                break;
            }

        if(ip == -1)
        {
            /*NL*/ if(screen) fprintf(screen,"step " BIGINT_FORMAT ": tags %d %d \n",
            /*NL*/         update->ntimestep,atom->tag[i1],atom->tag[i2]);
            error->one(FLERR,"Failed to operate on granular bond history during deletion2");
        }

        remove_bond(i2,ip);
    }
  }
}

inline void FixBondPropagateGran::remove_bond(int ilocal,int ibond)
{
    /*NL*///if(screen) fprintf(screen,"step " BIGINT_FORMAT ": atom tag %d removing bond with atom tag %d\n",
    /*NL*///        update->ntimestep,atom->tag[ilocal],atom->bond_atom[ilocal][ibond]);
    /*NL*///error->one(FLERR,"removing bond");
    int nbond = atom->num_bond[ilocal];
    atom->bond_atom[ilocal][ibond] = atom->bond_atom[ilocal][nbond-1];
    atom->bond_type[ilocal][ibond] = atom->bond_type[ilocal][nbond-1];

    for(int k = 0; k < atom->n_bondhist; k++)
           atom->bond_hist[ilocal][ibond][k] = atom->bond_hist[ilocal][nbond-1][k];

    atom->num_bond[ilocal]--;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixBondPropagateGran::write_restart(FILE *fp)
{

  //NP write a dummy value
  int n = 0;
  double list[1];
  list[n++] = 1.;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  //NP write data to atom arrays where it can then be stored
  //NP can be done this way bc modify writes before avec
  pre_exchange();
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixBondPropagateGran::restart(char *buf)
{

}
