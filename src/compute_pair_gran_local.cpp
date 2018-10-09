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
#include <string.h>
#include "compute_pair_gran_local.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "pair_gran.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "memory.h"
#include "error.h"
#include "pair_gran.h"
#include "domain.h"
#include "fix_heat_gran_conduction.h"
#include "fix_wall_gran.h"
#include "vector_liggghts.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::ComputePairGranLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all(FLERR,"Illegal compute pair/gran/local or wall/gran/local command");

  local_flag = 1;
  nmax = 0;
  array = NULL;

  // store everything by default expect heat flux
  posflag = velflag = idflag = fflag = tflag = hflag = aflag = 1;

  // do not store fn, ft, delta, heat flux by default
  fnflag = ftflag = deltaflag = hfflag = 0;

  // if further args, store only the properties that are listed
  if(narg > 3)
     posflag = velflag = idflag = fflag = tflag = hflag = aflag = 0;

  for (int iarg = 3; iarg < narg; iarg++)
  {
    //int i = iarg-3;
    if (strcmp(arg[iarg],"pos") == 0) posflag = 1;
    else if (strcmp(arg[iarg],"vel") == 0) velflag = 1;
    else if (strcmp(arg[iarg],"id") == 0) idflag = 1;
    else if (strcmp(arg[iarg],"force") == 0) fflag = 1;
    else if (strcmp(arg[iarg],"force_normal") == 0) fnflag = 1;
    else if (strcmp(arg[iarg],"force_tangential") == 0) ftflag = 1;
    else if (strcmp(arg[iarg],"torque") == 0) tflag = 1;
    else if (strcmp(arg[iarg],"history") == 0) hflag = 1;
    else if (strcmp(arg[iarg],"contactArea") == 0) aflag = 1;
    else if (strcmp(arg[iarg],"delta") == 0) deltaflag = 1;
    else if (strcmp(arg[iarg],"heatFlux") == 0) hfflag = 1;
    else error->compute_error(FLERR,this,"Invalid keyword");
  }

  // default: pair data
  wall = 0;

  reference_exists = 0;

  fixwall = NULL;
  fixheat = NULL;
  pairgran = NULL;

  if(update->ntimestep > 0 && !modify->fix_restart_in_progress())
    error->compute_error(FLERR,this,"Need to define this compute before first run");
}

/* ---------------------------------------------------------------------- */

ComputePairGranLocal::~ComputePairGranLocal()
{
  memory->destroy(array);

  if(reference_exists == 0) return;

  //NP TODO should do this in a pre_delete() step with uncomputeflag - as for fixes (unfixflag)

  if(wall == 0) pairgran->unregister_compute_pair_local(this);
  else fixwall->unregister_compute_wall_local(this);

  if(fixheat) fixheat->unregister_compute_pair_local(this);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::post_create()
{
  // check if wall data requested
  if(strcmp(style,"wall/gran/local") == 0) wall = 1;

  // initialize once as dump parses in constructor for length of per-atom data
  //NP take care to avoid two requests which would mess up neigh lists
  init_cpgl(false);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init()
{
    init_cpgl(false); /*NL*/ //HERE
    /*NL*/ //if (screen) fprintf(screen,"ComputePairGranLocal::init()\n");
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_cpgl(bool requestflag)
{
  int ifix, n_wall_fixes;
  FixWallGran *fwg;

  newton_pair = force->newton_pair;

  if (idflag && atom->tag_enable == 0)
      error->all(FLERR,"Compute pair/gran/local requested to compute IDs, this requires atoms have IDs.");

  // if available from previous run, unregister
  if(pairgran && reference_exists)
      pairgran->unregister_compute_pair_local(this);

  if(fixwall && reference_exists)
      fixwall->unregister_compute_wall_local(this);

  if(fixheat && reference_exists)
      fixheat->unregister_compute_pair_local(this);

  fixwall = NULL;
  fixheat = NULL;
  pairgran = NULL;

  //NP request neigh list for pair
  // register heat transfer fix for pair data
  if(wall == 0)
  {
      if (force->pair == NULL)
        error->all(FLERR,"No pair style is defined for compute pair/gran/local");

      pairgran = (PairGran*)force->pair_match("gran",0);

      if (pairgran == NULL)
        error->all(FLERR,"No valid granular pair style found for use with compute pair/gran/local");

      if (pairgran->cplenable() == 0)
        error->all(FLERR,"Pair style does not support compute pair/gran/local");

      //NP register with the pair style
      //NP is done twice before first run, out of post_create() + init()
      //NP no problem since always unregister when init_cpgl() called
      pairgran->register_compute_pair_local(this,dnum);

      if(requestflag)
      {
          //NP need an occasional gran neighbor list
          int irequest = neighbor->request((void *) this);
          neighbor->requests[irequest]->pair = 0;
          neighbor->requests[irequest]->compute = 1;
          neighbor->requests[irequest]->half = 0;
          neighbor->requests[irequest]->gran = 1;
          //neighbor->requests[irequest]->granhistory = 1;
          neighbor->requests[irequest]->occasional = 1;
      }

      // register heat transfer fix if applicable
      if(hfflag)
      {
          for(ifix = 0; ifix < modify->nfix; ifix++)
          {
              if(strcmp(modify->fix[ifix]->style,"heat/gran/conduction") == 0)
              {
                  fixheat = static_cast<FixHeatGranCond*>(modify->fix[ifix]);
              }
          }
          if(!fixheat) error->warning(FLERR,"Compute pair/gran/local can not calculate heat flux values since no fix heat/gran/conduction not compute them");

          // group of this compute and heat transfer fix must be same so same number of pairs is computed
          if(groupbit != fixheat->groupbit) error->all(FLERR,"Compute pair/gran/local group and fix heat/gran/conduction group cannot be different");
          fixheat->register_compute_pair_local(this);
      }
  }
  // register with granular wall, only accept mesh walls
  else
  {
      if(fixwall) fixwall->unregister_compute_wall_local(this);
      fixwall = NULL;

      n_wall_fixes = modify->n_fixes_style("wall/gran");

      for(ifix = 0; ifix < n_wall_fixes; ifix++)
      {
          fwg = static_cast<FixWallGran*>(modify->find_fix_style("wall/gran",ifix));
          if(fwg->is_mesh_wall()) fixwall = fwg;
      }

      if(!fixwall) error->all(FLERR,"Compute wall/gran/local requires a fix of type wall/gran using one or more mesh walls. This fix has come before the compute in the script");
      fixwall->register_compute_wall_local(this,dnum);
  }

  // at this point we know that ptr is valid
  reference_exists = 1;

  if(hflag && dnum == 0) error->all(FLERR,"Compute pair/gran/local or wall/gran/local can not calculate history values since pair or wall style does not compute them");
  // standard values: pos1,pos2,id1,id2,extra id for mesh wall,force,torque,contact area

  nvalues = posflag*6 + velflag*6 + idflag*3 + fflag*3 + fnflag*3 + ftflag*3 + tflag*3 + hflag*dnum + aflag + deltaflag + hfflag;
  size_local_cols = nvalues;

  /*NL*/ //if (screen) fprintf(screen,"nvalues = %d\n",nvalues);
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::reference_deleted()
{
    /*NL*/ //if (screen) fprintf(screen,"BLAHHH\n");
    reference_exists = 0;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::compute_local()
{
  invoked_local = update->ntimestep;

  if(!reference_exists) error->one(FLERR,"Compute pair/gran/local or wall/gran/local reference does no longer exist (pair or fix deleted)");

  // count local entries and compute pair info

  if(wall == 0) ncount = count_pairs();        // # pairs is ensured to be the same for pair and heat
  else          ncount = count_wallcontacts(); // # wall contacts ensured to be same for wall/gran and heat

  /*NL*/ //if (screen) fprintf(screen,"ncount %d\n",ncount);

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;

  /*NL*///if (screen) fprintf(screen,"ncount %d\n",ncount);
  /*NL*///if(ncount)error->all(FLERR,"ncount");

  // get pair data
  if(wall == 0)
  {
      ipair = 0;
      if(pairgran == NULL) error->one(FLERR,"null");
      pairgran->compute_pgl(0,0);

      // get heat flux data
      if(fixheat)
      {
          ipair = 0;
          fixheat->cpl_evaluate(this);
      }
  }
  // get wall data
  else
  {
      ipair = 0;
      // this also calls heat transfer if necessary
      fixwall->post_force_pgl();
  }
}

/* ----------------------------------------------------------------------
   count pairs on this proc
------------------------------------------------------------------------- */

int ComputePairGranLocal::count_pairs()
{
  int i,j,m,n,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  //NP invoke neighbor list (will copy or build if necessary)
  //NP cannot invoke here, since will copy outdated data from
  //NP fix contact history, so rather access pair
  //neighbor->build_one(list->index); /*NL*/ //HERE

  list = pairgran->list; /*NL*/ //HERE
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I or J are not in group

  m = n = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j >= nall) j %= nall;

      if (!(mask[j] & groupbit)) continue;
      if (newton_pair == 0 && j >= nlocal && atom->tag[i] <= atom->tag[j]) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= (radius[i]+radius[j])*(radius[i]+radius[j])) continue;
      m++;
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   count wall contacts on this proc
------------------------------------------------------------------------- */

int ComputePairGranLocal::count_wallcontacts()
{
    // account for fix group
    return fixwall->n_contacts_local(groupbit);
}

/* ----------------------------------------------------------------------
   add data from particle-particle contact on this proc
------------------------------------------------------------------------- */

void ComputePairGranLocal::add_pair(int i,int j,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist)
{
    double del[3],r,rsq,radi,radj,contactArea;
    double *xi,*xj,xi_w[3],xj_w[3],*vi,*vj;
    int nlocal;

    //NP return if pair not meant to be included (was not counted in the count() procedure
    if (!(atom->mask[i] & groupbit)) return;
    if (!(atom->mask[j] & groupbit)) return;

    nlocal = atom->nlocal;

    if (newton_pair == 0 && j >= nlocal && atom->tag[i] <= atom->tag[j]) return;

    xi = atom->x[i];
    xj = atom->x[j];
    vi = atom->v[i];
    vj = atom->v[j];

    //NP simply make a local copy of the data
    int n = 0;
    if(posflag)
    {
        vectorCopy3D(xi,xi_w);
        vectorCopy3D(xj,xj_w);
        domain->remap(xi_w);
        domain->remap(xj_w);
        vectorToBuf3D(xi_w,array[ipair],n);
        vectorToBuf3D(xj_w,array[ipair],n);
    }
    if(velflag)
    {
        vectorToBuf3D(vi,array[ipair],n);
        vectorToBuf3D(vj,array[ipair],n);
    }
    if(idflag)
    {
        array[ipair][n++] = static_cast<double>(atom->tag[i]);
        array[ipair][n++] = static_cast<double>(atom->tag[j]);
        if(i < nlocal && j < nlocal)
            array[ipair][n++] = 0.;
        else
        {
            if(domain->is_periodic_ghost(i) || domain->is_periodic_ghost(j))
                array[ipair][n++] = 1.;
            else
                array[ipair][n++] = 0.;
        }
    }
    if(fflag)
    {
        array[ipair][n++] = fx;
        array[ipair][n++] = fy;
        array[ipair][n++] = fz;
    }
    if(fnflag)
    {
        double fc[3],fn[3],normal[3];

        vectorConstruct3D(fc,fx,fy,fz);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);

        array[ipair][n++] = fn[0];
        array[ipair][n++] = fn[1];
        array[ipair][n++] = fn[2];
    }
    if(ftflag)
    {
        double fc[3],fn[3],ft[3],normal[3];

        vectorConstruct3D(fc,fx,fy,fz);
        vectorSubtract3D(atom->x[j],atom->x[i],normal);
        vectorNormalize3D(normal);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);
        vectorSubtract3D(fc,fn,ft);

        array[ipair][n++] = ft[0];
        array[ipair][n++] = ft[1];
        array[ipair][n++] = ft[2];
    }
    if(tflag)
    {
        array[ipair][n++] = tor1;
        array[ipair][n++] = tor2;
        array[ipair][n++] = tor3;
    }
    if(hflag)
    {
        for(int d = 0; d < dnum; d++)
           array[ipair][n++] = hist[d];
    }
    if(aflag)
    {
        radi = atom->radius[i];
        radj = atom->radius[j];
        vectorSubtract3D(atom->x[i],atom->x[j],del);
        rsq = vectorMag3DSquared(del);
        r = sqrt(rsq);
        contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/rsq;
        array[ipair][n++] = contactArea;
        /*NL*///if (screen) fprintf(screen,"adding %f rsq %f\n",contactArea,rsq);
    }
    if(deltaflag)
    {
        {
          radi = atom->radius[i];
          radj = atom->radius[j];
          vectorSubtract3D(atom->x[i],atom->x[j],del);
          array[ipair][n++] = radi+radj-vectorMag3D(del);

        }
    }

    ipair++;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::add_heat(int i,int j,double hf)
{
    if (newton_pair == 0 && j >= atom->nlocal && atom->tag[i] <= atom->tag[j]) return;

    //NP return if pair not meant to be included (was not counted in the count() procedure
    if (!(atom->mask[i] & groupbit)) return;
    if (!(atom->mask[j] & groupbit)) return;

    if(hfflag)
    {
        // heat flux is always last value
        array[ipair][nvalues-1] = hf;
        /*NL*///if (screen) fprintf(screen,"adding heat flux %f nvalues %d\n", hf,nvalues);
    }
    else error->one(FLERR,"Illegal situation in ComputePairGranLocal::add_heat");

    // inc counter again, since reset in compute_local() after getting pair data
    ipair++;
}

/* ----------------------------------------------------------------------
   add data from particle-wall contact on this proc
------------------------------------------------------------------------- */

void ComputePairGranLocal::add_wall_1(int iFMG,int idTri,int iP,double *contact_point,double *v_wall)
{
    if (!(atom->mask[iP] & groupbit)) return;

    //NP simply make a local copy of the data
    int n = 0;

    if(posflag)
    {
        array[ipair][n++] = contact_point[0];
        array[ipair][n++] = contact_point[1];
        array[ipair][n++] = contact_point[2];
        n += 3;
    }

    if(velflag)
    {
        if(v_wall)
        {
            array[ipair][n++] = v_wall[0];
            array[ipair][n++] = v_wall[1];
            array[ipair][n++] = v_wall[2];
        }
        else
        {
            array[ipair][n++] = 0.;
            array[ipair][n++] = 0.;
            array[ipair][n++] = 0.;
        }
        n += 3;
    }

    if(idflag)
    {
        array[ipair][n++] = static_cast<double>(iFMG);
        array[ipair][n++] = static_cast<double>(idTri);
        n += 1;
        /*NL*/ //if (screen) fprintf(screen,"having iFmG %d and idtri %d\n",iFMG,idTri);
    }

    if(fflag)
    {
        n += 3;
    }

    if(fnflag)
    {
        array[ipair][n++] = contact_point[0];
        array[ipair][n++] = contact_point[1];
        array[ipair][n++] = contact_point[2];
    }

    if(ftflag)
    {
        array[ipair][n++] = contact_point[0];
        array[ipair][n++] = contact_point[1];
        array[ipair][n++] = contact_point[2];
    }
    /*
    if(tflag)
    {
        n += 3;
    }
    if(hflag)
    {
        n += dnum;
    }
    if(aflag)
    {
        n++;
    }
    if(deltaflag)
    {
        n++;
    }*/
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::add_wall_2(int i,double fx,double fy,double fz,double tor1,double tor2,double tor3,double *hist,double rsq, double *normal)
{
    double contactArea;

    //NP return if pair not meant to be included (was not counted in the count() procedure
    if (!(atom->mask[i] & groupbit)) return;

    //NP simply make a local copy of the data
    int n = 0;

    if(posflag)
    {
        n += 3;
        array[ipair][n++] = atom->x[i][0];
        array[ipair][n++] = atom->x[i][1];
        array[ipair][n++] = atom->x[i][2];
    }
    if(velflag)
    {
        n += 3;
        array[ipair][n++] = atom->v[i][0];
        array[ipair][n++] = atom->v[i][1];
        array[ipair][n++] = atom->v[i][2];
    }
    if(idflag)
    {
        n += 2;
        array[ipair][n++] = static_cast<double>(atom->tag[i]);
        /*NL*/ //if (screen) fprintf(screen,"having tag %d\n",atom->tag[i]);
    }
    if(fflag)
    {
        array[ipair][n++] = fx;
        array[ipair][n++] = fy;
        array[ipair][n++] = fz;
    }
    if(fnflag)
    {
        double fc[3],fn[3];

        vectorConstruct3D(fc,fx,fy,fz);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);

        array[ipair][n++] = fn[0];
        array[ipair][n++] = fn[1];
        array[ipair][n++] = fn[2];
    }
    if(ftflag)
    {
        double fc[3],fn[3],ft[3];

        vectorConstruct3D(fc,fx,fy,fz);
        double fnmag = vectorDot3D(fc,normal);
        vectorScalarMult3D(normal,fnmag,fn);
        vectorSubtract3D(fc,fn,ft);

        array[ipair][n++] = ft[0];
        array[ipair][n++] = ft[1];
        array[ipair][n++] = ft[2];
    }
    if(tflag)
    {
        array[ipair][n++] = tor1;
        array[ipair][n++] = tor2;
        array[ipair][n++] = tor3;
    }
    if(hflag)
    {
        for(int d = 0; d < dnum; d++)
           array[ipair][n++] = hist[d];
    }
    if(aflag)
    {
        contactArea = (atom->radius[i]*atom->radius[i]-rsq)*M_PI;
        array[ipair][n++] = contactArea;
        /*NL*///if (screen) fprintf(screen,"adding %f rsq %f\n",contactArea,rsq);
    }
    if(deltaflag)
    {
        array[ipair][n++] = atom->radius[i]-sqrt(rsq);
    }

    // wall_1 and wall_2 are always called

    ipair++;
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::add_heat_wall(int ip,double hf)
{
    if (!(atom->mask[ip] & groupbit)) return;

    if(hfflag)
    {
        // heat flux is always last value
        // use ipair -1 , add_heat_wall is not always called
        array[ipair-1][nvalues-1] = hf;
    }
    else error->one(FLERR,"Illegal situation in ComputePairGranLocal::add_heat_wall");
}

/* ---------------------------------------------------------------------- */

void ComputePairGranLocal::reallocate(int n)
{
  // grow vector or array and indices array

  while (nmax < n) nmax += DELTA;

  memory->destroy(array);
  memory->create(array,nmax,nvalues,"pair/local:array");
  array_local = array;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePairGranLocal::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}
