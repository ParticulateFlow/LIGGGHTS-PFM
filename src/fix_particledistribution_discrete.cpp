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

#include <stdlib.h>
#include <assert.h>
#include "atom.h"
#include "neighbor.h"
#include "update.h"
#include "fix_template_sphere.h"
#include "fix_particledistribution_discrete.h"
#include "modify.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "comm.h"
#include "force.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::FixParticledistributionDiscrete(LAMMPS *lmp, int narg, char **arg) :
  FixParticledistribution(lmp, narg, arg)
{
  restart_global = 1;

  // random number generator, same for all procs

  if (narg < 7)
    error->all(FLERR,"Illegal fix particledistribution/discrete command, not enough arguments");
  seed = atoi(arg[3]) + comm->me;
  random = new RanPark(lmp,seed);
  ntemplates = atoi(arg[4]);
  if(ntemplates < 1)
    error->all(FLERR,"Illegal fix particledistribution/discrete command, illegal number of templates");

  templates = new FixTemplateSphere*[ntemplates];
  distweight = new double[ntemplates];
  cumweight = new double[ntemplates];
  parttogen = new int[ntemplates];
  distorder = new int[ntemplates];

  int iarg = 5;

  int itemp=0;

  if(narg != iarg+2*ntemplates)
    error->all(FLERR,"Illegal fix particledistribution/discrete command, # of templates does not match # of arguments");

  // parse further args
  do {
    if(itemp == ntemplates) break;
    if(narg < iarg+1)
        error->all(FLERR,"Illegal fix particledistribution/discrete command, not enough arguments");
    int ifix = modify->find_fix(arg[iarg]);

    if(ifix < 0)
        error->all(FLERR,"Illegal fix particledistribution/discrete command, invalid ID for fix particletemplate provided");

    if(strncmp(modify->fix[ifix]->style,"particletemplate/",16))
        error->all(FLERR,"Illegal fix particledistribution/discrete command, fix is not of type particletemplate");

    templates[itemp] = static_cast<FixTemplateSphere*>(modify->fix[ifix]);
    distweight[itemp] = atof(arg[iarg+1]);
    if (distweight[itemp] < 0) error->all(FLERR,"Illegal fix particledistribution/discrete command, invalid weight");
    itemp++;
    iarg += 2;
  } while (iarg < narg);

  // check for double use of template which is not allowed
  for(int i = 0; i < ntemplates; i++)
      for(int j = 0; j < i; j++)
        if(templates[i] == templates[j])
            error->all(FLERR,"Illegal fix particledistribution/discrete command, cannot use the same template twice");

  // normalize distribution
  double weightsum = 0;
  for(int i = 0; i < ntemplates; i++)
    weightsum += distweight[i];

  if(comm->me == 0 && fabs(weightsum-1.) > 0.00001)
    error->warning(FLERR,"particledistribution/discrete: sum of distribution weights != 1, normalizing distribution");

  for(int i = 0; i < ntemplates; i++)
    distweight[i]/=weightsum;

  if(comm->me == 0 && screen)
  {
      fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on mass%%:\n",this->id);
      for(int i = 0; i < ntemplates; i++)
        fprintf(screen,"    %s: d=%e (max. bounding sphere) mass%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
  }

  // convert distribution from mass% to number%
  for(int i=0;i<ntemplates; i++)
    distweight[i]=distweight[i]/templates[i]->massexpect();

  weightsum=0;
  for(int i=0;i<ntemplates; i++)
    weightsum+=distweight[i];

  for(int i=0;i<ntemplates; i++)
    distweight[i]/=weightsum;

  if(comm->me == 0 && screen)
  {
      fprintf(screen,"Fix particledistribution/discrete (id %s): distribution based on number%%:\n",this->id);
      for(int i = 0; i < ntemplates; i++)
        fprintf(screen,"    %s: d=%e (max. bounding sphere) number%%=%f%%\n",templates[i]->id,2.*templates[i]->max_r_bound(),100.*distweight[i]);
  }

  //NP calculate cumulative distribution
  cumweight[0] = distweight[0];
  for(int i = 1; i < ntemplates; i++)
    cumweight[i] = distweight[i]+cumweight[i-1];

  //NP calculate mass and volume expectancy
  volexpect=0.;massexpect=0.;

  for(int i = 0; i < ntemplates; i++)
  {
      volexpect  += templates[i]->volexpect()  * distweight[i];
      massexpect += templates[i]->massexpect() * distweight[i];
  }

  //get min/maxtype
  maxtype = 0;
  mintype = 10000;
  for(int i = 0; i < ntemplates; i++)
  {
    if(templates[i]->maxtype() > maxtype)
      maxtype = templates[i]->maxtype();
    if(templates[i]->mintype() < mintype)
      mintype = templates[i]->mintype();
  }

  // check which template has the most spheres
  maxnspheres = 0;
  for(int i = 0; i < ntemplates;i++)
    if(templates[i]->number_spheres() > maxnspheres)
      maxnspheres=templates[i]->number_spheres();

  // sort the distributions by insertion volume (in descending order)
  // use bubble sort
  for(int i = 0; i < ntemplates; i++)
    distorder[i]=i;

  bool swaped;
  int n = ntemplates;
  do
  {
      swaped = false;
      for(int i = 0; i < ntemplates-1; i++)
      {
          if(templates[distorder[i]]->volexpect() < templates[distorder[i+1]]->volexpect())
          {
            //swap
            int tmp = distorder[i];
            distorder[i] = distorder[i+1];
            distorder[i+1] = tmp;
            swaped = true;
          }
      }
      n--;
  } while(swaped && n > 0);

  /*NL*/ //if(comm->me == 0 && screen) {
  /*NL*/ //    fprintf(screen,"particledistribution/discrete (id %s): sorted distributions based on volume extectancy. beginning from large, the rank is as follows:\n",this->id);
  /*NL*/ //    for(int i=0;i<ntemplates;i++)
  /*NL*/ //      fprintf(screen,"    %s: d=%e (bounding sphere) number%%=%f%%\n",templates[distorder[i]]->id,2.*templates[distorder[i]]->pti->r_bound,100.*distweight[distorder[i]]);
  /*NL*/ //}

  //calc max radius and bounding sphere radius

  maxrad = maxrbound = 0.;
  minrad = 1000.;

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->max_r_bound() > maxrbound)
        maxrbound = templates[i]->max_r_bound();

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->max_rad() > maxrad)
        maxrad = templates[i]->max_rad();

  for(int i = 0; i < ntemplates;i++)
      if(templates[i]->min_rad() < minrad)
        minrad = templates[i]->min_rad();

}

/* ---------------------------------------------------------------------- */

FixParticledistributionDiscrete::~FixParticledistributionDiscrete()
{
    delete []templates;
    delete []distweight;
    delete []cumweight;
    delete []parttogen;
    delete []distorder;
}

/* ----------------------------------------------------------------------
   prepares the fix for a series of randomize_list() command
   also prepares templates
       - deletes their old lists if present and allocates new lists
   typically only called once before first insertion step

   allocates for max # particles

   can be called by multiple fix insert commands, so check first if max #
   particles to be inserted is exceeded and only re-allocate in this case
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::random_init_list(int ntotal)
{
    int parttogen_max_i, n_pti_max_requested;
    int nprocs = comm->nprocs;

    //NP 2 is worst case of random->uniform, may be applied by each template
    ntotal += 2 * ntemplates;

    // number of requested pti
    n_pti_max_requested = 0;

    //NP calc parttogen_max for each template - 1.01 stands for worst case of random->uniform
    //NP with exact_number = 1, all the round-off remainders may be attributed to one particle template
    //NP round-off comes from FixInsert::distribute_ninsert_this() and FixParticledistributionDiscrete::randomize_list()
    for(int i = 0; i < ntemplates; i++)
    {
        parttogen_max_i = static_cast<int>(static_cast<double>(ntotal) * distweight[i] + static_cast<double>(1.01)*(ntemplates+nprocs));
        n_pti_max_requested += parttogen_max_i;

        // re-allocated if need more ptis in this template than allocated so far
        if(parttogen_max_i > templates[i]->n_pti_max)
        {
            templates[i]->delete_ptilist();
            templates[i]->init_ptilist(parttogen_max_i);
        }
    }

    // re-allocate if need more total ptis in distribution than allocated so far
    pti_list.reserve(n_pti_max_requested);
}

/* ----------------------------------------------------------------------
   tell all templates to generate their pti_list, wire their pti_list to
   the list in this fix. returns number of particles to be inserted.
   typically called once per insertion step

   for exact_number = 1, truncate distribution so to exactly meet
                               requested # particles
   for exact_number = 0, use random gen to fulfil distribution
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::randomize_list(int ntotal,int insert_groupbit,int exact_number)
{
    ninsert = ntotal;
    ninserted = 0;

    // use random generator so long-time average of insertion will represent distribution correctly
    if(exact_number == 0)
    {

        for(int i = 0; i < ntemplates; i++)
        {
           parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i] + random->uniform());
           /*NL*///if (screen) fprintf(screen,"parttogen[%d] = %d, ninsert %d \n",i,parttogen[i],ninsert);
        }
    }
    // truncate distribution so # particles to insert is met exactly
    else
    {
        int ninsert_truncated = 0, j;
        double *remainder = new double[ntemplates], rsum, r;

        // distribute particles and calculate remainder
        for(int i = 0; i < ntemplates; i++)
        {
           parttogen[i] = static_cast<int>(static_cast<double>(ninsert) * distweight[i]);
           ninsert_truncated += parttogen[i];
           remainder[i] = static_cast<double>(ninsert) * distweight[i] - static_cast<double>(parttogen[i]);
           /*NL*/ //if (screen) fprintf(screen,"parttogen[i] %d remainder[i] %f\n",parttogen[i],remainder[i]);
        }

        int ninsert_gap = ninsert - ninsert_truncated;

        //NP check if sum of remainders is equal to gap
        /*NL*/ rsum = 0.;
        /*NL*/ for(int i = 0; i < ntemplates; i++) rsum+=remainder[i];
        /*NL*/ //if (screen) fprintf(screen,"ninsert_gap %d, remaindersum %f - SHOULD BE EQUAL\n",ninsert_gap,rsum);

        // distribute remaining ninsert_gap particles
        for(int i = 0; i < ninsert_gap; i++)
        {
            r = random->uniform() * static_cast<double>(ninsert_gap);
            j = 0;
            rsum = remainder[0];

            while(j < (ntemplates-1) && rsum < r)
            {
                rsum += remainder[j];
                j++;
            }
            /*NL*/ //if (screen) fprintf(screen,"chose template %d\n",j);
            parttogen[j]++;
        }

        delete []remainder;
    }

    /*NL*/// MPI_Barrier(world); if (screen) fprintf(screen,"FixParticledistributionDiscrete::randomize_list on proc %d\n",comm->me);

    /*NL*/ //if(comm->me == 0 && screen) { fprintf(screen,"randomizing particles to generate out of the distribution\n");
    /*NL*/ //for(int i=0;i<ntemplates;i++)
    /*NL*/ //    fprintf(screen,"dist %d: statistically %f particles should be generated, chose to generate %d particles\n",
    /*NL*/ //    i,static_cast<double>(ninsert)*distweight[i],parttogen[i]); }

    // count total particle number to be inserted, let templates generate a pti_list
    ninsert = 0;
    for(int i = 0; i < ntemplates; i++)
    {
        ninsert += parttogen[i];
        templates[i]->randomize_ptilist(parttogen[i],groupbit | insert_groupbit);
    }

    // wire lists, make sure in correct order (large to small particles)

    for(int i = 0; i < ntemplates; i++)
    {
        int chosendist = distorder[i];
        for (int j = 0; j < parttogen[chosendist]; j++)
        {
            pti_list.push_back(templates[chosendist]->pti_list[j]);
        }
    }

    if(pti_list.size() != static_cast<pti_list_type::size_type>(ninsert))
        error->all(FLERR,"Internal error in FixParticledistributionDiscrete::randomize_list");

    ninserted = ninsert;
    return ninsert;
}

/* ----------------------------------------------------------------------
   preparations before insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::pre_insert(int n, FixPropertyAtom *fp, double val, int idx, int ival, int iidx)
{
    FixParticledistribution::pre_insert();

    // set fix property as desired by fix insert
    // loop to n, not n_pti
    if(fp)
    {
        for(int i = 0; i < n; i++)
        {
            pti_list[i]->fix_properties.push_back(fp);
            std::vector<double> value(1,val);
            pti_list[i]->fix_property_values.push_back(value);
        }
    }
    else if(idx >= 0)
    {
        for(int i = 0; i < n; i++)
        {
            pti_list[i]->property_index = idx;
            pti_list[i]->fix_property_dvalue = val;
        }
    }
    else if(iidx >= 0)
    {
        for(int i = 0; i < n; i++)
        {
            pti_list[i]->property_iindex = iidx;
            pti_list[i]->fix_property_ivalue = ival;
        }
    }

    // save bond history of ghost particles from neighbor to atom and delete it in neighbor;
    // this needs to be done prior to insertion of new particles since that will destroy the
    // neighbor/ghost information;
    // neighbor operates exclusively on local atom indices, while atom stores tags of bond partners
    // bond history that does not involve ghost atoms is handled in FixBondPropagateGran::pre_exchange()
    if (atom->molecular && atom->n_bondhist)
    {
        int i1,i2,ip,n;

        int **bondlist = neighbor->bondlist;
        double **bondhistlist = neighbor->bondhistlist;
        int nbondlist = neighbor->nbondlist;

        int **bond_atom = atom->bond_atom;
        int *num_bond = atom->num_bond;
        double ***bond_hist = atom->bond_hist;
        int n_bondhist = atom->n_bondhist; // number of bond history values

        int nlocal = atom->nlocal;
        int *tag = atom->tag;

        int newton_bond = force->newton_bond;

        for (n = 0; n < nbondlist; n++) {

          if(bondlist[n][3]) continue; //do not copy broken bonds

          i1 = bondlist[n][0]; // local index
          i2 = bondlist[n][1]; // local index

          if ((newton_bond || i1 < nlocal) && i2 >= nlocal)
          {
              ip = -1;
              for(int k = 0; k < num_bond[i1]; k++)
              {
                  if(bond_atom[i1][k] == tag[i2])
                  {
                      ip = k;
                      break;
                  }
              }

              if(ip == -1)
              {
                  /*NL*/ if(screen) fprintf(screen,"[%d] step " BIGINT_FORMAT ": nlocal %d ilocal %d %d tags %d %d mol %d %d\n",
                  /*NL*/         comm->me,update->ntimestep,nlocal,i1,i2,atom->tag[i1],atom->tag[i2], atom->molecule[i1], atom->molecule[i2]);
                  error->one(FLERR,"Failed to operate on granular bond history during copy i1");
              }

              for(int k = 0; k < n_bondhist; k++)
                 bond_hist[i1][ip][k] = bondhistlist[n][k];
          }

          if ((newton_bond || i2 < nlocal) && i1 >= nlocal)
          {
              ip = -1;
              for(int k = 0; k < num_bond[i2]; k++)
              {
                  if(bond_atom[i2][k] == tag[i1])
                  {
                      ip = k;
                      break;
                  }
              }

              if(ip == -1)
              {
                  /*NL*/ if(screen) fprintf(screen,"[%d] step " BIGINT_FORMAT ": nlocal %d ilocal %d %d tags %d %d mol %d %d\n",
                  /*NL*/         comm->me,update->ntimestep,nlocal,i1,i2,atom->tag[i1],atom->tag[i2], atom->molecule[i1], atom->molecule[i2]);
                  error->one(FLERR,"Failed to operate on granular bond history during copy i2");
              }

              for(int k = 0; k < n_bondhist; k++)
                 bond_hist[i2][ip][k] = -bondhistlist[n][k];
          }
        }

        // delete neighbor bondlist involving ghosts
        for (n = nbondlist-1; n >=0; --n)
        {
            i1 = bondlist[n][0];
            i2 = bondlist[n][1];

            if(i1 < nlocal && i2 < nlocal) continue;

            if(n < nbondlist-1)
            {
                bondlist[n][0] = bondlist[nbondlist-1][0];
                bondlist[n][1] = bondlist[nbondlist-1][1];
                bondlist[n][2] = bondlist[nbondlist-1][2];
                bondlist[n][3] = bondlist[nbondlist-1][3];

                for(int k = 0; k < n_bondhist; k++)
                   bondhistlist[n][k] = bondhistlist[nbondlist-1][k];
            }
            --nbondlist;
        }

        neighbor->nbondlist = nbondlist;
    }
}

/* ----------------------------------------------------------------------
   set particle properties - only pti needs to know which properties to set
   loop to n, not pti_list.size(), since not all particles may have been inserted
------------------------------------------------------------------------- */

int FixParticledistributionDiscrete::insert(int n)
{
    int ninserted_spheres_local = 0;
    assert(static_cast<pti_list_type::size_type>(n) <= pti_list.size());
    for(int i = 0; i < n; i++)
    {
        /*NL*/ //if (screen) fprintf(screen,"inserting pti #%d\n",i);
        ninserted_spheres_local += pti_list[i]->insert();
    }
    return ninserted_spheres_local;
}

/* ----------------------------------------------------------------------
   wrap up insertion
------------------------------------------------------------------------- */

void FixParticledistributionDiscrete::finalize_insertion(bool pti_delete_flag)
{
  for(int i = 0; i < ntemplates; i++)
    templates[i]->finalize_insertion();

  if(pti_delete_flag){
    for(pti_list_type::iterator it = pti_list.begin();it != pti_list.end();it++){
      if(*it != 0)
        delete *it;
    }
  }

  pti_list.clear();
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::min_rad(int type) const
{
    //get minrad
    double minrad_type = 1000.;
    for(int i = 0; i < ntemplates;i++)
    {
      if(
            (type >= templates[i]->mintype() && type <= templates[i]->maxtype()) &&
            (templates[i]->min_rad() < minrad_type)
        )
        minrad_type = templates[i]->min_rad();
    }

    return minrad_type;
}

/* ----------------------------------------------------------------------*/

double FixParticledistributionDiscrete::max_rad(int type) const
{
    //get maxrad
    double maxrad_type = 0.;
    for(int i = 0; i < ntemplates;i++)
    {
      if(
          (type >= templates[i]->mintype() && type <= templates[i]->maxtype()) &&
          (templates[i]->max_rad() > maxrad_type)
        )
        maxrad_type = templates[i]->max_rad();
    }

    return maxrad_type;
}

/* ----------------------------------------------------------------------
   return true if any of the templates is of type multisphere
------------------------------------------------------------------------- */

bool FixParticledistributionDiscrete::has_multisphere()
{
  for(int i=0;i<ntemplates;i++){
    // strstr returns pointer to beginning of substring if substring is found
    if(strstr(templates[i]->style,"multisphere") != 0)
      return true;
  }
  return false;
}

/* ----------------------------------------------------------------------*/

ParticleToInsert* FixParticledistributionDiscrete::get_random_particle(int insert_groupbit)
{
  // first, choose distribution
  double r = random->uniform();
  int i_template=0;
  while(cumweight[i_template] < r)
    i_template++;

  // avoid floating point issues
  if(i_template == ntemplates)
    i_template = ntemplates-1;

  return templates[i_template]->get_single_random_pti(groupbit | insert_groupbit);
}
