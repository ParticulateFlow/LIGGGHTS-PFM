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

#include <stdio.h>
#include <string.h>
#include "modify.h"
#include "style_compute.h"
#include "style_fix.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "compute.h"
#include "group.h"
#include "update.h"
#include "fix_property_atom.h"
#include "fix_property_global.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   //NP modified C.K.
   add a fix property
------------------------------------------------------------------------- */

FixPropertyGlobal* Modify::add_fix_property_global(int narg,char **arg,const char *caller)
{
    if(narg < 5) error->all(FLERR,"Not enough arguments to add a fix property");
    add_fix(narg,arg);
    return static_cast<FixPropertyGlobal*>(find_fix_property(arg[3],"property/global",arg[4],0,0,caller));
}

FixPropertyAtom* Modify::add_fix_property_atom(int narg,char **arg,const char *caller)
{
    if(narg < 5) error->all(FLERR,"Not enough arguments to add a fix property");
    add_fix(narg,arg);
    return static_cast<FixPropertyAtom*>(find_fix_property(arg[3],"property/atom",arg[4],0,0,caller));
}

/* ----------------------------------------------------------------------
   find a fix scalar transport equation
------------------------------------------------------------------------- */

FixScalarTransportEquation* Modify::find_fix_scalar_transport_equation_strict(const char *equation_id)
{

    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,"transportequation/scalar") == 0)
      {
          FixScalarTransportEquation *fix_ste = static_cast<FixScalarTransportEquation*>(fix[ifix]);
          if(fix_ste->match_equation_id(equation_id)) return fix_ste;
      }
    return NULL;
}

FixScalarTransportEquation* Modify::find_fix_scalar_transport_equation(const char *equation_id)
{

    for(int ifix = 0; ifix < nfix; ifix++)
      if(dynamic_cast<FixScalarTransportEquation*>(fix[ifix]))
      {
          FixScalarTransportEquation *fix_ste = static_cast<FixScalarTransportEquation*>(fix[ifix]);
          if(fix_ste->match_equation_id(equation_id)) return fix_ste;
      }
    return NULL;
}

/* ----------------------------------------------------------------------
   find a fix with the requested style
------------------------------------------------------------------------- */

Fix* Modify::find_fix_style_strict(const char *style, int rank)
{
    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,style) == 0)
      {
          if(rank > 0) rank --;
          else return fix[ifix];
      }
    return NULL;
}

Fix* Modify::find_fix_style(const char *style, int rank)
{
    int len = strlen(style);
    for(int ifix = 0; ifix < nfix; ifix++)
      if(strncmp(fix[ifix]->style,style,len) == 0)
      {
          if(rank > 0) rank --;
          else return fix[ifix];
      }
    return NULL;
}

/* ----------------------------------------------------------------------
   find a compute with the requested style
------------------------------------------------------------------------- */

Compute* Modify::find_compute_style_strict(const char *style, int rank)
{
    for(int icompute = 0; icompute < ncompute; icompute++)
      if(strcmp(compute[icompute]->style,style) == 0)
      {
          if(rank > 0) rank --;
          else return compute[icompute];
      }
    return NULL;
}


/* ----------------------------------------------------------------------
   find a fix by ID, return NULL if not found
------------------------------------------------------------------------- */

Fix* Modify::find_fix_id(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) return NULL;
  return fix[ifix];
}

/* ----------------------------------------------------------------------
   find a compute by ID, return NULL if not found
------------------------------------------------------------------------- */

Compute* Modify::find_compute_id(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,compute[icompute]->id) == 0) break;
  if (icompute == ncompute) return NULL;
  return compute[icompute];
}

/* ----------------------------------------------------------------------
   find a fix by ID and style, return NULL if not found
------------------------------------------------------------------------- */

Fix* Modify::find_fix_id_style(const char *id,const char* style)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0 && strncmp(style,fix[ifix]->style,strlen(style)) == 0) break;
  if (ifix == nfix) return NULL;
  return fix[ifix];
}

/* ----------------------------------------------------------------------
   return number of fixes/computes with requested style
------------------------------------------------------------------------- */

int Modify::n_fixes_style(const char *style)
{
    int n_fixes,len;

    n_fixes = 0;
    len = strlen(style);

    for(int ifix = 0; ifix < nfix; ifix++)
      if(strncmp(fix[ifix]->style,style,len) == 0)
          n_fixes++;

    return n_fixes;
}

int Modify::n_computes_style(const char *style)
{
    int n_computes,len;

    n_computes = 0;
    len = strlen(style);

    for(int icompute = 0; icompute < ncompute; icompute++)
      if(strncmp(compute[icompute]->style,style,len) == 0)
          n_computes++;

    return n_computes;
}

int Modify::n_fixes_style_strict(const char *style)
{
    int n_fixes = 0;

    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,style) == 0)
          n_fixes++;

    return n_fixes;
}

/* ----------------------------------------------------------------------
   //NP modified C.K.

   checks if I am the first fix of a specified style
------------------------------------------------------------------------- */

bool Modify::i_am_first_of_style(Fix *fix_to_check)
{
    for(int ifix = 0; ifix < nfix; ifix++)
    {
        if(strcmp(fix[ifix]->style,fix_to_check->style) == 0)
        {
            if(fix_to_check == fix[ifix]) return true;
            return false;
        }
    }

    // no fix at all of this style found
    return false;
}

/* ----------------------------------------------------------------------
   //NP modified C.K.

   returns the index of first/last fix of specified style
------------------------------------------------------------------------- */

int Modify::index_first_fix_of_style(const char *style)
{
    for(int ifix = 0; ifix < nfix; ifix++)
        if(strcmp(fix[ifix]->style,style) == 0)
            return ifix;

    return -1;
}

int Modify::index_last_fix_of_style(const char *style)
{
    int idx = -1;

    for(int ifix = 0; ifix < nfix; ifix++)
        if(strcmp(fix[ifix]->style,style) == 0)
            idx = ifix;

    return idx;
}

/* ----------------------------------------------------------------------
   //NP modified C.K.

   returns the index of first fix that implements a specified function
------------------------------------------------------------------------- */

int Modify::my_index(Fix *fixptr)
{

    for(int ifix = 0; ifix < nfix; ifix++)
      if(fix[ifix] == fixptr)
        return ifix;

    return -1;
}

/* ----------------------------------------------------------------------
   //NP modified C.K.

   returns the index of first fix that implements a specified function
------------------------------------------------------------------------- */

int Modify::index_first_fix_with_function(const int FUNCTION, bool integrate)
{

    for(int ifix = 0; ifix < nfix; ifix++)
      if((!integrate || fix[ifix]->time_integrate) && (fmask[ifix] & FUNCTION))
        return ifix;

    return -1;
}

/* ----------------------------------------------------------------------
   //NP modified C.K.

   find a property registered by a fix property/global or fix property/atom
   check if it is of desired style
   return the index in the fix array
------------------------------------------------------------------------- */

Fix* Modify::find_fix_property(const char *varname,const char *style,const char *svmstyle,int len1,int len2,const char *caller,bool errflag)
{
  int ifix;
  char errmsg[500];
  Fix *fix_i = NULL;

  if((strcmp(style,"property/global")) && (strncmp(style,"property/atom",13)))
    error->all(FLERR,"Valid styles for find_fix_property are 'property/global' and 'property/atom'");

  if((len1 < 0) || (len2 < 0))
    error->all(FLERR,"Lengths for find_fix_property not valid");

  if(strncmp(style,"property/atom",13) == 0)
  {
    for(ifix = 0; ifix < nfix; ifix++)
    {
      if(strncmp(fix[ifix]->style,"property/atom",13) == 0 && dynamic_cast<FixPropertyAtom*>(fix[ifix]))
      {
        fix_i = static_cast<FixPropertyAtom*>(fix[ifix])->check_fix(varname,svmstyle,len1,len2,caller,errflag);
      }
      // check_fix returns either this or NULL
      if(fix_i) return fix_i;
    }
  }
  else if(strcmp(style,"property/global") == 0)
  {
    for(ifix = 0; ifix < nfix; ifix++)
    {
      if(strcmp(fix[ifix]->style,"property/global") == 0 && dynamic_cast<FixPropertyGlobal*>(fix[ifix]))
      {
        fix_i = static_cast<FixPropertyGlobal*>(fix[ifix])->check_fix(varname,svmstyle,len1,len2,caller,errflag);
      }
      // check_fix returns either this or NULL
      if(fix_i) return fix_i;
    }
  }

  // no fix found
  if(errflag)
  {
      sprintf(errmsg,"Could not locate a fix/property storing value(s) for %s as requested by %s.",varname,caller);
      error->all(FLERR,errmsg);
  }
  return NULL;
}


/* ----------------------------------------------------------------------
   return if fix restart in progress
------------------------------------------------------------------------- */

int Modify::fix_restart_in_progress()
{
    return  nfix_restart_global || nfix_restart_peratom;
}

/* ----------------------------------------------------------------------
   check if restart data available for this fix
------------------------------------------------------------------------- */

bool Modify::have_restart_data(Fix *f)
{
  //NP from Modify::add_fix

  // check if Fix is in restart_global list

  for (int i = 0; i < nfix_restart_global; i++)
    if (strcmp(id_restart_global[i],f->id) == 0 && strcmp(style_restart_global[i],f->style) == 0)
      return true;

  // check if Fix is in restart_peratom list

  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strcmp(id_restart_peratom[i],f->id) == 0 &&        strcmp(style_restart_peratom[i],f->style) == 0)
      return true;

  return false;
}

bool Modify::have_restart_data_style(const char* _style)
{
  // check if Fix is in restart_global list

  for (int i = 0; i < nfix_restart_global; i++)
    if (strncmp(style_restart_global[i],_style,strlen(_style)) == 0)
      return true;

  // check if Fix is in restart_peratom list

  for (int i = 0; i < nfix_restart_peratom; i++)
    if (strncmp(style_restart_peratom[i],_style,strlen(_style)) == 0)
      return true;

  return false;
}

int Modify::n_restart_data_global_style(const char* _style)
{
  int nhave = 0;

  // check if Fix is in restart_global list

  for (int i = 0; i < nfix_restart_global; i++)
    if (strncmp(style_restart_global[i],_style,strlen(_style)) == 0)
      nhave++;

  return nhave;
}

char* Modify::id_restart_data_global_style(const char* _style,int _rank)
{
  // check if Fix is in restart_global list

  for (int i = 0; i < nfix_restart_global; i++)
    if (strncmp(style_restart_global[i],_style,strlen(_style)) == 0)
    {
          if(_rank > 0) _rank --;
          else return id_restart_global[i];
    }

  return 0;
}

/* ----------------------------------------------------------------------
   let fixes extend bounding box
------------------------------------------------------------------------- */

void Modify::box_extent(double &xlo,double &xhi,double &ylo,double &yhi,double &zlo,double &zhi)
{
  //NP from Modify::add_fix

  // check if Fix is in restart_global list

  for (int i = 0; i < nfix; i++)
    fix[i]->box_extent(xlo,xhi,ylo,yhi,zlo,zhi);
}

/* ----------------------------------------------------------------------
   return min particle radius
------------------------------------------------------------------------- */

void Modify::max_min_rad(double &maxrad,double &minrad)
{
    maxrad = 0.;
    minrad = 1000.;
    int nlocal = atom->nlocal;
    double *radius = atom->radius;
    int ntypes = atom->ntypes;

    for (int i = 0; i < nfix; i++) {
      for (int j = 1; j <= ntypes; j++) {
        /*NL*/ //if (screen) fprintf(screen,"fix %s ntypes %d mr %f\n",modify->fix[i]->style,ntypes,modify->fix[i]->max_rad(j));
        maxrad = MathExtraLiggghts::max(maxrad,fix[i]->max_rad(j));
        if(modify->fix[i]->min_rad(j) > 0.)
            minrad = MathExtraLiggghts::min(minrad,fix[i]->min_rad(j));
      }
    }

    if (radius) {
      for (int i = 0; i < nlocal; i++) {
        maxrad = MathExtraLiggghts::max(maxrad,radius[i]);
        minrad = MathExtraLiggghts::min(minrad,radius[i]);
      }
    }

    MPI_Min_Scalar(minrad,world);
    MPI_Max_Scalar(maxrad,world);
}
