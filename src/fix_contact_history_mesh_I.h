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

#ifndef LMP_CONTACT_HISTORY_MESH_I_H
#define LMP_CONTACT_HISTORY_MESH_I_H

  //NP coded in the header, so they can be inlined by other classes

/*NL*/ #define DEBUG_P_TAG 607

  /* ---------------------------------------------------------------------- */

  /*NL*/inline void FixContactHistoryMesh::debug(int iP, int idTri)
  /*NL*/{
  /*NL*/ if(screen && 10024 == update->ntimestep) {
  /*NL*/        int iDeb = atom->map(DEBUG_P_TAG);
  /*NL*/        int nn = static_cast<int>(round(fix_nneighs_->vector_atom[iDeb]));
  /*NL*/        fprintf(screen,"***contact particle id %d with tri ID %d at step " BIGINT_FORMAT ", nn %d iDeb %d\n",atom->tag[iP],idTri,update->ntimestep,nn,iDeb);
  /*NL*/        for(int kk = 0; kk < nn; kk++)
  /*NL*/           fprintf(screen,"       * neigh %d of part id %d: tri id %d\n",kk,DEBUG_P_TAG ,partner_[iDeb][kk]); }
  /*NL*/}

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::handleContact(int iP, int idTri, double *&history)
  {
    /*NL*/ //if(screen && DEBUG_P_TAG == atom->tag[iP]) fprintf(screen,"***contact with tri ID %d (index %d) at step " BIGINT_FORMAT "\n",idTri,iP,update->ntimestep);

    // check if contact with iTri was there before
    // if so, set history to correct location and return
    if(haveContact(iP,idTri,history))
      return true;

    /*NL*/// if (screen) fprintf(screen,"   new contact - adding\n");

    // else new contact - add contact if did not calculate contact with coplanar neighbor already
    //NP this can be detected via delflag == false
    if(coplanarContactAlready(iP,idTri))
        // did not add new contact
        return false;
    else
    {
        addNewTriContactToExistingParticle(iP,idTri,history);

        // check if one of the contacts of previous steps is coplanar with iTri
        //NP can be seen via delflag == false
        // if so, copy history
        // also check if this contact has delflag = false, i.e. has been executed already
        // this step. If so, signalize not to execute this contact (return false)
        checkCoplanarContactHistory(iP,idTri,history);
        return true;
    }
  }



  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::swap(int ilocal,int ineigh, int jneigh, bool keepflag_swap)
  {
      //NP swap data of ineigh and jneigh

      int id_temp;

      id_temp                  = partner_[ilocal][ineigh];
      partner_[ilocal][ineigh] = partner_[ilocal][jneigh];
      partner_[ilocal][jneigh] = id_temp;

      vectorCopyN(&(contacthistory_[ilocal][ineigh*dnum_]),swap_,                                   dnum_);
      vectorCopyN(&(contacthistory_[ilocal][jneigh*dnum_]),&(contacthistory_[ilocal][ineigh*dnum_]),dnum_);
      vectorCopyN(swap_,                                   &(contacthistory_[ilocal][jneigh*dnum_]),dnum_);

      if(keepflag_swap)
      {
          const bool keepflag_temp  = keepflag_[ilocal][ineigh];
          keepflag_[ilocal][ineigh] = keepflag_[ilocal][jneigh];
          keepflag_[ilocal][jneigh] = keepflag_temp;
      }
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::haveContact(int iP, int idTri, double *&history)
  {
    int *tri = partner_[iP];
    const int nneighs = fix_nneighs_->get_vector_atom_int(iP);

    /*NL*/ //if (screen) fprintf(screen,"nneighs %d dnum_ %d\n",nneighs,dnum_);

    for(int i = 0; i < nneighs; i++)
    {
        if(tri[i] == idTri)
        {
            if(dnum_ > 0) history = &(contacthistory_[iP][i*dnum_]);
            keepflag_[iP][i] = true;
            return true;
        }
    }
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline bool FixContactHistoryMesh::coplanarContactAlready(int iP, int idTri)
  {
    const int nneighs = fix_nneighs_->get_vector_atom_int(iP);
    for(int i = 0; i < nneighs; i++)
    {
      /*NL*/// if (screen) fprintf(screen," step %d: iP %d, partner %d: %d, coplanar %s delflag %s, mesh_->map(tri[i]) %d\n",
      /*NL*///                update->ntimestep,iP,i,tri[i],mesh_->areCoplanarNodeNeighs(tri[i],idTri)?"y":"n",delflag[iP][i]?"y":"n",mesh_->map(tri[i]));
      /*NL*/ //if (screen) fprintf(screen,"ip %d i %d idTri %d tri[i] %d\n",iP,i,idTri,tri[i]);
      //NP do only if old partner owned or ghost on this proc
      int idPartnerTri = partner_[iP][i];

      if(idPartnerTri >= 0 && idPartnerTri != idTri && mesh_->map(idPartnerTri) >= 0 && mesh_->areCoplanarNodeNeighs(idPartnerTri,idTri))
      {
        /*NL*/ //if (screen) fprintf(screen," step %d: idTri %d, delflag %s \n",update->ntimestep,idTri,delflag[iP][i]?"true":"false");

        // other coplanar contact handled already - do not handle this contact
        if(keepflag_[iP][i]) return true;
      }
    }

    // no coplanar contact found - handle this contact
    return false;
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::checkCoplanarContactHistory(int iP, int idTri, double *&history)
  {
    int *tri = partner_[iP];
    const int nneighs = fix_nneighs_->get_vector_atom_int(iP);

    for(int i = 0; i < nneighs; i++)
    {
      //NP do only if old partner owned or ghost on this proc
      if(tri[i] >= 0 && tri[i] != idTri && mesh_->map(tri[i]) >= 0 && mesh_->areCoplanarNodeNeighs(tri[i],idTri))
      {
          //NP this is illegal since coplanarContactAlready() should have avoided this
          /*NL*/ if(keepflag_[iP][i]) error->one(FLERR,"internal error");

          // copy contact history
          if(dnum_ > 0) vectorCopyN(&(contacthistory_[iP][i*dnum_]),history,dnum_);
          /*NL*/// if (screen) fprintf(screen,"Found coplanar contact, old contact hist %f %f %f\n",
          /*NL*///                   contacthistory[iP][i][0],contacthistory[iP][i][1],contacthistory[iP][i][2]);
          /*NL*/// if (screen) fprintf(screen,"Found coplanar contact, new contact hist %f %f %f\n",
          /*NL*///                   history[0],history[1],history[2]);
          /*NL*/ //error->one(FLERR,"end");
      }
    }
  }

  /* ---------------------------------------------------------------------- */

  inline void FixContactHistoryMesh::addNewTriContactToExistingParticle(int iP, int idTri, double *&history)
  {
      //NP error if # contacts larger than # of neighs
      //NP cannot store contact history in this case since allocation insufficient
      //NP should not happen since only neighbor can become contact
      const int nneighs = fix_nneighs_->get_vector_atom_int(iP);
      int iContact = -1;

      if(-1 == idTri)
        error->one(FLERR,"internal error");

      if(npartner_[iP] >= nneighs)
      {
        /*NL*/ //if (screen)
        /*NL*/ //{
        /*NL*/ //  fprintf(screen,"step " BIGINT_FORMAT ": for particle tag %d, proc %d npartner %d nneighs %d newcontact id %d\n",update->ntimestep,atom->tag[iP],comm->me,npartner_[iP],nneighs,idTri);
        /*NL*/ //  for(int kk = 0; kk < nneighs; kk++)
        /*NL*/ //    fprintf(screen,"  neigh %d: tri id %d\n",kk,partner_[iP][kk]);
        /*NL*/ //}
        error->one(FLERR,"internal error");
      }

      //NP look for first free storage position in array
      for(int ineigh = 0; ineigh < nneighs; ineigh++)
      {
          if(-1 == partner_[iP][ineigh])
          {
              iContact = ineigh;
              break;
          }
      }

      if(iContact >= nneighs)
        error->one(FLERR,"internal error");

      partner_[iP][iContact] = idTri;
      keepflag_[iP][iContact] = true;

      if(dnum_ > 0)
      {
          history = &(contacthistory_[iP][iContact*dnum_]);
          vectorZeroizeN(history,dnum_);
      }
      else
          history = 0;

      npartner_[iP]++;

      /*NL*/ //if(screen && DEBUG_P_TAG == atom->tag[iP]) {
      /*NL*/ //   fprintf(screen,"step " BIGINT_FORMAT " adding contact # %d: tri id %d, npartner %d nneighs %d\n",update->ntimestep,npartner_[iP],idTri,npartner_[iP],nneighs);
      /*NL*/ //       for(int kk = 0; kk < nneighs; kk++)
      /*NL*/ //         fprintf(screen,"    neigh %d: tri id %d\n",kk,partner_[iP][kk]); }
  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistoryMesh::n_contacts()
  {
    int ncontacts = 0, nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
           ncontacts += npartner_[i];
    return ncontacts;
  }

  /* ---------------------------------------------------------------------- */

  inline int FixContactHistoryMesh::n_contacts(int contact_groupbit)
  {
    int ncontacts = 0, nlocal = atom->nlocal;
    int *mask = atom->mask;

    for(int i = 0; i < nlocal; i++)
        if(mask[i] & contact_groupbit)
           ncontacts += npartner_[i];
    return ncontacts;
  }

#endif
