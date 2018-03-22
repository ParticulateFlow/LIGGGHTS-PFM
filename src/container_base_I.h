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

/* ----------------------------------------------------------------------
   Contributing authors:
   Christoph Kloss (JKU Linz, DCS Computing GmbH, Linz)
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_CONTAINER_BASE_I_H
#define LMP_CONTAINER_BASE_I_H

  /* ----------------------------------------------------------------------
   definition of reference frames and comm types
  ------------------------------------------------------------------------- */

  // reference frame types
  // invariant: invariant to scaling, translation, rotation
  // trans invariant: invariant to translation, not invariant to scaling, rotation
  // trans+rot invariant: invariant to translation, rotation, not invariant to scaling
  // general: not invariant to scaling, translation, rotation

  enum{ REF_FRAME_UNDEFINED,
        REF_FRAME_INVARIANT,
        REF_FRAME_SCALE_TRANS_INVARIANT,
        REF_FRAME_TRANS_ROT_INVARIANT,
        REF_FRAME_TRANS_INVARIANT,
        REF_FRAME_GENERAL};

  // communication types

  enum{ // communication invoked manually
        COMM_TYPE_MANUAL,
        // only exchange and borders comm
        COMM_EXCHANGE_BORDERS,
        // forward comm every step
        COMM_TYPE_FORWARD,
        // forward comm based on reference frame setting
        // ie if mesh rotates, egdeVecs are communicated
        //NP does exchange, borders with buffer-initialized values
        COMM_TYPE_FORWARD_FROM_FRAME,
        // reverse comm every step
        //NP does exchange and borders with 0-initialized values
        COMM_TYPE_REVERSE,
        // no comm at all
        //NP does exchange and borders with 0-initialized values
        COMM_TYPE_NONE,
        // undefined state for error check
        COMM_TYPE_UNDEFINED};  // communication types

  // restart types

  enum{ RESTART_TYPE_UNDEFINED,
        RESTART_TYPE_YES,
        RESTART_TYPE_NO};


  /* ----------------------------------------------------------------------
   decide if property is pushed or pulled at all
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decidePackUnpackOperation(int operation,bool scale,bool translate, bool rotate) const
  {
      // return true for manual communication, such as for node_, node_orig_
      // etc in MultiNodeMeshParallel
      if(COMM_TYPE_MANUAL == communicationType_)
        return true;

      //NP check for restart
      if(OPERATION_RESTART == operation)
      {
          if(restartType_ == RESTART_TYPE_YES)
            return true;
          return false;
      }

      //NP communication in exchange() and borders() steps is always performed
      if(OPERATION_COMM_BORDERS == operation ||
         OPERATION_COMM_EXCHANGE == operation )
        return true;

      //NP no fw or reverse comm for COMM_TYPE_NONE, but do exchange or borders
      if(COMM_TYPE_NONE == communicationType_)
        return false;

      if(OPERATION_COMM_REVERSE == operation &&
         COMM_TYPE_REVERSE == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == operation &&
         COMM_TYPE_FORWARD == communicationType_)
        return true;

      if(OPERATION_COMM_FORWARD == operation &&
         COMM_TYPE_FORWARD_FROM_FRAME == communicationType_)
      {
         if(scale && !isScaleInvariant())
           return true;
         if(translate && !isTranslationInvariant())
           return true;
         if(rotate && !isRotationInvariant())
           return true;

         return false;
      }

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   decide if operation performs data communication
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decideCommOperation(int operation) const
  {
      //NP have to decide at unpack if data is initialized with 0
      //NP or pulled from buffer

      //NP e.g. at exchange:
      //NP      forces would be initialized with 0
      //NP      positions would be initialized from buffer data

      //NP restart always pulls from buffer
      if(operation == OPERATION_RESTART)
          return true;

      //NP forward and reverse comm always pull from buffer
      //NP (thats why they are done)
      if(operation == OPERATION_COMM_FORWARD ||
         operation == OPERATION_COMM_REVERSE )
        return true;


      //NP exchange() and borders()
      if(operation == OPERATION_COMM_BORDERS ||
         operation == OPERATION_COMM_EXCHANGE )
      {
          //NP comm none and comm reverse dont pull from buffer
          if(communicationType_ == COMM_TYPE_NONE ||
             communicationType_ == COMM_TYPE_REVERSE)
             return false;

          //NP all others do
          return true;
      }

      // default
      return true;
  }

  /* ----------------------------------------------------------------------
   decide if unpack creates new element or overwrites existing data
  ------------------------------------------------------------------------- */

  inline bool ContainerBase::decideCreateNewElements(int operation) const
  {
      //NP have to decide at unpack if new elements are created or
      //NP existing ones are over-written

      //NP restart always creates new elements
      if(operation == OPERATION_RESTART)
          return true;

      //NP exchange() and borders() always create new elements
      if(operation == OPERATION_COMM_BORDERS ||
         operation == OPERATION_COMM_EXCHANGE )
        return true;

      //NP forward and reverse comm never create new elements
      if(operation == OPERATION_COMM_FORWARD ||
         operation == OPERATION_COMM_REVERSE )
        return false;

      // default
      return false;
  }

  /* ----------------------------------------------------------------------
   fast test for reference frame
   note that rotation is only carried out for LEN_VEC==3
  ------------------------------------------------------------------------- */

    bool ContainerBase::isScaleInvariant() const
    {
       return ( refFrame_ == REF_FRAME_INVARIANT ||
                refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT);
    }

    bool ContainerBase::isTranslationInvariant() const
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT ||
                 refFrame_ == REF_FRAME_SCALE_TRANS_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_INVARIANT);
    }

    bool ContainerBase::isRotationInvariant() const
    {
        return ( refFrame_ == REF_FRAME_INVARIANT ||
                 refFrame_ == REF_FRAME_TRANS_ROT_INVARIANT ||
                 lenVec() != 3);
    }

  /* ----------------------------------------------------------------------
   ID operations
  ------------------------------------------------------------------------- */

  inline void ContainerBase::id(char *_id)
  {
      strcpy(_id,id_);
  }

  inline bool ContainerBase::matches_id(const char *_id)
  {
      if(strcmp(_id,id_) == 0) return true;
      return false;
  }

#endif
