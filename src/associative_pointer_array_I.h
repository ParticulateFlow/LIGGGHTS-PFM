/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-2015 DCS Computing GmbH, Linz
   Copyright 2015-     JKU Linz

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
   Richard Berger (JKU Linz)
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_I_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_I_H

  /* ----------------------------------------------------------------------
   constructors, destructor
  ------------------------------------------------------------------------- */

  template<typename T>
  AssociativePointerArray<T>::AssociativePointerArray()
  {
  }

  template<typename T>
  AssociativePointerArray<T>::~AssociativePointerArray()
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      delete it->second;
  }

  /* ----------------------------------------------------------------------
   add for per-element and per-mesh properties
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::add(const char *_id, const char* _comm, const char* _ref, const char *_restart, int _scalePower)
  {
    content_[_id] = static_cast<T*>(new U(_id,_comm,_ref,_restart,_scalePower));

    return static_cast<U*>(content_[_id]);
  }

  /* ----------------------------------------------------------------------
   delete properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::remove(const char *_id)
  {
    content_iterator it = content_.find(_id);
    if(it != content_.end())
    {
      delete it->second;
      content_.erase(it);
    }
  }

  /* ----------------------------------------------------------------------
   check if all have the same length
  ------------------------------------------------------------------------- */

  template<typename T>
  bool AssociativePointerArray<T>::sameLength(int _len)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      if(it->second->size() != _len)
        return false;
    return true;
  }

  /* ----------------------------------------------------------------------
   get pointer to property
  ------------------------------------------------------------------------- */

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerById(const char *_id)
  {
    content_iterator it = content_.find(_id);
    if(it != content_.end())
      return dynamic_cast<U*>(it->second);
    else
      return 0;
  }

  template<typename T>
  T* AssociativePointerArray<T>::getBasePointerById(const char *_id)
  {
    content_iterator it = content_.find(_id);
    if(it != content_.end())
      return it->second;
    else
      return 0;
  }

  template<typename T> template<typename U>
  U* AssociativePointerArray<T>::getPointerByIndex(int i)
  {
    if(i >= size() || i < 0) return 0;
    else return dynamic_cast<U*>(content_[i]);
  }

  template<typename T>
  T* AssociativePointerArray<T>::getBasePointerByIndex(int i) const
  {
    if(i >= size() || i < 0) return 0;
    else return content_[i];
  }

  /* ----------------------------------------------------------------------
   memory management
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::grow(int to)
  {
    int by;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
    {
      by = to - it->second->size();
      if(by > 0)
        it->second->addUninitialized(by);
    }
  }

  template<typename T>
  int AssociativePointerArray<T>::size() const
  {
    return content_.size();
  }

  /* ----------------------------------------------------------------------
   copy data from element from to element to
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::copyElement(int from, int to)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->copy(from,to);
  }

  /* ----------------------------------------------------------------------
   add an element and initialize its properties with 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::addUninitializedElement()
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->addUninitialized(1);
  }

  /* ----------------------------------------------------------------------
   add an element and initialize its properties with 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::addZeroElement()
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->addZero();
  }

  /* ----------------------------------------------------------------------
   delete element n
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteElement(int n)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->del(n);
  }

  /* ----------------------------------------------------------------------
   delete forward properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteForwardElement(int n,bool scale,bool translate,bool rotate)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->delForward(n,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   delete restart properties of element i
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteRestartElement(int n,bool scale,bool translate,bool rotate)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->delRestart(n,scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   delete restart properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::deleteRestartGlobal(bool scale,bool translate,bool rotate)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->delRestart(scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   clear reverse properties, i.e. reset all of them to 0
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::clearReverse(bool scale,bool translate,bool rotate)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->clearReverse(scale,translate,rotate);
  }

  /* ----------------------------------------------------------------------
   has id
  ------------------------------------------------------------------------- */

  template<typename T>
  bool AssociativePointerArray<T>::hasId(const char *_id)
  {
      content_iterator it = content_.find(_id);
      return (it != content_.end());

  }

  /* ----------------------------------------------------------------------
   store original value for reset
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::storeOrig(AssociativePointerArray &orig)
  {
    for(content_iterator it = orig.content_.begin(); it != orig.content_.end(); ++it)
        it->second->setFromContainer(content_[it->first]);
  }

  template<typename T>
  void AssociativePointerArray<T>::storeOrig(const char *_id, AssociativePointerArray &orig)
  {
    for(content_iterator it = orig.content_.begin(); it != orig.content_.end(); ++it)
      if(content_[it->first]->matches_id(_id))
        it->second->setFromContainer(content_[it->first]);
  }

  /* ----------------------------------------------------------------------
   reset to original value
  ------------------------------------------------------------------------- */

  template<typename T>
  bool AssociativePointerArray<T>::reset(AssociativePointerArray &orig)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->setFromContainer(orig.content_[it->first]);

    return true;
  }

  template<typename T>
  bool AssociativePointerArray<T>::reset(const char *_id, AssociativePointerArray &orig)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      if(it->second->matches_id(_id))
        it->second->setFromContainer(orig.content_[it->first]);

    return true;
  }

  template<typename T>
  void AssociativePointerArray<T>::setToDefault(int n)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      if(it->second->useDefault())
        it->second->setToDefault(n);
  }

  /* ----------------------------------------------------------------------
   move, rotate scale all properties
  ------------------------------------------------------------------------- */

  template<typename T>
  void AssociativePointerArray<T>::rotate(const double *dQ)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->rotate(dQ);
  }

  template<typename T>
  void AssociativePointerArray<T>::scale(double factor)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->scale(factor);
  }

  template<typename T>
  void AssociativePointerArray<T>::move(const double *delta)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->move(delta);
  }

  template<typename T>
  void AssociativePointerArray<T>::moveElement(int n, const double *delta)
  {
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      it->second->moveElement(n,delta);
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for all elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::bufSize(int operation,bool scale,bool translate,bool rotate) const
  {
    int buf_size = 0;
    for(content_const_iterator it = content_.begin(); it != content_.end(); ++it)
      buf_size += it->second->bufSize(operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nsend += it->second->pushToBuffer(&(buf[nsend]),operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nrecv += it->second->popFromBuffer(&(buf[nrecv]),operation,scale,translate,rotate);
    return nrecv;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for list of elements
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      buf_size += it->second->elemListBufSize(n,operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nsend += it->second->pushElemListToBuffer(n,list,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nrecv += it->second->popElemListFromBuffer(first,n,&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nrecv += it->second->pushElemListToBufferReverse(first,n,&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nsend += it->second->popElemListFromBufferReverse(n,list,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  /* ----------------------------------------------------------------------
   buf size, push, pop for single element
  ------------------------------------------------------------------------- */

  template<typename T>
  int AssociativePointerArray<T>::elemBufSize(int operation,bool scale,bool translate,bool rotate)
  {
    int buf_size = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      buf_size += it->second->elemBufSize(operation,scale,translate,rotate);
    return buf_size;
  }

  template<typename T>
  int AssociativePointerArray<T>::pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nsend = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nsend += it->second->pushElemToBuffer(n,&buf[nsend],operation,scale,translate,rotate);
    return nsend;
  }

  template<typename T>
  int AssociativePointerArray<T>::popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate)
  {
    int nrecv = 0;
    for(content_iterator it = content_.begin(); it != content_.end(); ++it)
      nrecv += it->second->popElemFromBuffer(&buf[nrecv],operation,scale,translate,rotate);
    return nrecv;
  }

#endif
