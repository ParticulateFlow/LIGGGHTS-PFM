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
   Daniel Queteschiner (JKU Linz)
------------------------------------------------------------------------- */

#ifndef LMP_ASSOCIATIVE_POINTER_ARRAY_H
#define LMP_ASSOCIATIVE_POINTER_ARRAY_H

#include <string>
#include <map>
#include "memory.h"

namespace LAMMPS_NS
{
  #define ID_LEN 100

template<typename T>
class AssociativePointerArray
{
      public:
        AssociativePointerArray();
        ~AssociativePointerArray();

        template <typename U>
        U* add(const char *_id, const char* _comm, const char* _ref, const char *_restart,int _scalePower = 1);

        void remove(const char *_id);

        template <typename U>
        U* getPointerById(const char *_id);

        T* getBasePointerById(const char *_id);

        template <typename U>
        U* getPointerByIndex(int i);

        T* getBasePointerByIndex(int i) const;

        void grow(int to);

        int size() const;

        bool sameLength(int _len);

        inline void copyElement(int from, int to);
        inline void addUninitializedElement();
        inline void addZeroElement();
        inline void deleteElement(int n);
        inline void deleteForwardElement(int n,bool scale,bool translate,bool rotate);
        inline void deleteRestartElement(int n,bool scale,bool translate,bool rotate);
        inline void deleteRestartGlobal(bool scale,bool translate,bool rotate);

        inline void clearReverse(bool scale,bool translate,bool rotate);

        inline void storeOrig(class AssociativePointerArray &orig);
        inline void storeOrig(const char *_id,class AssociativePointerArray &orig);
        inline bool reset(class AssociativePointerArray &orig);
        inline bool reset(const char *_id,class AssociativePointerArray &orig);
        inline void setToDefault(int n);

        void rotate(const double *dQ);
        void move(const double *delta);
        void moveElement(int i,const double *delta);
        void scale(double factor);

        inline int bufSize(int operation,bool scale,bool translate,bool rotate) const;
        inline int pushToBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemListBufSize(int n,int operation,bool scale,bool translate,bool rotate);
        inline int pushElemListToBuffer(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBuffer(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int pushElemListToBufferReverse(int first, int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemListFromBufferReverse(int n, int *list, double *buf, int operation,bool scale,bool translate, bool rotate);

        inline int elemBufSize(int operation,bool scale,bool translate,bool rotate);
        inline int pushElemToBuffer(int n, double *buf, int operation,bool scale,bool translate, bool rotate);
        inline int popElemFromBuffer(double *buf, int operation,bool scale,bool translate, bool rotate);

        bool hasId(const char *_id);

      private:

        std::map<std::string, T*> content_;
        typedef typename std::map<std::string, T*>:: iterator content_iterator;
        typedef typename std::map<std::string, T*>:: const_iterator content_const_iterator;
};

  // *************************************
  #include "associative_pointer_array_I.h"
  // *************************************

} /* LAMMPS_NS */
#endif /* ASSOCIATIVEPOINTERARRAY_H_ */
