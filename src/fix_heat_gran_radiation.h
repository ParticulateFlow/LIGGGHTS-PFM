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

#ifdef FIX_CLASS

FixStyle(heat/gran/radiation,FixHeatGranRad)

#else

#ifndef LMP_FIX_HEATGRAN_RADIATION_H
#define LMP_FIX_HEATGRAN_RADIATION_H

#include "fix_heat_gran.h"
#include "random_mars.h"

namespace LAMMPS_NS {

	class FixHeatGranRad : public FixHeatGran {

	public:
		FixHeatGranRad(class LAMMPS *, int, char **);
		~FixHeatGranRad();
    void pre_delete(bool){};

		int setmask();
		void init();
		void post_force(int);

	private:
    bool intersectRaySphere(double*, double*, double*, double, double&, double*);
    void randDir(double*, double*);
    void randOnSphere(double*, double, double*, double*);

    double TB;   // background temperature
    double Qtot; // total energy of all particles (heat)
    double Qr;   // energy of one ray

    class FixPropertyGlobal* fix_emissivity;
    double *emissivity;

    double sigma; // stefan bolzmann constant
    class RanMars rGen;
	};
}

#endif
#endif