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

namespace LAMMPS_NS {

	class FixHeatGranRad : public FixHeatGran {

	public:
		FixHeatGranRad(class LAMMPS *, int, char **);
		~FixHeatGranRad();

		int setmask();
		void init();
		void post_force(int);

	private:
		double TB;	//NP background temperature
	};
}

#endif
#endif