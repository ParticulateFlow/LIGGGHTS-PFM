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

/*
Contributing author: Stefan Amberger (JKU, DCS Computing)
*/

#include "neighbor.h"
#include <cmath>

using namespace LAMMPS_NS;

/*
	bin2XYZ writes the dimentional indices of bin i into the variables
	ix, iy and iz
	these are indices for local bins
*/
void Neighbor::bin2XYZ(int bin, int &ix, int &iy, int &iz){

	ix = mbinxlo + (bin % mbinx);
	iy = mbinylo + (((bin - (ix-mbinxlo))/(mbinx)) % mbiny);
	iz = mbinzlo + (bin - (ix-mbinxlo) - (iy-mbinylo)*mbinx)/(mbinx*mbiny);

	/*NL*/ //CHRISTOPH ix = (bin % (mbiny*mbinx)) % mbinx + mbinxlo;
    /*NL*/ //CHRISTOPH iy = static_cast<int>(round(static_cast<double>(((bin - mbinxlo) % (mbiny*mbinx)))/static_cast<double>(mbinx))) + mbinylo;
    /*NL*/ //CHRISTOPH iz = static_cast<int>(round(static_cast<double>((bin - mbinxlo - mbinylo*mbinx))/static_cast<double>((mbiny*mbinx))))+ mbinzlo;


	//NP int yzpart;

	//NP ix = i % mbinx;
	//NP yzpart = (i - ix) / mbinx;
	//NP iy = yzpart % mbiny;
	//NP iz = (yzpart - iy) / mbiny;
}
//NP void Neighbor::bin2XYZ(int i, int &ix, int &iy, int &iz){
//NP 	int yzpart;

//NP 	ix = i % mbinx;
//NP 	yzpart = (i - ix) / mbinx;
//NP 	iy = yzpart % mbiny;
//NP 	iz = (yzpart - iy) / mbiny;
//NP }

/*
	XYZ2bin returns the sequential index corresponding to the directional
	indices ix, iy and iz
	these are indices for local bins
*/
int Neighbor::XYZ2bin(int ix, int iy, int iz){
	return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
}

/* ---------------------------------------------------------------------- */

void Neighbor::binBorders(int ibin, double &xlo, double &xhi, double &ylo, double &yhi, double &zlo, double &zhi){
	/*NL*/ bool debugflag = false;
	int ix, iy, iz;

	bin2XYZ(ibin, ix, iy, iz);

	xlo = bboxlo[0] + (ix)*binsizex;
	xhi = bboxlo[0] + (ix + 1)*binsizex;
	ylo = bboxlo[1] + (iy)*binsizey;
	yhi = bboxlo[1] + (iy + 1)*binsizey;
	zlo = bboxlo[2] + (iz)*binsizez;
	zhi = bboxlo[2] + (iz + 1)*binsizez;

	/*NL*/ if(debugflag) printf("xrange: [%f %f]\n", xlo,xhi);
	/*NL*/ if(debugflag) printf("yrange: [%f %f]\n", ylo,yhi);
	/*NL*/ if(debugflag) printf("zrange: [%f %f]\n", zlo,zhi);
}


/*
	return -1 if resulting bin would be of out boundary
	return sequential index of requested offset bin of
	bin i otherwise
*/

int Neighbor::binHop(int i, int x, int y, int z){
	/*NL*/ bool debugflag = false;
	int ix, iy, iz;

	bin2XYZ(i, ix, iy, iz);
	/*NL*/ if(debugflag) printf("current bin: %d: %d %d %d\n", i,ix,iy,iz);
	/*NL*/ if(debugflag) printf("x,y,z: %d %d %d\n", x,y,z);
	/*NL*/ if(debugflag) printf("mbinx, mbinxlo: %d %d\n", mbinx,mbinxlo);
	/*NL*/ if(debugflag) printf("mbiny, mbinylo: %d %d\n", mbiny,mbinylo);
	/*NL*/ if(debugflag) printf("mbinz, mbinzlo: %d %d\n", mbinz,mbinzlo);
	/*NL*/ if(debugflag) printf("sx, sy, sz: %d %d %d\n", sx,sy,sz);

	/*NL*/ /*ORIGINAL*/// if (ix + x < 0 || ix + x >= mbinx || 0 > iy + y || iy + y >= mbiny || 0 > iz + z || iz + z >= mbinz)
	/*NL*/ /*SECOND*///if ((ix + x < mbinxlo) || ((ix + x - mbinxlo) >= mbinx) || (iy + y < mbinylo) || ((iy + y - mbinylo) >= mbiny) || (iz + z < mbinzlo) || ((iz + z - mbinzlo) >= mbinz))
	/*NL*/ /*THIRD*/// if ((ix + x - sx < mbinxlo) || ((ix + x + sx - mbinxlo) >= mbinx) || (iy + y - sy < mbinylo) || ((iy + y + sy - mbinylo) >= mbiny) || (iz + z - sz < mbinzlo) || ((iz + z + sz - mbinzlo) >= mbinz))
	if ((ix + x < mbinxlo) || ((ix + x - mbinxlo) >= mbinx) || (iy + y < mbinylo) || ((iy + y - mbinylo) >= mbiny) || (iz + z < mbinzlo) || ((iz + z - mbinzlo) >= mbinz))
		return -1;

	ix += x;
	iy += y;
	iz += z;
	/*NL*/ if(debugflag) printf("next bin: %d: %d %d %d\n", XYZ2bin(ix, iy, iz),ix,iy,iz);
	return XYZ2bin(ix, iy, iz);
}
