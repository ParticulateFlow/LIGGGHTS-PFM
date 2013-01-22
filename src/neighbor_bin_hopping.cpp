#include "neighbor.h"
#include <cmath>

using namespace LAMMPS_NS;

/*
	bin2XYZ writes the dimentional indices of bin i into the variables
	ix, iy and iz
	these are indices for local bins
*/
void Neighbor::bin2XYZ(int bin, int &ix, int &iy, int &iz){

	ix = (bin % (mbiny*mbinx)) % mbinx + mbinxlo;
    iy = static_cast<int>(round(static_cast<double>(((bin - mbinxlo) % (mbiny*mbinx)))/static_cast<double>(mbinx))) + mbinylo;
    iz = static_cast<int>(round(static_cast<double>((bin - mbinxlo - mbinylo*mbinx))/static_cast<double>((mbiny*mbinx))))+ mbinzlo;


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
	int ix, iy, iz;

	bin2XYZ(ibin, ix, iy, iz);

	xlo = bboxlo[0] + (mbinxlo + ix)*binsizex;
	xhi = bboxlo[0] + (mbinxlo + ix + 1)*binsizex;
	ylo = bboxlo[1] + (mbinylo + iy)*binsizey;
	yhi = bboxlo[1] + (mbinylo + iy + 1)*binsizey;
	zlo = bboxlo[2] + (mbinzlo + iz)*binsizez;
	zhi = bboxlo[2] + (mbinzlo + iz + 1)*binsizez;
}


/*
	return -1 if resulting bin would be of out boundary
	return sequential index of requested offset bin of
	bin i otherwise
*/

int Neighbor::binHop(int i, int x, int y, int z){
	int ix, iy, iz;

	bin2XYZ(i, ix, iy, iz);

	/*ORIGINAL*/// if (ix + x < 0 || ix + x >= mbinx || 0 > iy + y || iy + y >= mbiny || 0 > iz + z || iz + z >= mbinz)
	if (ix + x < mbinxlo || (ix + x - mbinxlo) >= mbinx || iy + y < mbinylo || (iy + y - mbinylo) >= mbiny || iz + z < mbinzlo || (iz + z - mbinzlo) >= mbinz)
		return -1;

	ix += x;
	iy += y;
	iz += z;

	return XYZ2bin(ix, iy, iz);
}
