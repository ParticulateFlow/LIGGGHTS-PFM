#include "neighbor.h"

using namespace LAMMPS_NS;

/*
	bin2XYZ writes the dimentional indices of bin i into the variables
	ix, iy and iz
	these are indices for local bins
*/
void Neighbor::bin2XYZ(int i, int &ix, int &iy, int &iz){
	int yzpart;

	ix = i % mbinx;
	yzpart = (i - ix) / mbinx;
	iy = yzpart % mbiny;
	iz = (yzpart - iy) / mbiny;
}

/*
	XYZ2bin returns the sequential index corresponding to the directional
	indices ix, iy and iz
	these are indices for local bins
*/
int Neighbor::XYZ2bin(int ix, int iy, int iz){
	return iz*mbiny*mbinx + iy*mbinx + ix;
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

	if (0 > ix + x || ix + x >= mbinx || 0 > iy + y || iy + y >= mbiny || 0 > iz + z || iz + z >= mbinz)
		return -1;

	ix += x;
	iy += y;
	iz += z;

	return XYZ2bin(ix, iy, iz);
}
