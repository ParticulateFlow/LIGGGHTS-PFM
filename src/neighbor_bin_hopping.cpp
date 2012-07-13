#include "neighbor.h"

using namespace LAMMPS_NS;

/*
	bin2XYZ writes the dimentional indices of bin i into the variables
	ix, iy and iz
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
*/
int Neighbor::XYZ2bin(int ix, int iy, int iz){
	return iz*mbiny*mbinx + iy*mbinx + ix;
}

/*
	functions binUpxxx(int i)

	return -1 if resulting bin would be of out boundary
	return sequential index of neighboring (up, xxx-direction) bin of
	bin i otherwise

	functions binDownxxx(int i)

	return -1 if resulting bin would be of out boundary
	return sequential index of neighboring (down, xxx-direction) bin of
	bin i otherwise
*/

int Neighbor::binUpX(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (ix == mbinx - 1)
		return -1;

	return XYZ2bin(ix + 1, iy, iz);
}

int Neighbor::binDownX(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (ix == 0)
		return -1;

	return XYZ2bin(ix - 1, iy, iz);
}

int Neighbor::binUpY(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (iy == mbiny - 1)
		return -1;

	return XYZ2bin(ix, iy + 1, iz);
}

int Neighbor::binDownY(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (iy == 0)
		return -1;

	return XYZ2bin(ix, iy - 1, iz);
}

int Neighbor::binUpZ(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (iz == mbinz - 1)
		return -1;

	return XYZ2bin(ix, iy, iz + 1);
}

int Neighbor::binDownZ(int i){
	int ix,iy,iz;
	
	bin2XYZ(i, ix, iy, iz);

	if (iz == 0)
		return -1;

	return XYZ2bin(ix, iy, iz - 1);
}

