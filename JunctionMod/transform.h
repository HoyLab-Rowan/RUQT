//transform.h
#ifndef TRANSFORM_H
#define TRANSFORM_H
#include <vector>
#include <string>
#include "molecule.h"
class Transform
{
	
	public:
	void ConvertCoord(Molecule&, char, char);
	void CoordEdit(Molecule&, double, double, double);
	void SelectiveEdit(Molecule&, double, double, double, std::vector<int>&);
	void Center(Molecule&, int);
	void AlignAxis(Molecule&, char, int, int);
	void BondRotate(Molecule&, double, int, int, std::vector<int>&);
	/*
	Functions labeled ODSP are oddly specific types of transformations.
	These transformations are typically not made as a child for another function
	but are neccessary for specific observations in research.
	*/
	
	void DeviceSwitch_ODSP(Molecule&, int *, double, double, std::vector<int>&, std::vector<int>&, std::string);
};
#endif
