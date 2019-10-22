//electrode.h
#ifndef ELECTRODE_H
#define ELECTRODE_H
#include <iostream>
#include <string>
#include "molecule.h"
#include "transform.h"

class Electrode
{
	Transform transformer;
	void JuncLocate(Molecule&,int&,int&);
	int NearestH(Molecule&,int,char);
public:
	void LinearElectrode(Molecule&,int,std::string,int *,double,double);
	void HydroRepElectrode(Molecule&,int,std::string,int *,double,double);
	void TipSquareElectrode(Molecule&,std::string,int *,double,double);
	void PyramidElectrode(Molecule&,int,std::string,int *,double,double);
};
#endif
