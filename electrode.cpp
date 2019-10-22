//electrode.cpp
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "molecule.h"
#include "transform.h"
#include "electrode.h"
void Electrode::JuncLocate(Molecule& examMole, int& atom1, int& atom2)
{
	short int ithJunc = 0;
	atom1 = 0;
	atom2 = 0;
	for (int i=0;i<examMole.GetAtomCount();i++)
	{
		if (examMole.GetAtomSym(i) == "S" && ithJunc == 0){atom1 = i; ithJunc++;}
		else if (examMole.GetAtomSym(i) == "S" && ithJunc == 1){atom2 = i;}
	}
	return;
}

int Electrode::NearestH(Molecule& examMole, int ithAtom, char Decision)
{
	double atom1x = examMole.GetAtomCoord(ithAtom, 'x');
	double atom1y = examMole.GetAtomCoord(ithAtom, 'y');
	double atom1z = examMole.GetAtomCoord(ithAtom, 'z');
	double atom2x, atom2y, atom2z;
	double NearestDist = 0.0;
	double Dist;
	int NearHnum; 
	for (int i=0;i<examMole.GetAtomCount();i++)
	{
		if (examMole.GetAtomSym(i) == "H")
		{
			atom2x = examMole.GetAtomCoord(i,'x');
			atom2y = examMole.GetAtomCoord(i,'y');
			atom2z = examMole.GetAtomCoord(i,'z');
			Dist = sqrt(pow(atom2x-atom1x,2.0)+pow(atom2y-atom1y,2.0)+pow(atom2z-atom1z,2.0));
			if (NearestDist == 0.0) { NearestDist = Dist; NearHnum = i;}
			if (NearestDist > Dist) { NearestDist = Dist; NearHnum = i;}
		}
	}
	if (Decision == 'r'){examMole.RemoveAtom(NearHnum); return 0;}
	if (Decision == 'f'){return NearHnum;}
	
}

void Electrode::LinearElectrode(Molecule& examMole, int MCount, std::string MType, int * juncAtoms, double BL1, double BL2)
{
	int junc1 = *juncAtoms; int junc2 = *(juncAtoms+2);
	int Hydro1 = *(juncAtoms+1); int Hydro2 = *(juncAtoms+3);
	const double junc1coord = examMole.GetAtomCoord(junc1,'z');
	const double junc2coord = examMole.GetAtomCoord(junc2,'z');
	transformer.AlignAxis(examMole,'z',junc1,junc2);
	examMole.AddAtom(MType,0.0,0.0,junc1coord-BL1);	
	examMole.AddAtom(MType,0.0,0.0,junc2coord+BL1);
	for (int i = 2; i <= MCount; i++)
	{
		examMole.AddAtom(MType,0.0,0.0,junc1coord-BL1-(BL2*(i-1)));
		examMole.AddAtom(MType,0.0,0.0,junc2coord+BL1+(BL2*(i-1)));
	}
	if (Hydro2 > Hydro1) {examMole.RemoveAtom(Hydro2); examMole.RemoveAtom(Hydro1);}
	else {examMole.RemoveAtom(Hydro1); examMole.RemoveAtom(Hydro2);}

	return;
}

void Electrode::HydroRepElectrode(Molecule& examMole, int MCount, std::string MType, int * juncAtoms, double BL1, double BL2)
{
	int junc1 = *juncAtoms; int junc2 = *(juncAtoms+2);
	int Hydro1 = *(juncAtoms+1); int Hydro2 = *(juncAtoms+3);

	transformer.AlignAxis(examMole,'z',junc1,junc2);
	double S1x = examMole.GetAtomCoord(junc1,'x'); double S1y = examMole.GetAtomCoord(junc1,'y'); double S1z = examMole.GetAtomCoord(junc1,'z');
	double H1x = examMole.GetAtomCoord(Hydro1,'x'); double H1y = examMole.GetAtomCoord(Hydro1,'y'); double H1z = examMole.GetAtomCoord(Hydro1,'z');
	double SCoeff1 = BL1/sqrt(pow(S1x-H1x,2.0)+pow(S1y-H1y,2.0)+pow(S1z-H1z,2.0));
	double Au1x = SCoeff1*(H1x-S1x)+S1x; double Au1y = SCoeff1*(H1y-S1y)+S1y; double Au1z = SCoeff1*(H1z-S1z)+S1z;
	double S2x = examMole.GetAtomCoord(junc2,'x'); double S2y = examMole.GetAtomCoord(junc2,'y'); double S2z = examMole.GetAtomCoord(junc2,'z');
	double H2x = examMole.GetAtomCoord(Hydro2,'x'); double H2y = examMole.GetAtomCoord(Hydro2,'y'); double H2z = examMole.GetAtomCoord(Hydro2,'z');
	double SCoeff2 = BL1/sqrt(pow(S2x-H2x,2.0)+pow(S2y-H2y,2.0)+pow(S2z-H2z,2.0));
	double Au2x = SCoeff2*(H2x-S2x)+S2x; double Au2y = SCoeff2*(H2y-S2y)+S2y; double Au2z = SCoeff2*(H2z-S2z)+S2z;

	examMole.AddAtom(MType,Au1x,Au1y,Au1z);
	examMole.AddAtom(MType,Au2x,Au2y,Au2z);
	for (int i = 2; i<= MCount; i++){ examMole.AddAtom(MType,Au1x,Au1y,Au1z-(BL2*(i-1))); examMole.AddAtom(MType,Au2x,Au2y,Au2z+(BL2*(i-1))); }
	if (Hydro2 > Hydro1) {examMole.RemoveAtom(Hydro2); examMole.RemoveAtom(Hydro1);}
	else {examMole.RemoveAtom(Hydro1); examMole.RemoveAtom(Hydro2);}
	return;
}

void Electrode::TipSquareElectrode(Molecule& examMole, std::string MType, int * juncAtoms, double BL1, double BL2) 
{
	int junc1 = *juncAtoms; int junc2 = *(juncAtoms+1);
	int Hydro1 = *(juncAtoms+2); int Hydro2 = *(juncAtoms+3);
	transformer.AlignAxis(examMole,'z',junc1,junc2);
	double S2D = examMole.GetAtomCoord(junc2,'z');

	examMole.AddAtom(MType,0.0,0.0,-BL1);
	examMole.AddAtom(MType,BL2/2.0,BL2/2.0,-BL1-(BL2/pow(2.0,0.5)));
	examMole.AddAtom(MType,-BL2/2.0,BL2/2.0,-BL1-(BL2/pow(2.0,0.5)));
	examMole.AddAtom(MType,BL2/2.0,-BL2/2.0,-BL1-(BL2/pow(2.0,0.5)));
	examMole.AddAtom(MType,-BL2/2.0,-BL2/2.0,-BL1-(BL2/pow(2.0,0.5)));

	examMole.AddAtom(MType,0.0,0.0,BL1+S2D);
	examMole.AddAtom(MType,BL2/2.0,BL2/2.0,BL1+(BL2/pow(2.0,0.5)+S2D));
	examMole.AddAtom(MType,-BL2/2.0,BL2/2.0,BL1+(BL2/pow(2.0,0.5)+S2D));
	examMole.AddAtom(MType,BL2/2.0,-BL2/2.0,BL1+(BL2/pow(2.0,0.5)+S2D));
	examMole.AddAtom(MType,-BL2/2.0,-BL2/2.0,BL1+(BL2/pow(2.0,0.5)+S2D));
	if (Hydro2 > Hydro1) {examMole.RemoveAtom(Hydro2); examMole.RemoveAtom(Hydro1);}
	else {examMole.RemoveAtom(Hydro1); examMole.RemoveAtom(Hydro2);}
	return;
}

void Electrode::PyramidElectrode(Molecule& examMole, int layers, std::string MType, int * juncAtoms, double BL1, double BL2)
{
	const double V1 = pow(3.0,0.5)/3.0; const double V2 = pow(2.0/3.0,0.5); const double V3 = pow(3.0,0.5)/2.0;
	int junc1 = *juncAtoms; int junc2 = *(juncAtoms+1);
	int Hydro1 = *(juncAtoms+2); int Hydro2 = *(juncAtoms+3);
	transformer.AlignAxis(examMole,'z',junc1,junc2);
	double S2D = examMole.GetAtomCoord(junc2,'z');

	for (int i=1; i<=layers; i++)
	{
		examMole.AddAtom(MType,0.0,(i-1)*BL2*V1,-(i-1)*BL2*V2-BL1);
		for (int j=1; j<i; j++)
		{
			examMole.AddAtom(MType,j*0.5*BL2,(i-1)*BL2*V1-j*BL2*V3,-(i-1)*BL2*V2-BL1);
			for (int k=1; k<=j; k++)
			{
				examMole.AddAtom(MType,j*0.5*BL2-BL2*k,(i-1)*BL2*V1-j*BL2*V3,-(i-1)*BL2*V2-BL1);
			}
		}
	}
	for (int i=1; i<=layers; i++)
	{
		examMole.AddAtom(MType,0.0,(i-1)*BL2*V1,(i-1)*BL2*V2+BL1+S2D);
		for (int j=1; j<i; j++)
		{
			examMole.AddAtom(MType,j*0.5*BL2,(i-1)*BL2*V1-j*BL2*V3,(i-1)*BL2*V2+BL1+S2D);
			for (int k=1; k<=j; k++)
			{
				examMole.AddAtom(MType,j*0.5*BL2-BL2*k,(i-1)*BL2*V1-j*BL2*V3,(i-1)*BL2*V2+BL1+S2D);
			}
		}
	}
	if (Hydro2 > Hydro1) {examMole.RemoveAtom(Hydro2); examMole.RemoveAtom(Hydro1);}
	else {examMole.RemoveAtom(Hydro1); examMole.RemoveAtom(Hydro2);}
	return;
}
