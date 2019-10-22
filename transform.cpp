//transform.cpp
#include <string>
#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>
#include "molecule.h"
#include "transform.h"
void Transform::ConvertCoord(Molecule& examMole, char fromType, char toType)
{
	if (fromType == 'C' && toType == '1')//Cartesian to Cylindrical_X
	{
		double r, theta; //x will remain the same. r and theta depend on y and z
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			r = sqrt(pow(examMole.GetAtomCoord(i,'y'),2)
				+pow(examMole.GetAtomCoord(i,'z'),2));
			theta = atan2(examMole.GetAtomCoord(i,'y'),
					examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'y',r);
			examMole.SetAtomCoord(i,'z',theta);
		}
	}
	
	else if (fromType == 'C' && toType == '2')//Cartesian to Cylindrical_Y
	{
		double r, theta;
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			r = sqrt(pow(examMole.GetAtomCoord(i,'x'),2)
				+pow(examMole.GetAtomCoord(i,'z'),2));
			theta = atan2(examMole.GetAtomCoord(i,'x'),
					examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'x',examMole.GetAtomCoord(i,'y'));
			examMole.SetAtomCoord(i,'y',r);
			examMole.SetAtomCoord(i,'z',theta);
		}	
	}

	else if (fromType == 'C' && toType == '3')//Cartesian to Cylindrical_Z
	{
		double r, theta;
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			r = sqrt(pow(examMole.GetAtomCoord(i,'x'),2)
				+pow(examMole.GetAtomCoord(i,'y'),2));
			theta = atan2(examMole.GetAtomCoord(i,'y'),
					examMole.GetAtomCoord(i,'x'));
			examMole.SetAtomCoord(i,'x',examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'y',r);
			examMole.SetAtomCoord(i,'z',theta);
		}	
	}
	
	else if (fromType == '1' && toType == 'C')//Cylindrical_X to Cartesian
	{
		double y, z;
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			y = examMole.GetAtomCoord(i,'y')
				*sin(examMole.GetAtomCoord(i,'z'));
			z = examMole.GetAtomCoord(i,'y')
				*cos(examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'y',y);
			examMole.SetAtomCoord(i,'z',z);
		}
	}
	
	else if (fromType == '2' && toType == 'C')//Cylindrical_Y to Cartesian
	{
		double x, z;
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			x = examMole.GetAtomCoord(i,'y')
				*sin(examMole.GetAtomCoord(i,'z'));
			z = examMole.GetAtomCoord(i,'y')
				*cos(examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'y',examMole.GetAtomCoord(i,'x'));
			examMole.SetAtomCoord(i,'x',x);
			examMole.SetAtomCoord(i,'z',z);
		}
	}
	
	else if (fromType == '3' && toType == 'C')//Cylindrical_Z to Cartesian
	{
		double x, y;
		for (int i=0; i<examMole.GetAtomCount(); i++)
		{
			y = examMole.GetAtomCoord(i,'y')
				*sin(examMole.GetAtomCoord(i,'z'));
			x = examMole.GetAtomCoord(i,'y')
				*cos(examMole.GetAtomCoord(i,'z'));
			examMole.SetAtomCoord(i,'z',examMole.GetAtomCoord(i,'x'));
			examMole.SetAtomCoord(i,'x',x);
			examMole.SetAtomCoord(i,'y',y);
		}
	}
	return;
}



void Transform::CoordEdit(Molecule& examMole, double xTrans, double yTrans, double zTrans)
{
	double tempXcoord, tempYcoord, tempZcoord;
	for(int i=0; i<examMole.GetAtomCount(); i++)
	{
		tempXcoord=examMole.GetAtomCoord(i,'x')+xTrans;
		tempYcoord=examMole.GetAtomCoord(i,'y')+yTrans;
		tempZcoord=examMole.GetAtomCoord(i,'z')+zTrans;
		examMole.SetAtomCoord(i,'x',tempXcoord);
		examMole.SetAtomCoord(i,'y',tempYcoord);
		examMole.SetAtomCoord(i,'z',tempZcoord);
	}
	return;
}



void Transform::SelectiveEdit(Molecule& examMole, double xTrans, double yTrans, double zTrans, std::vector<int>& atomSel)
{//0 0 theta atomSel
	double tempXcoord, tempYcoord, tempZcoord; int j;
	for(int i=0; i<atomSel.size(); i++)
	{
		j=atomSel[i];
		//tempXcoord=examMole.GetAtomCoord(j,'x')+xTrans;
		//tempYcoord=examMole.GetAtomCoord(j,'y')+yTrans;
		tempZcoord=examMole.GetAtomCoord(j,'z')+zTrans;
		//examMole.SetAtomCoord(j,'x',tempXcoord);
		//examMole.SetAtomCoord(j,'y',tempYcoord);
		examMole.SetAtomCoord(j,'z',tempZcoord);
	}
}



void Transform::Center(Molecule& examMole, int A1)
{
	double xTrans, yTrans, zTrans;
	xTrans = -(examMole.GetAtomCoord(A1,'x'));
	yTrans = -(examMole.GetAtomCoord(A1,'y'));
	zTrans = -(examMole.GetAtomCoord(A1,'z'));
	CoordEdit(examMole,xTrans,yTrans,zTrans);
	return;
}



void Transform::AlignAxis(Molecule& examMole, char axis, int A1, int A2)
{
	const double PI = 3.14159265358979323846;
	if (axis == 'z' || axis == 'Z')
	{
		Center(examMole, A1);
		ConvertCoord(examMole,'C','1');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z')));
		ConvertCoord(examMole,'1','C');
		ConvertCoord(examMole,'C','2');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z')));
		ConvertCoord(examMole,'2','C');
	}
	if (axis == 'y' || axis == 'Y')
	{
		Center(examMole,A1);
		ConvertCoord(examMole,'C','1');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z'))+(PI/2));
		ConvertCoord(examMole,'1','C');
		ConvertCoord(examMole,'C','3');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z'))+(PI/2));
		ConvertCoord(examMole,'3','C');
	}
	if (axis == 'x' || axis == 'X')
	{
		Center(examMole, A1);
		ConvertCoord(examMole,'C','2');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z'))-(PI/2));
		ConvertCoord(examMole,'2','C');
		ConvertCoord(examMole,'C','3');
		CoordEdit(examMole,0,0,-(examMole.GetAtomCoord(A2,'z')));
		ConvertCoord(examMole,'3','C');
	}
	return;
}

void Transform::BondRotate(Molecule& examMole, double theta, int atom1, int atom2, std::vector<int>& atomSel)
{
	AlignAxis(examMole, 'x', atom1, atom2);
	ConvertCoord(examMole,'C','1');
	SelectiveEdit(examMole,0,0,theta,atomSel);
	ConvertCoord(examMole,'1','C');
	return;
}

void Transform::DeviceSwitch_ODSP(Molecule& examMole, int * axiAtoms, double angleDisp, double angleInc,std::vector<int>& atomSelRot, std::vector<int>& atomSelRec, std::string filename )
{
	int counter = 1; int i = 0;
	std::stringstream strCount;
	for (double angleTot=angleInc; angleTot<=angleDisp; angleTot=angleTot+angleInc)
	{
		strCount << counter;
		BondRotate(examMole,  angleInc, *(axiAtoms)  , *(axiAtoms+1), atomSelRot);
		BondRotate(examMole, -angleInc, *(axiAtoms+2), *(axiAtoms+3), atomSelRec);
		examMole.PrintInfo(filename + "_" + strCount.str() + ".gamess","gamess");
		counter++; strCount.str("");
	}
	return;
}
