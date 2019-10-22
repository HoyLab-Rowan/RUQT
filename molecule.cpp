//molecule.cpp
#include "molecule.h"
#include <cstring>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cstdlib>
void Molecule::Initiate(std::string filename, std::string filetype)
{
	
	if (filetype == "xyz")
	{
		std::string strLine;
		std::string stringtok; //needed to read tokLine
		char * tokLine;
		char * charLine;
		int counterLine=0;
		int counterTok=0;
		int strlength;
		std::ifstream inputFile(filename.c_str());
		
		if (!inputFile){std::cout << "File did not open\n"; exit(1);}
		while (getline(inputFile, strLine))
		{
			if (counterLine == 0) {atomCount = atoi(strLine.c_str());}
			else if (counterLine > 1)
			{
				strlength = strLine.length();
				charLine = new char[strlength+1]();
				strcpy(charLine, strLine.c_str());
				
				atomVec.push_back(atom());
				tokLine = strtok(charLine," ");
				
				while (tokLine != NULL)
				{
					stringtok = tokLine;
					if (counterTok == 0){atomVec[counterLine-2].atomSym = stringtok;}
					if (counterTok == 1){atomVec[counterLine-2].coord[0] = atof(tokLine);}
					if (counterTok == 2){atomVec[counterLine-2].coord[1] = atof(tokLine);}
					if (counterTok == 3){atomVec[counterLine-2].coord[2] = atof(tokLine);}
					tokLine = strtok (NULL, " ");
					counterTok++;
				}
				counterTok = 0;
				delete charLine;
			}
			counterLine++;
		}
		inputFile.close();
	}
	return;
}

double Molecule::GetAtomNum(std::string sym)
{
	if (sym == "H"){return 1.0;}
	if (sym == "C"){return 6.0;}
	if (sym == "N"){return 7.0;}
	if (sym == "O"){return 8.0;}
	if (sym == "F"){return 9.0;}
	if (sym == "Al" || sym == "AL"){return 13.0;}
	if (sym == "Si" || sym == "SI"){return 14.0;}
	if (sym == "S"){return 16.0;}
	if (sym == "Cl" || sym == "CL"){return 17.0;}
	if (sym == "Au" || sym == "AU"){return 79.0;}
	return 0;
}

void Molecule::PrintInfo(std::string filename, std::string filetype)
{
	if (filetype == "xyz")
	{
		std::ofstream outputFile;
		outputFile.open(filename.c_str());
		if (outputFile)
		{
			outputFile << atomCount << "\n\n";
			for (int i = 0; i < atomCount; i++)
			{
				outputFile << std::fixed;//keep in loop?
				outputFile << atomVec[i].atomSym << "    ";
				outputFile << std::setprecision(9) << atomVec[i].coord[0] 
					<< "    " << atomVec[i].coord[1] << "    "
				       	<< atomVec[i].coord[2] << '\n';
			}
		}
	}
	if (filetype == "gamess")
	{
		std::ofstream outputFile;
		outputFile.open(filename.c_str());
		if (outputFile)
		{
			for (int i = 0; i < atomCount; i++)
			{
				outputFile << std::fixed;//keep in loop?
				outputFile << atomVec[i].atomSym << "    ";
				outputFile << std::setprecision(1) << GetAtomNum(atomVec[i].atomSym) << "  ";
				outputFile << std::setprecision(9) << atomVec[i].coord[0] 
					<< "    " << atomVec[i].coord[1] << "    "
				       	<< atomVec[i].coord[2] << '\n';
			}
		}
	}
	return;
}

void Molecule::SetAtomCoord(int ithAtom,char toEdit, double coordValue)
{
	if(toEdit == 'x'){atomVec[ithAtom].coord[0] = coordValue;}
	else if(toEdit == 'y'){atomVec[ithAtom].coord[1] = coordValue;}
	else if(toEdit == 'z'){atomVec[ithAtom].coord[2] = coordValue;}
	return;
}

void Molecule::SetAtomSym(int ithAtom, std::string newAtom){atomVec[ithAtom].atomSym = newAtom;}

double Molecule::GetAtomCoord(int ithAtom, char toGet)
{
	if (toGet == 'x'){return atomVec[ithAtom].coord[0];}
	else if (toGet == 'y'){return atomVec[ithAtom].coord[1];}
	else if (toGet == 'z'){return atomVec[ithAtom].coord[2];}
	return 0.0;
}

std::string Molecule::GetAtomSym(int ithAtom){return atomVec[ithAtom].atomSym;}

void Molecule::AddAtom(std::string AtSym, double coordX, double coordY, double coordZ)
{
	const int CURRSIZE = atomVec.size();
	atomVec.push_back(atom());
	atomVec[CURRSIZE].atomSym = AtSym;
	atomVec[CURRSIZE].coord[0] = coordX;
	atomVec[CURRSIZE].coord[1] = coordY;
	atomVec[CURRSIZE].coord[2] = coordZ;
	atomCount += 1;
	return;
}

void Molecule::RemoveAtom(int ithAtom)
{
	atomVec.erase(atomVec.begin()+ithAtom);
	atomCount -= 1;
	return;
}
int Molecule::GetAtomCount() {return atomCount;}

void Molecule::SwapStruct(atom& atom1, atom& atom2)
{
	atom tmpatom;
	tmpatom.atomSym  = atom1.atomSym;
	tmpatom.coord[0] = atom1.coord[0];
	tmpatom.coord[1] = atom1.coord[1];
	tmpatom.coord[2] = atom1.coord[2];
	atom1.atomSym  = atom2.atomSym;
	atom1.coord[0] = atom2.coord[0];
	atom1.coord[1] = atom2.coord[1];
	atom1.coord[2] = atom2.coord[2];
	atom2.atomSym  = tmpatom.atomSym;
	atom2.coord[0] = tmpatom.coord[0];
	atom2.coord[1] = tmpatom.coord[1];
	atom2.coord[2] = tmpatom.coord[2];
}

void Molecule::SortByCoord(char axis)
{
	int i, j, a;
	double tmpval;
	switch (axis)
	{
		case 'x': a=0; break;
		case 'X': a=0; break;
		case 'y': a=1; break;
		case 'Y': a=1; break;
		case 'z': a=2; break;
		case 'Z': a=2; break;
		default: break;
	}
	for (i=1;i<atomVec.size();i++)
	{
		j=i;
		tmpval=atomVec[i].coord[a];
		while (j>0 && tmpval<atomVec[j-1].coord[a])
		{
			SwapStruct(atomVec[j],atomVec[j-1]);
			j--;
		}
	}
}
