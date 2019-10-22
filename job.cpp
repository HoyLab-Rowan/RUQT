//job.cpp
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include "molecule.h"
#include "transform.h"
#include "electrode.h"
#include "job.h"

void Job::PresentSelection(std::string moleName)
{
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			std::cout << "Molecule: " << moleName << "\n";
			for (int j=0; j<moleVec[i].jobMole.GetAtomCount(); j++)
			{
				std::cout << j << "). " << moleVec[i].jobMole.GetAtomSym(j) << " ";
				std::cout << moleVec[i].jobMole.GetAtomCoord(j,'x') << " ";
				std::cout << moleVec[i].jobMole.GetAtomCoord(j,'y') << " ";
				std::cout << moleVec[i].jobMole.GetAtomCoord(j,'z') << "\n";
			}
			break;
		}
	}
	return;
}
void Job::INIT(std::string filename, std::string filetype, std::string moleName)
{
	moleVec.push_back(Molecules());
	moleVec[moleVec.size()-1].name = moleName;
	moleVec[moleVec.size()-1].jobMole.Initiate(filename, filetype);
	return;
}

void Job::PRINT(std::string filename, std::string filetype, std::string moleName)
{
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			moleVec[i].jobMole.PrintInfo(filename, filetype);
		}
	}
	return;
}

void Job::SORT(std::string moleName, std::string axis)
{
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			char * axischar = new char[1]();
			strcpy(axischar,axis.c_str());
			moleVec[i].jobMole.SortByCoord(*axischar);
		}
	}
}

void Job::ALIGN(std::string moleName, std::string STRaxis)
{
	char axis = STRaxis.at(0); std::string atom1; std::string atom2;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Align phase!\nWhat atoms would you like to select?";
			std::cout << "\nAtom 1 (centered): "; std::cin >> atom1;
			std::cout << "\nAtom 2 (aligned): "; std::cin >> atom2;
			std::cout << "\n\n\n";
			int A1 = atoi(atom1.c_str()); int A2 = atoi(atom2.c_str());
			TF.AlignAxis(moleVec[i].jobMole,axis,A1,A2);
		}
	}
	return;
}

void Job::DEVICE_SWITCH(std::string moleName, std::string angleDisp, std::string angleInc)
{
	char * tokChar; char * charLine;
	std::string atom1, atom2, atom3, atom4, atomlist1, atomlist2, filenamechosen;
	std::vector<int> atomSelRot, atomSelRec;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Device Switch phase!\nLets define rotation 1";
			std::cout << "\nAxis 1; Atom 1: "; std::cin >> atom1;
			std::cout << "\nAxis 1; Atom 2: "; std::cin >> atom2;
			std::cout << "\nAtoms affected? (separate by commas only):"; std::cin >> atomlist1;
			std::cout << "\nNow lets define rotation 2";
			std::cout << "\nAxis 2; Atom 1: "; std::cin >> atom3;
			std::cout << "\nAxis 2; Atom 2: "; std::cin >> atom4;
			std::cout << "\nAtoms affected? (separate by commas only):"; std::cin >> atomlist2;
			std::cout << "\nWhat will the name of the rotated geometries be:"; std::cin >> filenamechosen;
			int axiAtoms[4] = {atoi(atom1.c_str()),atoi(atom2.c_str()),atoi(atom3.c_str()),atoi(atom4.c_str())};
			charLine = new char[atomlist1.length()+1]();
			strcpy(charLine,atomlist1.c_str());
			tokChar = strtok(charLine,",");
			while (tokChar != NULL)
			{
				atomSelRot.push_back(atoi(tokChar));
				tokChar = strtok(NULL,",");
			}
			delete [] charLine;
			charLine = new char[atomlist2.length()+2]();
			strcpy(charLine,atomlist2.c_str());
			tokChar = strtok(charLine,",");
			while (tokChar != NULL)
			{
				atomSelRec.push_back(atoi(tokChar));
				tokChar = strtok(NULL,",");
			}
			TF.DeviceSwitch_ODSP(moleVec[i].jobMole, axiAtoms, atof(angleDisp.c_str()), atof(angleInc.c_str()), atomSelRot, atomSelRec, filenamechosen);
		}
	}
	return;
}

void Job::HYDROREP_ELECTRODE(std::string moleName)
{
	std::string junc1, junc2, Hydro1, Hydro2, MType, A_Mbl, M_Mbl,MCountstr;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Hydrogen replaced electrode phase!\nLets get the first electrode!";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc1;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro1;
			std::cout << "\nNow lets get the second electrode";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc2;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro2;
			std::cout << "\nAtomic symbol of the metal: "; std::cin >> MType;
			std::cout << "\nAtom-Electrode bond length (S-Au = 2.62): "; std::cin >> A_Mbl;
			std::cout << "\nAtom-Electrode bond length (Au-Au = 2.88): "; std::cin >> M_Mbl;
			std::cout << "\n# of Electrode atoms on each electrode: "; std::cin >> MCountstr;
			int juncAtoms[4] = {atoi(junc1.c_str()),atoi(Hydro1.c_str()),atoi(junc2.c_str()),atoi(Hydro2.c_str())};
			double BL1 = atof(A_Mbl.c_str()); double BL2 = atof(M_Mbl.c_str());
			int MCount = atoi(MCountstr.c_str());
			ELE.HydroRepElectrode(moleVec[i].jobMole, MCount, MType, juncAtoms, BL1,BL2);
		} 
	}
	return;
}

void Job::LINEAR_ELECTRODE(std::string moleName)
{
	std::string junc1, junc2, Hydro1, Hydro2, MType, A_Mbl, M_Mbl,MCountstr;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Hydrogen replaced electrode phase!\nLets get the first electrode!";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc1;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro1;
			std::cout << "\nNow lets get the second electrode";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc2;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro2;
			std::cout << "\nAtomic symbol of the metal: "; std::cin >> MType;
			std::cout << "\nAtom-Electrode bond length (S-Au = 2.62): "; std::cin >> A_Mbl;
			std::cout << "\nAtom-Electrode bond length (Au-Au = 2.88): "; std::cin >> M_Mbl;
			std::cout << "\n# of Electrode atoms on each electrode: "; std::cin >> MCountstr;
			int juncAtoms[4] = {atoi(junc1.c_str()),atoi(Hydro1.c_str()),atoi(junc2.c_str()),atoi(Hydro2.c_str())};
			double BL1 = atof(A_Mbl.c_str()); double BL2 = atof(M_Mbl.c_str());
			int MCount = atoi(MCountstr.c_str());
			ELE.LinearElectrode(moleVec[i].jobMole, MCount, MType, juncAtoms, BL1,BL2);
		} 
	}
	return;
}

void Job::TIPSQUARE_ELECTRODE(std::string moleName)
{
	std::string junc1,junc2,MType,A_Mbl,M_Mbl,Hydro1,Hydro2;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Square Tipped electrode phase!\nLets get the first electrode!";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc1;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro1;
			std::cout << "\nNow lets get the second electrode";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc2;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro2;
			std::cout << "\nAtomic symbol of the metal: "; std::cin >> MType;
			std::cout << "\nAtom-Electrode bond length (S-Au = 2.62): "; std::cin >> A_Mbl;
			std::cout << "\nAtom-Electrode bond length (Au-Au = 2.88): "; std::cin >> M_Mbl;
			int juncAtoms[4] = {atoi(junc1.c_str()),atoi(junc2.c_str()),atoi(Hydro1.c_str()),atoi(Hydro2.c_str())};
			double BL1 = atof(A_Mbl.c_str()); double BL2 = atof(M_Mbl.c_str());
			ELE.TipSquareElectrode(moleVec[i].jobMole, MType, juncAtoms, BL1,BL2);
		} 
	}
	return;
}

void Job::PYRAMID_ELECTRODE(std::string moleName)
{
	std::string junc1, junc2, MType, A_Mbl, M_Mbl, layerIter, Hydro1, Hydro2;
	for (int i=0; i<moleVec.size(); i++)
	{
		if (moleName == moleVec[i].name)
		{
			PresentSelection(moleVec[i].name);
			std::cout << "\nMoleMod is in the Pyramid electrode phase!\nLets get the first electrode!";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc1;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro1;
			std::cout << "\nNow lets get the second electrode";
			std::cout << "\nAtom bonded to gold: "; std::cin >> junc2;
			std::cout << "\nHydrogen to be replaced: "; std::cin >> Hydro2;
			std::cout << "\nAtomic symbol of the metal: "; std::cin >> MType;
			std::cout << "\nHow many layers will the pyramid have:"; std::cin >> layerIter;
			std::cout << "\nAtom-Electrode bond length (S-Au = EDIT): "; std::cin >> A_Mbl;
			std::cout << "\nAtom-Electrode bond length (Au-Au = EDIT): "; std::cin >> M_Mbl;
			int juncAtoms[4] = {atoi(junc1.c_str()),atoi(junc2.c_str()),atoi(Hydro1.c_str()),atoi(Hydro2.c_str())};
			int layers = atoi(layerIter.c_str());
			double BL1 = atof(A_Mbl.c_str()); double BL2 = atof(M_Mbl.c_str());
			ELE.PyramidElectrode(moleVec[i].jobMole, layers, MType, juncAtoms, BL1,BL2);
		} 
	}
	return;
}

void Job::run(std::string mmFile)
{
	std::cout << "Currently working with: " << mmFile << "\n";
	std::string strLine;
	std::string command; char * tokChar;
	char * charLine;
	int strlength;
	std::ifstream inputFile(mmFile.c_str());
	
	if (!inputFile){std::cout << mmFile << " was not found\n"; return;}
	
	while (getline(inputFile, strLine))
	{
		strlength = strLine.length();
		charLine = new char[strlength+1]();
		strcpy(charLine, strLine.c_str());
		tokChar = strtok(charLine," "); command = tokChar;
		if (command == "INIT")
		{
			tokChar = strtok(NULL," "); std::string INIT_1 = tokChar;
			tokChar = strtok(NULL," "); std::string INIT_2 = tokChar;
			tokChar = strtok(NULL," "); std::string INIT_3 = tokChar;
			INIT(INIT_1, INIT_2, INIT_3);
		}
		if (command == "PRINT")
		{
			tokChar = strtok(NULL," "); std::string PRINT_1 = tokChar;
			tokChar = strtok(NULL," "); std::string PRINT_2 = tokChar;
			tokChar = strtok(NULL," "); std::string PRINT_3 = tokChar;
			PRINT(PRINT_1, PRINT_2, PRINT_3);
		}
		if (command == "ALIGN")
		{
			tokChar = strtok(NULL," "); std::string ALIGN_1 = tokChar;
			tokChar = strtok(NULL," "); std::string ALIGN_2 = tokChar;
			ALIGN(ALIGN_1, ALIGN_2);
		}
		if (command == "DEVICE_SWITCH")
		{
			tokChar = strtok(NULL," "); std::string DS_1 = tokChar;
			tokChar = strtok(NULL," "); std::string DS_2 = tokChar;
			tokChar = strtok(NULL," "); std::string DS_3 = tokChar;
			DEVICE_SWITCH(DS_1,DS_2,DS_3);
		}
		if (command == "HR_ELEC")
		{
			tokChar = strtok(NULL," "); std::string HRE_1 = tokChar;
			HYDROREP_ELECTRODE(HRE_1);
		}
		if (command == "L_ELEC")
		{
			tokChar = strtok(NULL," "); std::string LE_1 = tokChar;
			LINEAR_ELECTRODE(LE_1);
		}
		if (command == "SORT")
		{
			tokChar = strtok(NULL," "); std::string S_1 = tokChar;
			tokChar = strtok(NULL," "); std::string S_2 = tokChar;
			SORT(S_1,S_2);
		}
		if (command == "TSQ_ELEC")
		{
			tokChar = strtok(NULL," "); std::string TSQE_1 = tokChar;
			TIPSQUARE_ELECTRODE(TSQE_1);
		}
		if (command == "PYR_ELEC")
		{
			tokChar = strtok(NULL," "); std::string PYRE_1 = tokChar;
			PYRAMID_ELECTRODE(PYRE_1);
		}
		delete charLine;
	}
	return;
}
