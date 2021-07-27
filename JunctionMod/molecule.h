// molecule.h
#ifndef MOLECULE_H
#define MOLECULE_H
#include <iostream>
#include <string>
#include <vector>
class Molecule
{
	int atomCount;
	struct atom {
		std::string atomSym;	
		double coord[3];
	};
	std::vector<atom> atomVec; //Lists all the atoms in no specified order (from file)
	void SwapStruct(atom&, atom&); //swap position of two atoms in previous vector ^
	double GetAtomNum(std::string); //Gets atomic number of atom
public:
	void Initiate(std::string, std::string); //
	void PrintInfo(std::string, std::string);
	void SetAtomCoord(int,char,double);
	void SetAtomSym(int, std::string);
	double GetAtomCoord(int, char);
	std::string GetAtomSym(int);
	void AddAtom(std::string, double, double, double);
	void RemoveAtom(int);
	int GetAtomCount();
	void SortByCoord(char);
};
#endif
