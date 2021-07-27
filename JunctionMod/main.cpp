//main.cpp
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <vector>
#include "molecule.h"
#include "transform.h"
#include "electrode.h"
#include "job.h"
void PrintHelp()
{
std::cout << "Init: Collects data from file to create a molecule in code\n\n";
std::cout << "INIT fileName fileType moleculeName\n\n";
std::cout << "do not include any special characters except _ in moleculeName,\n";
std::cout << "and no spaces. The moleculeName is arbitrary, but must remain\n";
std::cout << "consistent. Only filetype available is xyz.\n\n";
std::cout << "Example: INIT benzene.xyz xyz ben\n\n";
std::cout << "----------------------------------------------------------------\n\n";
std::cout << "Print: Prints out current coordinates of molecule in code\n\n";
std::cout << "PRINT filename fileType moleculeName\n\n";
std::cout << "The molecule name must be exactly the way it was Initiated\n";
std::cout << "It is reccomended to choose a new file name to keep the original\n";
std::cout << "geometry. fileTypes available are xyz and gamess\n\n";
std::cout << "Example: INIT benzene.xyz xyz ben\n\n";
std::cout << "----------------------------------------------------------------\n\n";
std::cout << "Sort:\n\n";
std::cout << "SORT moleculeName axis\n\n";
std::cout << "Changes the order of atoms in code to be organized by spatial\n";
std::cout << "arrangements. Very handy when a command triggers a display.\n\n";
std::cout << "Example: SORT ben x\n\n";
std::cout << "----------------------------------------------------------------\n\n";
std::cout << "Align:\n\n";
std::cout << "ALIGN moleculeName axis\n\n";
std::cout << "Translates and rotates molecule until the first atom is 0,0,0\n";
std::cout << "and the second atom intersects the desired axis. Triggers a\n";
std::cout << "display, for additional input during program.\n\n";
std::cout << "Example: ALIGN ben x\n\n";
std::cout << "----------------------------------------------------------------\n\n";
std::cout << "Device Switch:\n\n";
std::cout << "DEVICE_SWITCH moleculeName angularDisplacement angularIncrement\n\n";
std::cout << "Rotates a desired bond in the molecule for analysis of different\n";
std::cout << "conformations. The angular Displacement represents how much the\n";
std::cout << "is rotated in total, and the angular Increment is the difference\n";
std::cout << "the dihedral angle changes with each molecule. This function\n";
std::cout << "prints results without the PRINT command.\n\n";
std::cout << "Example: DEVICE_SWITCH ben .785 .157\n\n";
std::cout << "----------------------------------------------------------------\n\n";
std::cout << "Electrode Family:\n\n";
std::cout << "???_ELEC moleculeName\n\n";
std::cout << "The ??? that preceeds _ELEC depends on the desired electrode.\n";
std::cout << "Prefixes include: HR, Hydrogen Replaced; L, Linear;\n";
std::cout << "TSQ, Tip SQuare; PYR, PYRamid.\n";
std::cout << "This input type only requires the moleculeName, then requests\n";
std::cout << "the additional input later.\n";
std::cout << "Example: L_ELEC ben\n\n";
return;
}

int main(int argc, char **argv){
//
	if (argc == 1)
	{
		std::cout << "Too few arguments, include .mm file\n";
		std::cout << "If you would like a list of JunctionMod Line Inputs\nthen type 'junctionmod help'\n";
		std::cout << "If you would like to keep this list of functions\nthen type 'junctionmod help > FILENAME.txt'\n";
		return 1;
		}

	else 
	{
		Job jobExecuter;
		std::string argument = argv[1];
	
              
                 std::cout << "        __                 __  _                __  _______  ____ \n";
                 std::cout << "       / /_  ______  _____/ /_(_)___  ____     /  |/  / __ \/ __ \ \n";
                 std::cout << "  __  / / / / / __ \/ ___/ __/ / __ \/ __ \   / /|_/ / / / / / / / \n";
                 std::cout << " / /_/ / /_/ / / / / /__/ /_/ / /_/ / / / /  / /  / / /_/ / /_/ /  \n";
                 std::cout << " \____/\__,_/_/ /_/\___/\__/_/\____/_/ /_/  /_/  /_/\____/_____/  \n";
                                                                 

                                                                 

		if (argument == "help"){PrintHelp(); exit(1);}
		std::cout << "          ************************STARTING MODIFICATIONS************************          \n\n\n";
		for (int i=1; i<argc; i++)
		{
			argument = argv[i];
			jobExecuter.run(argument);
		}
	}
	std::cout << "\n          ************************MODIFICATION COMPLETED************************          \n\n";
//
	return 0;
}
