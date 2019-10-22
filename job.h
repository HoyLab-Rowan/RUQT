//job.h
#ifndef JOB_H
#define JOB_H
#include <iostream>
#include <cstring>
#include <vector>
#include "molecule.h"
#include "transform.h"
#include "electrode.h"

class Job
{
	struct Molecules {
		Molecule jobMole;
		std::string name;
	};
	std::vector<Molecules> moleVec;
	Transform TF;
	Electrode ELE;
	
	void PresentSelection(std::string);
	void INIT(std::string,std::string,std::string);
	void PRINT(std::string,std::string,std::string);
	void SORT(std::string,std::string);
	void ALIGN(std::string,std::string);
	void DEVICE_SWITCH(std::string,std::string,std::string);
	void HYDROREP_ELECTRODE(std::string);
	void TIPSQUARE_ELECTRODE(std::string);
	void PYRAMID_ELECTRODE(std::string);
	void LINEAR_ELECTRODE(std::string);
public:
	void run(std::string);//only public function
};
#endif
