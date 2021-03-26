# RUQT

***Rowan University Quantum Transport (RUQT)***

RUQT is an in-progress non-equilibrium Green's function (NEGF) based software package focused on studying electron correlation effects in charge transport problems. It is focused on integrating N- and 2-electron methods into the NEGF transport formalism and is the home of the NEGF-RDM transport method. As of March 2021, it is capable of performing Landauer NEGF calculations using electronic structure data acquired from separate Hartree-Fock, Density Functional Theory, and Parametric 2-Electron Reduced Density Matrix Theory calculations. The code itself is written in Fortran90+ while the supporting scripts are either in python or Maple. The recommended supporting code JunctionMod is included and is compiled by default.

RUQT currently can read data from the following electronic structure software packages (more methods/package support to come in future):

HF and DFT: Q-Chem, PYSCF, GAMESS, Maple Quantum Chemistry Toolbox, and Molcas
p2-RDM: Maple Quantum Chemistry and GAMESS (Not publically available)

RUQT v1.0 is capable of performing non-self-consistent current and transmission calculations with metal wide band limit electrodes (fixed Fermi level). 
Additional electrode and calculation types to be added later in 2021.

***Installation***

Required Libraries: Intel MKL

Compiler: GNU or ifort

1. Download the code and open the Makefile. 
2. Set the MKLROOT variable to the location of your MKL library. 
3. Compiles using Intel fortran compilers by default. Make sure that ifort can find your MKL libraries using the -mkl flag
3. To use gfortran uncomment the gfortran instead line with MKL library links and comment out the ifort line.
4. Close Makefile and enter the make command in the terminal.
5. This should produce an "RUQT_1r.x" executable. 
6. If needed, use "make clean" to clear all .mod and .o files from the source and object directories.

***How to Use***

To run calculations, use either the Maple or Python scripts located in the Hoy Research Group Github: https://github.com/HoyLab-Rowan. 

The Maple scripts are available for all methods, and Python scripts are available for Q-Chem (HF & DFT), PYSCF (HF & DFT), GAMESS (HF), and Molcas (DFT) calculations. 

The Maple scripts only require an .xyz file with the junction (molecule+electrode) geometry. These can be produced by the JunctionMod code is included in this package and can be used to construct/align molecular junctions. If you use another software package to generate the xyz make sure that the junction must be placed along the x-axis. 

Using the python scripts requires forming separate electronic structure calculations with specific keywords in the inputs in order to print the Fock, overlap, and T1/T2+MO energies needed by the NEGF-RDM code. Example inputs and outputs for Maple and GAMESS (HF & p2-RDM) can be found in the /examples folder of the source code.

For running p2-RDM calculations, using the Maple Quantum Chemistry Toolbox and related scripts is highly recommended.

The GAMESS p2-RDM code is acquired from the Mazziotti group at UChicago and is not publically available at this time. (Email Dr. Hoy at hoy@rowan.edu if you need access to the GAMESS p2-RDM code). 

***If you use this code, cite:***
1. Erik P. Hoy, David A. Mazziotti, and Tamar Seideman, “Development and application of a 2-electron reduced density matrix approach to electron transport via molecular junctions” J. Chem. Phys. 147, 184110 (2017).
