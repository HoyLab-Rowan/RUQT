#The following lines are for gfortran compilation with MKL libraries
#FXX = gfortran -o2 -ffast-math -ffree-line-length-none -m64 #-fbacktrace
#CXX = g++ -o2
#MKLROOT     = /csm_data/hoy_group/intel/mkl
#BLAS     =   -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
#INC      = -I${MKLROOT}/include
#LIB       = ${BLAS}

#ifort commands(if -mkl does not work, try explicit links above)
FXX       = ifort -O2 -axCORE-AVX2 -parallel -fPIC -mkl:parallel -static-intel -qopenmp-link static
CXX       = g++

CWD       = $(shell pwd)
OBJ       = $(CWD)/Object
JMOD      = $(CWD)/JunctionMod
EXEC      = RUQT_v1r.x


main: $(EXEC) junctionmod

$(EXEC): RUQT.f90 ${OBJ}/TypeMod.o  ${OBJ}/InterfaceMod.o ${OBJ}/Build_G_SD_Invert.o ${OBJ}/Build_B0_CISD.o

	$(FXX) ${INC} ${LIB} RUQT.f90 ${OBJ}/TypeMod.o ${OBJ}/InterfaceMod.o ${OBJ}/Build_G_SD_Invert.o ${OBJ}/Build_B0_CISD.o -o ${EXEC}

junctionmod: ${OBJ}/main.o ${OBJ}/molecule.o ${OBJ}/transform.o ${OBJ}/electrode.o ${OBJ}/job.o
	$(CXX) -o junctionmod ${OBJ}/main.o ${OBJ}/molecule.o ${OBJ}/transform.o ${OBJ}/electrode.o ${OBJ}/job.o

${OBJ}/main.o : ${JMOD}/main.cpp
	$(CXX) -c ${JMOD}/main.cpp -o ${OBJ}/main.o

${OBJ}/molecule.o : ${JMOD}/molecule.cpp ${JMOD}/molecule.h
	$(CXX) -c ${JMOD}/molecule.cpp -o ${OBJ}/molecule.o

${OBJ}/transform.o : ${JMOD}/transform.cpp ${JMOD}/transform.h
	$(CXX) -c ${JMOD}/transform.cpp -o ${OBJ}/transform.o

${OBJ}/electrode.o : ${JMOD}/electrode.cpp ${JMOD}/electrode.h
	$(CXX) -c ${JMOD}/electrode.cpp -o ${OBJ}/electrode.o

${OBJ}/job.o : ${JMOD}/job.cpp ${JMOD}/job.h
	$(CXX) -c ${JMOD}/job.cpp -o ${OBJ}/job.o

${OBJ}/TypeMod.o: TypeMod.f90
	$(FXX) -c TypeMod.f90 -o ${OBJ}/TypeMod.o

${OBJ}/InterfaceMod.o: InterfaceMod.f90
	$(FXX) -c InterfaceMod.f90 -o ${OBJ}/InterfaceMod.o

${OBJ}/Build_G_SD_Invert.o: Build_G_SD_Invert.f90
	$(FXX) -c Build_G_SD_Invert.f90 -o ${OBJ}/Build_G_SD_Invert.o

${OBJ}/Build_B0_CISD.o: Build_B0_CISD.f90
	$(FXX) -c Build_B0_CISD.f90 -o ${OBJ}/Build_B0_CISD.o

clean:
	rm $(OBJ)/*.o; rm *.mod
