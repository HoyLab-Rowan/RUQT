molemod : clean main.o molecule.o transform.o electrode.o job.o
	g++ -o molemod main.o molecule.o transform.o electrode.o job.o

main.o :
	g++ -c main.cpp

molecule.o : molecule.cpp molecule.h
	g++ -c molecule.cpp

transform.o : transform.cpp transform.h
	g++ -c transform.cpp

electrode.o : electrode.cpp electrode.h
	g++ -c electrode.cpp

job.o : job.cpp job.h
	g++ -c job.cpp

clean :
	rm -f *.o
