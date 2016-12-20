

objects = refinement_main.o InitialAntigen.o
refinement.out : $(objects)
	#g++ -L/usr/global/boost/1_51_0/lib -o refinement.out $(objects)
	g++ -o refinement.out $(objects)
refinement_main.o : refinement_main.cpp InitialAntigen.h MOLECULES.h
	g++ -I/usr/global/boost/1_51_0/include -c refinement_main.cpp -o refinement_main.o -std=c++0x
InitialAntigen.o : InitialAntigen.cpp InitialAntigen.h MOLECULES.h ENERGY_FUNCTIONS.h DOCKING_FUNCTIONS.h
	g++  -I/usr/global/boost/1_51_0/include -c InitialAntigen.cpp -o InitialAntigen.o -std=c++0x

clean :
	rm -f $(objects)
