
objects = humanization_main.o Deimmunization.o 

humanization.out : $(objects)
	g++ -o humanization.out $(objects)
humanization_main.o : humanization_main.cpp Deimmunization.h
	g++ -I/usr/global/boost/1_51_0/include -c humanization_main.cpp -o humanization_main.o
/usr/global/boost/1_51_0/include Deimmunization.o : Deimmunization.cpp Deimmunization.h
	g++ -c Deimmunization.cpp -o Deimmunization.o

clean :
	rm -f $(objects)
