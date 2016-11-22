###################################################################
#	Author: Ben Maxey
#	
#	Makefile for HPSC Project 4
###################################################################

#Flags
FLAG = -std=c++14

#Targets
all : test_int.exe test_adapt.exe

test_int.exe : test_int.o
	g++ $^ -o $@ $(FLAG)

test_adapt.exe : test_adapt.o
	g++ $^ -o $@ $(FLAG) 

%.o: %.cpp
	g++ -c $< $(FLAG)

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.o *.txt *.exe