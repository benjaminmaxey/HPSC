###############################################################################
#	Author: Ben Maxey
#
#	Compiles all source code for the programs vandermonde, test_newton, and
#	kepler.
###############################################################################

STD = -std=c++11

all: vandermonde test_newton kepler

vandermonde: vandermonde.o matrix.o
	g++ $^ -o $@ $(STD)

test_newton: test_newton.o newton.o
	g++ $^ -o $@ $(STD)

kepler: kepler.o newton.o matrix.o
	g++ $^ -o $@ $(STD)

vandermonde.o: vandermonde.cpp
	g++ -c vandermonde.cpp $(STD)

matrix.o: matrix.cpp
	g++ -c matrix.cpp $(STD)

test_newton.o: test_newton.cpp
	g++ -c test_newton.cpp $(STD)

newton.o: newton.cpp
	g++ -c newton.cpp $(STD)

kepler.o: kepler.cpp
	g++ -c kepler.cpp $(STD)

clean:
	rm -f *o
	rm -f *txt
	rm vandermonde
	rm test_newton
	rm kepler