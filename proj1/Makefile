all: proj1_a proj1_b

proj1_a: proj1_a.o nest.o matrix.o
	g++ $^ -o $@

proj1_b: proj1_b.o matrix.o
	g++ $^ -o $@

proj1_a.o: proj1_a.cpp
	g++ -c proj1_a.cpp

proj1_b.o: proj1_b.cpp
	g++ -c proj1_b.cpp

nest.o: nest.cpp
	g++ -c nest.cpp

matrix.o: matrix.cpp
	g++ -c matrix.cpp

clean:
	rm -f *o
	rm -f *txt