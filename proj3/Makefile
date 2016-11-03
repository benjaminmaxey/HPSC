FLAG = -std=c++14

all: test_Newtonform lagrange2D runge_uniform runge_chebyshev

test_Newtonform: test_Newtonform.o newton.o matrix.o
	g++ $^ $(FLAG) -o $@

lagrange2D: test_Lagrange2D.o matrix.o
	g++ $^ $(FLAG) -o $@

runge_uniform: Runge_uniform.o Lagrange2D.o matrix.o
	g++ $^ $(FLAG) -o $@

runge_chebyshev: Runge_Chebyshev.o Lagrange2D.o matrix.o
	g++ $^ $(FLAG) -o $@

%.o: %.cpp
	g++ -c $< $(FLAG)

clean:
	rm -f *.txt *.o test_Newtonform lagrange2D runge_uniform runge_chebyshev