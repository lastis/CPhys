#include "CPhys.h"
int main(int argc, const char *argv[])
{
	// TODO make unit test class
	// TODO implement move constructors
	// TODO implement more "saftey" functions for 
	// the functions acessible by the user
	int    N = 2;
	Vector a = Vector(N);
	Vector b = Vector(N);
	Matrix A = Matrix(N,N);
	Matrix B = Matrix(N,N);
	a = 2;
	b = 3;
	a = b;
	A.diag(a,-2);

	b = A.diag(-2);

	A = 4;
	B = 2;

	A += 2;
	A.print();
	
	return 0;
}
