#include "CPhys.h"
#include <time.h>
#include <iostream>
int main(int argc, const char *argv[])
{
	// TODO make unit test class
	// TODO implement move constructors
	// TODO implement more "saftey" functions for 
	// the functions acessible by the user
	int    N = 900;
	Matrix A1 = Matrix(N,N);
	double** pA1 = A1.getArrayPointer();

	Matrix A2 = Matrix(N,N);
	double** pA2 = A2.getArrayPointer();

	Vector A3 = Vector(N*N);
	double* pA3 = A3.getArrayPointer();


	clock_t start1 = clock(); 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			pA1[i][j] = j + 300;
		}
	}
	clock_t end1 = clock();

	clock_t start2 = clock(); 
	/* for (int i = 0; i < N; i++) { */
	/* 	for (int j = 0; j < N; j++) { */
	/* 		A2(i,j) = i+j; */
	/* 	} */
	/* } */
	for (int i = 0; i < N*N; i++) {
		pA3[i] = i + 300;
	}
	clock_t end2 = clock();

	double time1 = double(end1 - start1)/CLOCKS_PER_SEC;
	double time2 = double(end2 - start2)/CLOCKS_PER_SEC;
	std::cout << "Time1 = " << time1 << std::endl;
	std::cout << "Time2 = " << time2 << std::endl;
	std::cout << "Ratio = " << time2/time1 << std::endl;

	N = 2;
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
