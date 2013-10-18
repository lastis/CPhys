#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "CPhys.h"

using namespace std;
using namespace CPhys;

void EigVal::jacobiMethod(Matrix& A, Matrix& R, int N){
	using namespace mEigVal;
	double** ppR = R.getArrayPointer();
	double** ppA = A.getArrayPointer();
	int 	 k, l;
	int 	 iterations  	= 0;
	double	 eps 		= 10e-8;
	double 	 maxIterations 	= (double) N * N * N;
	// Prepare eigenvector matrix
	R.eye();
	// Decide which elements to rotate
	double 	maxoff		= maxoffdiag(ppA, &k, &l, N);
	// Do the rotations until we converge at a solution
	while (		maxoff > eps 	
		 && iterations < maxIterations ) {
		
		int num = iterations;
		if(num%20000 == 0){
			cout << "Computing - "<< iterations << " iterations ";
			cout << " max off diag = " << maxoff << endl;
		}
		rotate(ppA, ppR, k, l, N );
		maxoff = maxoffdiag(ppA, &k, &l, N );
		iterations++;
	}
	cout << "Number of iterations: " << iterations << "\n";
}

double VecOp::normalize(Vector& v, double dx){
	int     N   = v.getLength();
	double* pV  = v.getArrayPointer();
	double  sum = 0;
	double  Z   = 0;
	for (int i = 0; i < N; i++) {
		sum += pV[i]*pV[i]*dx;
	}
	Z = 1/sqrt(sum);
	for (int i = 0; i < N; i++) {
		pV[i] *= Z;
	}
	return Z;

}

void	MatOp::sortCol(Matrix& A, Vector& v){
	// We want to sort the colums in A by the elements in v.
	int rows = A.getN();
	int cols = A.getM();
	// Make a new matrix wit dim N = A.N+1 and M = A.M
	Matrix tmp = Matrix(rows+1,cols);
	// Insert v in the first row
	tmp.setRow(0,v);
	// Copy A to tmp after the first row
	double** ppTmp = tmp.getArrayPointer();
	double** ppA   =    A.getArrayPointer();
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ppTmp[i+1][j] = ppA[i][j];
		}
	}
	// Our matrix is row major so we need to transpose
	// to be able to sort by columns
	// Now **arr[0,1,2, ... , i] points to whole arrays 
	// of the columns
	tmp.t();
	// Sort by first elements of A (sorting by vector v)
	ppTmp = tmp.getArrayPointer();
	qsort((void*)ppTmp,cols,sizeof(double),&mMatOp::compareTwoRows);
	// Transpose back
	tmp.t();
	ppTmp = tmp.getArrayPointer();
	// Now copy the sorted elements to A
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			 ppA[i][j] = ppTmp[i+1][j];
		}
	}
	// Then sort v
	v.sort();
}

namespace{

	
	int 	mMatOp::compareTwoRows(const void* rowA, const void* rowB){
		return (**(double**) rowA - **(double**) rowB);
	}

	double 	mEigVal::maxoffdiag(double** A, int* k, int* l, int N){
		double max  = 0.0;
		double a_ij = 0.0;
		for (int i = 0; i < N; i++) {
			for (int j = i + 1; j < N; j++) {
				a_ij = fabs(A[i][j]);
				if (a_ij > max){
					max	= a_ij;
					*l 	= i;
					*k 	= j;
				}
			}
		}
		return max;
	}

	void mEigVal::rotate(double** A, double** R,int k, int l, int N){
		double s, c;
		double tau, t;
		tau = (A[l][l] -A[k][k])/(2*A[k][l]);
		if(tau > 0) t = -tau + sqrt(1+tau*tau);
		else 	    t = -tau - sqrt(1+tau*tau);
		c = 1 / sqrt(1+t*t);
		s = t * c;
		double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
		a_kk = A[k][k];
		a_ll = A[l][l];
		// changing the matrix elements with indices k and l
		A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
		A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
		A[k][l] = 0.0; // hard-coding of the zeros
		A[l][k] = 0.0;
		// and then we change the remaining elements
		for ( int i = 0; i < N; i++ ) {
			if ( i != k && i != l ) {
				a_ik 	= A[i][k];
				a_il 	= A[i][l];
				A[i][k]	= c*a_ik - s*a_il;
				A[k][i]	= A[i][k];
				A[i][l]	= c*a_il + s*a_ik;
				A[l][i]	= A[i][l];
			}
			// Finally, we compute the new eigenvectors
			r_ik 	= R[i][k];
			r_il 	= R[i][l];
			R[i][k]	= c*r_ik - s*r_il;
			R[i][l]	= c*r_il + s*r_ik;
		}
	}
}

