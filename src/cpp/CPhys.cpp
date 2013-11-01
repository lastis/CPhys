#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include "CPhys.h"


using namespace std;
using namespace CPhys;

void	ODE::solveEuler(Matrix&  u, Vector& 	   uInit, 
		        double  dt, ODE_Interface*   ode){
	// Call a more general solve method
	using namespace mODE;
	solve(u,uInit,dt,ode,&stepEuler);
}

void	ODE::solveRK4(Matrix&  u, Vector& 	 uInit, 
		      double  dt, ODE_Interface*   ode){
	// Call more general solve method
	using namespace mODE;
	solve(u,uInit,dt,ode,&stepRK4);
}

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



double 	VecOp::normalize(Vector& v, double dx){
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

	void	mODE::solve(Matrix&  u, Vector&        uInit,  
			    double& dt, ODE_Interface*   ode,
		      void (*stepFunc) (double*,double,double*       ,
			                double&,int&  ,ODE_Interface*)){
		// Number of equations
		int 	 eq    = u.getM();
		int	 N     = u.getN();
		double** ppU   = u.getArrayPointer();
		double*   pU   = new double[eq];
		double*   pV    = new double[eq];
		double*   pUOut = new double[eq];

		// Insert the inital values
		for (int i = 0; i < eq; i++) {
			pU[i]   = uInit(i);
			ppU[0][i] = pU[i];
		}
		// Start calculating with the given step function
		for (int i = 0; i < N-1; i++) {
			// Paramaters are;
			// uIn, current time, vIn, uOut, dt
			// the number of equations and the object with the
			// derivative function. 
			stepFunc(pU,i*dt,pUOut,dt,eq,ode);
			// We have now calculated step i+1
			for (int j = 0; j < eq; j++) {
				ppU[i+1][j] = pUOut[j];
				pU[j]     = pUOut[j];
			}
		}
	}


	void 	mODE::stepEuler(double* pU, double   t, double*        pUOut, 
			        double& dt, int&    eq, ODE_Interface*   ode){
		double* pV  = new double[eq];
		// Calculate pV
		ode->derivative(pU,t,pV);
		for (int i = 0; i < eq; i++) {
			pUOut[i] = pU[i] + pV[i]*dt;
		}
	}
	
	void 	mODE::stepRK4(double* pU, double   t, double* 	     pUOut, 
			      double& dt, int&    eq, ODE_Interface*   ode){
		double* pV  = new double[eq];
		double* k1  = new double[eq];
		double* k2  = new double[eq];
		double* k3  = new double[eq];
		double* k4  = new double[eq];
		double* tmp = new double[eq];
		double dt2 = dt/2;
		// The derivative functions calculates pV
		// Calculate k2
		ode->derivative(pU,t,pV);
		for (int i = 0; i < eq; i++) {
			k1[i] = pV[i]*dt;
		}
		// Calculate k2
		for (int i = 0; i < eq; i++) {
			tmp[i] = pU[i] + 0.5*k1[i];
		}
		ode->derivative(tmp,t+dt2,pV);
		for (int i = 0; i < eq; i++) {
			k2[i] = pV[i]*dt;
		}
		// Calculate k3
		for (int i = 0; i < eq; i++) {
			tmp[i] = pU[i] + 0.5*k2[i];
		}
		ode->derivative(tmp,t+dt2,pV);
		for (int i = 0; i < eq; i++) {
			k3[i] = pV[i]*dt;
		}
		// Calculate k4
		for (int i = 0; i < eq; i++) {
			tmp[i] = pU[i] + k3[i];
		}
		ode->derivative(tmp,t+dt,pV);
		for (int i = 0; i < eq; i++) {
			k4[i] = pV[i]*dt;
		}
		// Calculate new value
		for (int i = 0; i < eq; i++) {
			pUOut[i] = pU[i] + (1/6.0)*
				(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
		}
		
	}

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

