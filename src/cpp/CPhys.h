#include "Matrix.h"
namespace CPhys{
	namespace EigVal{
		void 	jacobiMethod(Matrix& A, Matrix& R, int N);
	}
	/*
	namespace Euler_ODE{
		void	step(Matrix& u, Matrix& v, double* step);
	}
	namespace RK4_ODE{
		void	step(Matrix& u, Matrix& v, double* step);
	}
	*/
}
// Make a unnamed namespace to hide functions from the user
// prefix m (abbr. for 'my') is needed
namespace{
	namespace mEigVal{
		double 	maxoffdiag        (double** A, int* k, int* l, int N);
		void	rotate(double** A, double** R, int  k, int  l, int N);
	}
}

