#include "Matrix.h"
#include "ODE_Interface.cpp"
namespace CPhys{
	namespace EigVal{
		void 	jacobiMethod(Matrix& A, Matrix& R, int N);
	}

	namespace MatOp{
		void	sortCol(Matrix& A, Vector& v);
	}
namespace VecOp{
		double	normalize   (Vector& v, double dx);
	}

	namespace ODE{
	      	// This function includes function pointer paramters
		void	solveEuler(Matrix&  u, Vector&        uInit,   
				   double  dt, ODE_Interface*   ode);

		void	solveRK4  (Matrix&  u, Vector& 	    uInit, 
				   double  dt, ODE_Interface*   ode);
	}

}
// Make a unnamed namespace to hide functions from the user
// prefix m (abbr. for 'my') is used to avoid name collision
namespace{
	namespace mMatOp{
		int 	compareTwoRows(const void* rowA, const void* rowB);
	}
	namespace mEigVal{
		double 	maxoffdiag        (double** A, int* k, int* l, int N);
		void	rotate(double** A, double** R, int  k, int  l, int N);
	}
	namespace mODE{
	      	// This function includes function pointer parameters
		void	solve(Matrix&  u, Vector&        uInit,  
			      double& dt, ODE_Interface*   ode,
		     	void (*stepFunc) (double*,double,double*       ,
			               	  double&,int&  ,ODE_Interface*));

		void 	stepEuler(double* pU, double   t, double*        pUOut, 
			          double& dt, int&    eq, ODE_Interface*   ode);

		void 	stepRK4(  double* pU, double  t, double* pUOut, 
			          double& dt, int&    eq, ODE_Interface*   ode);
	}
}

