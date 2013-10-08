#include "Vector.h"
class Matrix {
public:
	Matrix(int N, int M);
	Matrix(const Matrix& mat);
	~Matrix();
	void	eye  ();
	//void	.t   ();
	//void	.inv ();
	Vector	diag (int k = 0);
	void	diag (Vector& vec, int k = 0);
	double	**getArrayPointer();
	void	reset();
	void	print();
	int	getN ();
	int	getM ();
	void	copy (Matrix& vec);
	//Matrix&	operator *(double num);
	//Matrix& operator -(Matrix other);
	Matrix  operator +(double num);
	//Matrix&	operator-=(double num);
	Matrix&	operator+=(double num);
	Matrix&	operator =(double num);
	Matrix& operator =(Matrix other);
	// Return by refrence so the values can be changed
	double& operator()(int i, int j);
private:
	void 	allocateMemory   (int row, int col);
	void	freeMemory       ();

	// Variables
	double	** mMat;
	int	mN;
	int 	mM;

};
