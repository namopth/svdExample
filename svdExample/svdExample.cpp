#include "svdExample.h"

template <typename T> 
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

void jacobiSVD(const Matrix& input, Matrix& outU, Matrix& outD, Matrix& outV)
{
	static const float accErr = M_EPSILON;
	float converge = accErr + 1.0f;
	while (converge > accErr)
	{
		for (unsigned int j = 2; j < input.GetColSize(); j++)
		{
			for (unsigned int i = 0; i < j; i++)
			{
				// compute submatrix of transpose U multiplied by U

				// compute convergence

				// compute jacobi rotation matrix

				// update matrix U

				// update matrix V
			}
		}
	}
}

void main()
{
	std::cout << "===Singular Value Decomposition Example===" << std::endl;
	Matrix A(2, 2);
	A.At(0, 0) = 1.0f;
	A.At(1, 0) = 2.0f;
	A.At(0, 1) = 3.0f;
	A.At(1, 1) = 4.0f;
	//Matrix B(2, 2);
	//B.At(0, 0) = 2.0f;
	//B.At(1, 0) = 0.0f;
	//B.At(0, 1) = 1.0f;
	//B.At(1, 1) = 2.0f;
	std::cout << "Test Matrix A = " << std::endl << A.ToString();
	//std::cout << "Test Transpose of Matrix A = " << std::endl << A.Transpose().ToString();
	//std::cout << "Test Matrix B = " << std::endl << B.ToString();
	//std::cout << "Test Matrix A*B = " << std::endl << Matrix::Multiply(A,B).ToString();


	std::cout << "===Jacobi SVD===" << std::endl;
	Matrix AU(A.Identity());
	Matrix AD(A.Identity());
	Matrix AV(A.Identity());
	jacobiSVD(A, AU, AD, AV);
	std::cout << "U of Matrix A = " << std::endl << AU.ToString();
	std::cout << "D of Matrix A = " << std::endl << AD.ToString();
	std::cout << "V of Matrix A = " << std::endl << AV.ToString();
	std::cout << "Reconstruced Matrix A = " << std::endl << AU.Multiply(AD).Multiply(AV).ToString();

	std::cin.get();
}