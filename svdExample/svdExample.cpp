#include "svdExample.h"



void main()
{
	std::cout << "===Singular Value Decomposition Example===" << std::endl;
	Matrix A(2, 2);
	A.At(0, 0) = 3.0f;
	A.At(1, 0) = 3.0f;
	A.At(0, 1) = 3.0f;
	A.At(1, 1) = 3.0f;
	Matrix B(2, 2);
	B.At(0, 0) = 3.0f;
	B.At(1, 0) = 3.0f;
	B.At(0, 1) = 3.0f;
	B.At(1, 1) = 3.0f;
	std::cout << "Test Matrix A = " << std::endl << A.ToString();
	std::cout << "Test Matrix B = " << std::endl << B.ToString();
	std::cout << "Test Matrix A*B = " << std::endl << Matrix::Multiply(A,B).ToString();
	std::cin.get();
}