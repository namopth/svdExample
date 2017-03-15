#include "svdExample.h"

template <typename T> 
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

void jacobiSVD(const Matrix& input, Matrix& outU, Matrix& outD, Matrix& outV)
{
	static const float accErr = M_EPSILON;
	unsigned int matSize = input.GetColSize();
	outU = input;
	outV = outD = input.Identity();

	float converge = accErr + 1.0f;
	while (converge > accErr)
	{
		converge = 0.f;
		for (unsigned int j = 2; j < matSize; j++)
		{
			for (unsigned int i = 0; i < j; i++)
			{
				// compute submatrix of transpose U multiplied by U
				float alpha = 0.f, beta = 0.f, gamma = 0.f;
				for (unsigned int k = 0; k < matSize; k++)
				{
					alpha += outU.Get(i, k) * outU.Get(i, k);
					beta += outU.Get(j, k) * outU.Get(j, k);
					gamma += outU.Get(i, k) * outU.Get(j, k);
				}

				// compute convergence
				converge = max(converge, abs(gamma) / sqrt(alpha*beta));

				// compute jacobi rotation matrix
				float zeta = (beta - alpha) / (2.f*gamma);
				float t = sgn(zeta) / (abs(zeta) + sqrt(1 + zeta*zeta));
				float c = 1.f / sqrt(1.0f + t*t);
				float s = c * t;

				// update matrix U
				Matrix tempU = outU;
				for (unsigned int k = 0; k < matSize; k++)
				{
					outU.At(i, k) = c * tempU.Get(i, k) - s * tempU.Get(j, k);
					outU.At(j, k) = s * tempU.Get(i, k) + c * tempU.Get(j, k);
				}

				// update matrix V
				Matrix tempV = outV;
				for (unsigned int k = 0; k < matSize; k++)
				{
					outV.At(i, k) = c * tempV.Get(i, k) - s * tempV.Get(j, k);
					outV.At(j, k) = s * tempV.Get(i, k) + c * tempV.Get(j, k);
				}
			}
		}
	}

	for (unsigned int j = 0; j < matSize; j++)
	{
		float sqNorm = 0.f;
		for (unsigned int i = 0; i < matSize; i++)
		{
			sqNorm += outU.Get(j, i) * outU.Get(j, i);
		}
		outD.At(j,j) = sqrt(sqNorm);
		for (unsigned int i = 0; i < matSize; i++)
		{
			assert(outD.Get(j, j) > M_EPSILON);
			if (outD.Get(j, j) > M_EPSILON)
				outU.At(j, i) = outU.At(j, i) / outD.Get(j, j);
		}
	}

}

void main()
{
	std::cout << "===Singular Value Decomposition Example===" << std::endl;
	Matrix A(3, 3);
	A.At(0, 0) = 1.0f;
	A.At(1, 0) = 3.0f;
	A.At(2, 0) = 2.0f;
	A.At(0, 1) = 5.0f;
	A.At(1, 1) = 6.0f;
	A.At(2, 1) = 4.0f;
	A.At(0, 2) = 7.0f;
	A.At(1, 2) = 8.0f;
	A.At(2, 2) = 9.0f;
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
	std::cout << "Reconstructed Matrix A = " << std::endl << AU.Multiply(AD).Multiply(AV.Transpose()).ToString();

	std::cin.get();
}