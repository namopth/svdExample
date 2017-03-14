#pragma once

#include <math.h>
#include <algorithm>
#include <string>
#include <assert.h>

// SIMPLE and SUPER SLOW NM-dimension Matrix Processing Class
// Used to demonstrating a singular value decomposition algorithm

#define MATRIX_COLUMN_MAJOR // Matrix data storing layout definition

class Matrix
{
public:
	Matrix(const unsigned int col, const unsigned int row)
		: m_pData(nullptr)
		, m_uiCol(col)
		, m_uiRow(row)
	{
		m_pData = new float[m_uiCol*m_uiRow];
		for (unsigned int i = 0; i < m_uiCol * m_uiRow; i++)
			m_pData[i] = 0.f;
	}

	Matrix(const Matrix& rhs)
		: m_pData(nullptr)
		, m_uiCol(rhs.m_uiCol)
		, m_uiRow(rhs.m_uiRow)
	{
		m_pData = new float[m_uiCol*m_uiRow];

		for (unsigned int i = 0; i < m_uiCol * m_uiRow; i++)
			m_pData[i] = rhs.m_pData[i];
	}

	~Matrix()
	{
		if(m_pData)
			delete[] m_pData;
		m_pData = nullptr;
	}

	const float Get(const unsigned int col, const unsigned int row) const
	{
		assert(row < m_uiRow && col < m_uiCol);
#ifdef MATRIX_COLUMN_MAJOR
		return m_pData[row * m_uiCol + col];
#else
		return m_pData[col * m_uiRow + row];
#endif
	}

	float& At(const unsigned int col, const unsigned int row)
	{
		assert(row < m_uiRow && col < m_uiCol);
#ifdef MATRIX_COLUMN_MAJOR
		return m_pData[row * m_uiCol + col];
#else
		return m_pData[col * m_uiRow + row];
#endif
	}

	const std::string ToString() const
	{
		std::string result = "";
		for (unsigned int r = 0; r < m_uiRow; r++)
		{
			result += "[\t";
			for (unsigned int c = 0; c < m_uiCol; c++)
			{
				result += std::to_string(Get(c, r));
				result += "\t";
			}
			result += "]\n";
		}
		return result;
	}

	const unsigned int GetColSize() const { return m_uiCol; }
	const unsigned int GetRowSize() const { return m_uiRow; }

	static Matrix Identity(const unsigned int col, const unsigned int row)
	{
		Matrix result(col, row);
		for (unsigned int i = 0; i < (std::min)(col, row); i++)
			result.At(i, i) = 1.f;
		return result;
	}

	Matrix Identity()
	{
		return Identity(m_uiCol, m_uiRow);
	}

	static Matrix Multiply(const Matrix& lhs, const Matrix& rhs)
	{
		assert(lhs.GetColSize() == rhs.GetRowSize());
		Matrix result(lhs.GetRowSize(), rhs.GetColSize());
		for (unsigned int i = 0; i < result.GetColSize(); i++)
		{
			for (unsigned int j = 0; j < result.GetRowSize(); j++)
			{
				float sum = 0.f;
				for (unsigned int k = 0; k < lhs.GetColSize(); k++)
				{
					sum += lhs.Get(k, j) * rhs.Get(i, k);
				}
				result.At(i, j) = sum;
			}
		}
		return result;
	}

	Matrix Multiply(const Matrix& rhs)
	{
		return Multiply(*this, rhs);
	}

	Matrix operator*(const Matrix& rhs)
	{
		return Multiply(rhs);
	}

	static Matrix Transpose(const Matrix& rhs)
	{
		Matrix result(rhs);
		for (unsigned int i = 0; i < result.GetColSize(); i++)
		{
			for (unsigned int j = 0; j < result.GetRowSize(); j++)
			{
				result.At(i, j) = rhs.Get(j, i);
			}
		}
		return result;
	}

	Matrix Transpose()
	{
		return Transpose(*this);
	}

private:
	float* m_pData;
	unsigned int m_uiCol;
	unsigned int m_uiRow;
};