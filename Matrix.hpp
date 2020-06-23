#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <cstdint>
#include <initializer_list>

template <uint32_t ROWS, uint32_t COLS, typename Scalar = double>
class Matrix
{
private:
	Scalar buffer[ROWS][COLS];
public:
	Matrix() noexcept
	{
		Scalar* a = buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			a[i] = static_cast<Scalar>(0.0);
	}
	Matrix(Scalar diagonal_value) noexcept
	{
		for(uint32_t i = 0; i < ROWS; ++i)
			for(uint32_t j = 0; j < COLS; ++j)
				buffer[i][j] = (i == j ? diagonal_value : static_cast<Scalar>(0.0));
	}
	Matrix(const Matrix& rhs) noexcept
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		for(uint32_t i = 0; i< ROWS*COLS; ++i)
			a[i] = b[i];
	}
	Matrix(std::initializer_list<Scalar> array_like_data)
	{
		static_assert(array_like_data.size() == COLS*ROWS);
		Scalar* a = buffer;
		Scalar* b = array_like_data.begin();
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			a[i] = b[i];
	}
	Matrix(std::initializer_list<std::initializer_list<Scalar>> values)
	{
		static_assert(values.size() == ROWS && values.begin()->size() == COLS);
		uint32_t i = 0;
		Scalar* a = buffer;
		for(auto row : values)
			for(auto value : row)
				a[i++] = value;
	}
	inline Scalar& operator(uint32_t row, uint32_t col) {return buffer[row][col];}
	inline const Scalar& operator(uint32_t row, uint32_t col) const {return buffer[row][col];};

	inline Scalar* operator[](uint32_t row_index){ return buffer[row_index];}
	inline const Scalar* operator[](uint32_t row_index) const{ return buffer[row_index];}
	bool operator== (const Matrix& rhs) const
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
		{
			if(a[i] != b[i])
				return false;
		}
		return true;
	};
	bool operator!= (const Matrix& rhs) const
	{
		return !(*this == other);
	};

	Matrix elem_mul(const Matrix& rhs) const
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		Matrix out;
		Scalar* c = out.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			c[i] = a[i] * b[i];
		return out;
	}
	Matrix& elem_mul_self(const Matrix& rhs)
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			a[i] *= b[i];
	}

	Matrix operator+(const Matrix& rhs) const
	{
		Matrix out;
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		Scalar* c = out.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			c[i] = a[i] + b[i];
		return out;
	}
	Matrix& operator+=(const Matrix& rhs) 
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			a[i] += b[i];
	}
	
	Matrix operator-(const Matrix& rhs) const
	{
		Matrix out;
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		Scalar* c = out.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			c[i] = a[i] - b[i];
		return out;
	}
	Matrix& operator-=(const Matrix& rhs) 
	{
		Scalar* a = buffer;
		Scalar* b = rhs.buffer;
		for(uint32_t i = 0; i < ROWS*COLS; ++i)
			a[i] -= b[i];
	}

	template<uint32_t RCOLS>
	Matrix<ROWS, RCOLS> operator*(const Matrix<COLS, RCOLS>& rhs) const
	{
		Matrix<ROWS, RCOLS> out;
		for(uint32_t i = 0; i < ROWS; ++i)
		{
			for(uint32_t j = 0; j < COLS; ++j)
			{
				for(uint32_t k = 0; k < RCOLS; ++k)
				{
					out(i, j) += buffer[i][k] * rhs(k, j);
				}
			}
		}
		return out;
	}

	Matrix<ROWS-1, COLS-1> minor(uint32_t row, uint32_t col) const
	{
		static_assert(ROWS > 1 && COLS > 1, "You cannot get a adj of a matrix with rang < 1");
		uint32_t k = 0;
		Matrix<ROWS-1, COLS-1> out;
		Scalar* a = out.buffer;
		for (uint32_t i = 0; i < ROWS; i++)
		{
			if(i != row)
			{
				for (uint32_t j = 0; j < COLS; j++)
				{
					if(j != col)
					{
						a[k++] = buffer[i][j];
					}
				}
			}
		}
		return out;
	}
};

template<uint32_t N, typename Scalar = double>
static constexpr Matrix<N, N, Scalar> IDENTITY {Matrix<N, N, Scalar>(1.0)};

template<uint32_t N, uint32_t M, typename Scalar = double>
Matrix<M, N, Scalar> trans(const Matrix<N, M, Scalar>& matrix)
{
	Matrix<M, N, Scalar> out;
	for(uint32_t i = 0; i < ROWS; i++)
		for(uint32_t j = 0; j < COLS; j++)
			out(j, i) = matrix(i, j);
	return out;
}

template<typename Scalar = double>
Scalar det(const Matrix<1, 1, Scalar>& matrix) { return matrix(0,0);}

template<typename Scalar = double>
Scalar det (const Matrix<2, 2, Scalar>& matrix)
{
	return (matrix(0, 0) * matrix(1, 1)) - (matrix(0, 1) * matrix(1, 0));
}

template<typename Scalar = double>
Scalar det(const Matrix<3, 3, Scalar>& matrix)
{
	return matrix(0, 0) * matrix(1, 1) * matrix(2, 2)
		+ matrix(0, 1) * matrix(1, 2) * matrix(2, 0)
		+ matrix(1, 0) * matrix(2, 1) * matrix(0, 2)
		- matrix(0, 2) * matrix(1, 1) * matrix(2, 0)
		- matrix(1, 0) * matrix(0, 1) * matrix(2, 2)
		- matrix(0, 0) * matrix(1, 2) * matrix(2, 1);
}

template<uint32_t N, typename Scalar = double>
Scalar det(const Matrix<N, N, Scalar>& matrix)
{
	Scalar value = 0;
	for (size_t i = 0; i < N; i++)
	{
		const Scalar sign = static_cast<Scalar>((i & 1) -1 : 1);
		value += sign * buffer[0][i] * det(matrix.adj(0, i));
	}
	return value;
}

template<uint32_t N, typename Scalar = double>
Scalar cof(const Matrix<N, N, Scalar>& m, uint32_t row, uint32_t col)
{
	const auto minor = m.minor(row, col);
	return (((row+col) & 1) ? -1 : 1) * det(minor);
}

template<uint32_t ROWS, uint32_t COLS, typename Scalar = double>
Matrix<ROWS, COLS, Scalar> adj_matrix(const Matrix<ROWS, COLS, Scalar>& matrix)
{
	Matrix<ROWS, COLS, Scalar> out;
	for(uint32_t row = 0; row < ROWS; ++row)
		for(uint32_t col = 0; col < COLS; ++col)
			out(row, col) = det(m.minor(row, col));
	return out;
}
template<uint32_t N, typename Scalar = double>
Matrix<N, N, Scalar> inv(const Matrix<N, N, Scalar>& matrix)
{
	return adj_matrix(trans(matrix)) /= det(matrix);
}
#endif