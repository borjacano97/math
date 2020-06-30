#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <cstdint>
#include <initializer_list>
template <uint32_t ROWS, uint32_t COLS, typename Scalar = double>
class Matrix;

template<uint32_t N, typename Scalar = double>
using Column = Matrix<N, 1, Scalar>;
template<uint32_t N, typename Scalar = double>
using Row = Matrix<1, N, Scalar>;

template <uint32_t ROWS, uint32_t COLS, typename Scalar = double>
class Matrix
{
	class ColRef
	{
		friend class Matrix;
		Matrix& matrix;
		const int col_index;
		ColRef(const Matrix& matrix, int index) noexcept:
			matrix(matrix), col_index(index)
		{}
	public:
		ColRef(const ColRef&) = default;
		ColRef(ColRef&&) = default;

		Scalar& operator[](uint32_t index){ return m.at(index, col_index);}
		const Scalar& operator[](uint32_t index) const { return m.at(index, col_index);}
		
		constexpr uint32_t size() const { return ROWS;}
		uint32_t index() const {return col_index;}

		ColRef& operator+=(Scalar n) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) += n;
			return *this;
		}
		ColRef& operator-=(Scalar n) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) -= n;
			return *this;
		}
		ColRef& operator*=(Scalar n) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) *= n;
			return *this;
		}
		ColRef& operator/=(Scalar n) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) /= n;
			return *this;
		}
		friend ColRef& operator+=(const ColRef& col, Scalar n){return col += scalar;}
		friend ColRef& operator-=(const ColRef& col, Scalar n){return col -= scalar;}
		friend ColRef& operator*=(const ColRef& col, Scalar n){return col *= scalar;}
		friend ColRef& operator/=(const ColRef& col, Scalar n){return col /= scalar;}

		ColRef& operator+=(const Matrix<ROWS, 1, Scalar>& col_matrix) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) += col_matrix(row, 0);
			return *this;
		}
		ColRef& operator-=(const Matrix<ROWS, 1, Scalar>& col_matrix) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) -= col_matrix(row, 0);
			return *this;
		}
		ColRef& operator*=(const Matrix<ROWS, 1, Scalar>& col_matrix) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) *= col_matrix(row, 0);
			return *this;
		}
		ColRef& operator/=(const Matrix<ROWS, 1, Scalar>& col_matrix) const
		{
			for(uint32_t row = 0; row < ROWS; row++)
				m(row, col_index) /= col_matrix(row, 0);
			return *this;
		}
		friend Matrix<ROWS, 1, Scalar>& operator+=(Matrix<ROWS, 1, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t row = 0; row < ROWS; row++)
				col_matrix(row, col_index) += col_ref[row];
			return col_matrix;
		}
		friend Matrix<ROWS, 1, Scalar>& operator-=(Matrix<ROWS, 1, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t row = 0; row < ROWS; row++)
				col_matrix(row, col_index) -= col_ref[row];
			return col_matrix;
		}
		friend Matrix<ROWS, 1, Scalar>& operator*=(Matrix<ROWS, 1, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t row = 0; row < ROWS; row++)
				col_matrix(row, col_index) *= col_ref[row];
			return col_matrix;
		}
		friend Matrix<ROWS, 1, Scalar>& operator/=(Matrix<ROWS, 1, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t row = 0; row < ROWS; row++)
				col_matrix(row, col_index) /= col_ref[row];
			return col_matrix;
		}
		friend Matrix<ROWS, 1, Scalar> operator+(Matrix<ROWS, 1, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix += col_ref;
		}
		friend Matrix<ROWS, 1, Scalar> operator-(Matrix<ROWS, 1, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix -= col_ref;
		}
		friend Matrix<ROWS, 1, Scalar> operator*(Matrix<ROWS, 1, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix *= col_ref;
		}
		friend Matrix<ROWS, 1, Scalar> operator/(Matrix<ROWS, 1, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix /= col_ref;
		}
	};
	class RowRef
	{
		friend class Matrix;
		Matrix& matrix;
		const int row_index;
		RowRef(const Matrix& matrix, int index) noexcept:
			matrix(matrix), row_index(index)
		{}
	public:
		RowRef(const RowRef&) = default;
		RowRef(RowRef&&) = default;

		Scalar& operator[](uint32_t index){ return m.at(row_index, index);}
		const Scalar& operator[](uint32_t index) const { return m.at(row_index, index);}
		
		constexpr uint32_t size() const { return COLS;}
		uint32_t index() const {return row_index;}
		
		RowRef& operator+=(Scalar n) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) += n;
			return *this;
		}
		RowRef& operator-=(Scalar n) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) -= n;
			return *this;
		}
		RowRef& operator*=(Scalar n) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) *= n;
			return *this;
		}
		RowRef& operator/=(Scalar n) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) /= n;
			return *this;
		}

		friend RowRef& operator+=(const RowRef& row, Scalar n){return row += scalar;}
		friend RowRef& operator-=(const RowRef& row, Scalar n){return row -= scalar;}
		friend RowRef& operator*=(const RowRef& row, Scalar n){return row *= scalar;}
		friend RowRef& operator/=(const RowRef& row, Scalar n){return row /= scalar;}

		RowRef& operator+=(const Matrix<1, ROWS, Scalar>& row_matrix) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) += row_matrix(0, col);
			return *this;
		}
		RowRef& operator-=(const Matrix<1, ROWS, Scalar>& row_matrix) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) -= row_matrix(0, col);
			return *this;
		}
		RowRef& operator*=(const Matrix<1, ROWS, Scalar>& row_matrix) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) *= row_matrix(0, col);
			return *this;
		}
		RowRef& operator/=(const Matrix<1, ROWS, Scalar>& row_matrix) const
		{
			for(uint32_t col = 0; col < COLS; col++)
				m(row_index, col) /= row_matrix(0, col);
			return *this;
		}
		friend Matrix<1, COLS, Scalar>& operator+=(Matrix<1, COLS, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t col = 0; col < COLS; col++)
				row_matrix(row_index, col) += col_ref[col];
			return col_matrix;
		}
		friend Matrix<1, COLS, Scalar>& operator-=(Matrix<1, COLS, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t col = 0; col < COLS; col++)
				row_matrix(row_index, col) -= col_ref[col];
			return col_matrix;
		}
		friend Matrix<1, COLS, Scalar>& operator*=(Matrix<1, COLS, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t col = 0; col < COLS; col++)
				row_matrix(row_index, col) *= col_ref[col];
			return col_matrix;
		}
		friend Matrix<1, COLS, Scalar>& operator/=(Matrix<1, COLS, Scalar>& col_matrix, const ColRef& col_ref)
		{
			for(uint32_t col = 0; col < COLS; col++)
				row_matrix(row_index, col) /= col_ref[col];
			return col_matrix;
		}
		friend Matrix<1, COLS, Scalar> operator+(Matrix<1, COLS, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix += col_ref;
		}
		friend Matrix<1, COLS, Scalar> operator-(Matrix<1, COLS, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix -= col_ref;
		}
		friend Matrix<1, COLS, Scalar> operator*(Matrix<1, COLS, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix *= col_ref;
		}
		friend Matrix<1, COLS, Scalar> operator/(Matrix<1, COLS, Scalar> col_matrix, const ColRef& col_ref)
		{
			return col_matrix /= col_ref;
		}
	};
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

	inline ColRef getColRef(uint32_t col_index){ return ColRef(*this, col_index);}
	inline RowRef getRowRef(uint32_t row_index){ return RowRef(*this, row_index);}

	Matrix<ROWS, 1, Scalar> col(uint32_t index) const
	{
		Column out;
		for(uint32_t row = 0; row < ROWS; row++)
		{
			out(row, index) = at(row, index);
		}
		return out;
	}
	Matrix<1, COLS, Scalar> row(uint32_t index) const
	{
		Row out;
		for(uint32_t col = 0; col < COLS; col++)
		{
			out(index, col) = at(index, col);
		}
		return out;
	}
	void swapCols(uint32_t i, uint32_t j)
	{
		Scalar aux;
		for(uint32_t row = 0; row < ROWS; row++)
		{
			aux = at(row, i);
			at(row, i) = at(row, j);
			at(row, j) = aux;
		}
	}
	
	void swapRows(uint32_t i, uint32_t j)
	{
		Scalar aux;
		for(uint32_t col = 0; col < COLS; col++)
		{
			aux = at(i, col);
			at(i, col) = at(j, col);
			at(j, col) = aux;
		}
	}

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
bool may_have_inv(const Matrix<N, N, Scalar>& matrix)
{
	return det(matrix);
}
template<uint32_t N, typename Scalar = double>
Matrix<N, N, Scalar> inv(const Matrix<N, N, Scalar>& matrix)
{
	return adj_matrix(trans(matrix)) /= det(matrix);
}

template<uint32_t ROWS,uint32_t COLS, typename Scalar = double>
Matrix<ROWS, COLS, Scalar> triangulate(Matrix<ROWS, COLS, Scalar> matrix)
{
	for (uint32_t i = 0; i < ROWS-1; i++)
	{
		//Make the pivot 1, and transform the whole row as a result
		const Scalar raw_pivot = matrix(i, i);
		for(uint32_t k = 0; k < COLS; k++)
			matrix(i, k) /= raw_pivot;
		//for the remaining rows
		for (uint32_t j = i+1; j < ROWS; j++)
		{
			//The number to become 0
			const Scalar n = matrix(j, i);
			// Skip this row if already 0
			if(n == Scalar(0)) continue;

			// Fj - nFi
			for(uint32_t k = 0; k < COLS; k++)
				matrix(j, k) -= n * matrix(i, k);
		}
	}
	return matrix;
}
#endif