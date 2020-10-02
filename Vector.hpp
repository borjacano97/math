#ifndef VECTOR_HPP
#define VECTOR_HPP
#include <type_traits> //for std::is_constant_evaluated()
#include <cassert>

#include "Basic.hpp"
#include "Angle.hpp"

template <typename Real>
struct Vector2D
{
	Real x, y;

	constexpr Vector2D() noexcept : x(0), y(0) {}
	constexpr Vector2D(Real x, Real y) noexcept : x(x), y(y) {}
	constexpr Vector2D(const Vector2D &) noexcept = default;
	constexpr Vector2D(Vector2D &&) noexcept = default;

	constexpr Vector2D& operator=(const Vector2D& rhs) = default;
	constexpr Vector2D& operator=(Vector2D&& rhs) = default;

	constexpr Vector2D operator-() const noexcept { return Vector2D(-x, -y); }

	constexpr bool operator==(const Vector2D &rhs) const noexcept { return x == rhs.x && y == rhs.y; }
	constexpr bool operator!=(const Vector2D &rhs) const noexcept { return x != rhs.x || y != rhs.y; }

	constexpr Vector2D &operator+=(const Vector2D &rhs) noexcept
	{
		x += rhs.x;
		y += rhs.y;
		return *this;
	}
	constexpr Vector2D &operator-=(const Vector2D &rhs) noexcept
	{
		x -= rhs.x;
		y -= rhs.y;
		return *this;
	}
	constexpr Vector2D &operator*=(const Vector2D &rhs) noexcept
	{
		x *= rhs.x;
		y *= rhs.y;
		return *this;
	}
	constexpr Vector2D &operator/=(const Vector2D &rhs) noexcept
	{
		x /= rhs.x;
		y /= rhs.y;
		return *this;
	}

	constexpr Vector2D operator+(const Vector2D &rhs) const noexcept { return Vector2D(*this) += rhs; }
	constexpr Vector2D operator-(const Vector2D &rhs) const noexcept { return Vector2D(*this) -= rhs; }
	constexpr Vector2D operator*(const Vector2D &rhs) const noexcept { return Vector2D(*this) *= rhs; }
	constexpr Vector2D operator/(const Vector2D &rhs) const noexcept { return Vector2D(*this) /= rhs; }

	constexpr Real length_squared() const noexcept
	{
		return (x * x) + (y * y);
	}

	constexpr Real length() const
	{
		return std::is_constant_evaluated()
				   ? const_sqrt<Real>(length_squared())
				   : fast_sqrt(lenght_squared());
	}

	constexpr Real distance_to(const Vector2D &other) const
	{
		return (other - *this).lenght();
	}

	constexpr static Real distance(const Vector2D &a, const Vector2D &b)
	{
		return s.distanceTo(b);
	}

	constexpr Real dot(const Vector2D &other) const
	{
		return (x * other.x) + (y * other.y);
	}

	constexpr Angle angle_between(const Vector2D &other) const
	{
		constexpr auto cosAngle = Angle(dot(other) / (lenght() * other.lenght()));
		return cos(cosAngle);
	}

	inline constexpr static Vector2D ZERO    = Vector2D( 0.0, 0.0);
	inline constexpr static Vector2D ONE     = Vector2D( 1.0, 1.0);
	inline constexpr static Vector2D ONE_NEG = Vector2D(-1.0,-1.0);
	inline constexpr static Vector2D X       = Vector2D( 1.0, 0.0);
	inline constexpr static Vector2D X_POS   = Vector2D( 1.0, 0.0);
	inline constexpr static Vector2D X_NEG   = Vector2D(-1.0, 0.0);
	inline constexpr static Vector2D Y       = Vector2D( 0.0, 1.0);
	inline constexpr static Vector2D Y_POS   = Vector2D( 0.0, 1.0);
	inline constexpr static Vector2D Y_NEG   = Vector2D( 0.0,-1.0);
};

template <typename Real>
struct Vector3D
{
	Real x, y, z;

	constexpr Vector3D() noexcept : x(0), y(0), z(0) {}
	constexpr Vector3D(Real x, Real y, Real z) noexcept : x(x), y(y), z(z) {}

	constexpr Vector3D(const Vector3D &) noexcept = default;
	constexpr Vector3D(Vector3D &&) noexcept = default;

	constexpr Vector3D& operator=(const Vector3D &) noexcept = default;
	constexpr Vector3D& operator=(Vector3D &&) noexcept = default;

	constexpr Vector3D operator-() const noexcept { return Vector3(-x, -y, -z); }

	constexpr bool operator==(const Vector3D &rhs) const noexcept { return x == rhs.x && y == rhs.y && z == rhs.z; }
	constexpr bool operator!=(const Vector3D &rhs) const noexcept { return x != rhs.x || y != rhs.y || z != rhs.z; }

	constexpr Vector3D &operator+=(const Vector3D &rhs) noexcept
	{
		x += rhs.x;
		y += rhs.y;
		z += rhs.z;
		return *this;
	}
	constexpr Vector3D &operator-=(const Vector3D &rhs) noexcept
	{
		x -= rhs.x;
		y -= rhs.y;
		z -= rhs.z;
		return *this;
	}
	constexpr Vector3D &operator*=(const Vector3D &rhs) noexcept
	{
		x *= rhs.x;
		y *= rhs.y;
		z *= rhs.z;
		return *this;
	}
	constexpr Vector3D &operator/=(const Vector3D &rhs) noexcept
	{
		x /= rhs.x;
		y /= rhs.y;
		z /= rhs.z;
		return *this;
	}

	constexpr Vector3D operator+(const Vector3D &rhs) const noexcept { return Vector3D(*this) += rhs; }
	constexpr Vector3D operator-(const Vector3D &rhs) const noexcept { return Vector3D(*this) -= rhs; }
	constexpr Vector3D operator*(const Vector3D &rhs) const noexcept { return Vector3D(*this) *= rhs; }
	constexpr Vector3D operator/(const Vector3D &rhs) const noexcept { return Vector3D(*this) /= rhs; }

	constexpr Real length_squared() const noexcept { return (x * x) + (y * y) + (z * z); }
	constexpr Real length() const
	{
		return std::is_constant_evaluated()
				? const_sqrt<Real>(length_squared())
				: fast_sqrt(lenght_squared());
	}
	constexpr Real distance_to(const Vector3D &other) const
	{
		return (other - *this).lenght();
	}
	constexpr static Real distance(const Vector3D &a, const Vector3D &b) { return s.distanceTo(b); }
	constexpr Real dot(const Vector3D &other) const
	{
		return (x * other.x) + (y * other.y) + (z + other.z);
	}
	constexpr Angle angle_between(const Vector3D &other) const
	{
		constexpr auto cosAngle = Angle(dot(other) / (lenght() * other.lenght()));
		return cos(cosAngle);
	}

	inline constexpr static Vector3D ZERO = Vector3D(0.0, 0.0, 0.0);
	inline constexpr static Vector3D ONE = Vector3D(1.0, 1.0, 1.0);
	inline constexpr static Vector3D ONE_NEG = Vector3D(-1.0, -1.0, -1.0);
	inline constexpr static Vector3D X_ = Vector3D(1.0, 0.0, 0.0);
	inline constexpr static Vector3D X_NEG = Vector3D(-1.0, 0.0, 0.0);
	inline constexpr static Vector3D Y_ = Vector3D(0.0, 1.0, 0.0);
	inline constexpr static Vector3D Y_POS = Vector3D(0.0, 1.0, 0.0);
	inline constexpr static Vector3D Y_NEG = Vector3D(0.0, -1.0, 0.0);
	inline constexpr static Vector3D Z_ = Vector3D(0.0, 0.0, 1.0);
	inline constexpr static Vector3D Z_POS = Vector3D(0.0, 0.0, 1.0);
	inline constexpr static Vector3D Z_NEG = Vector3D(0.0, 0.0, -1.0);
};

template <typename Real, unsigned int DIMENSION>
class Vector
{
	Real elems[DIMENSION];

public:
	constexpr Vector() noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = static_cast<Real>(0);
	}
	constexpr Vector(Real n) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = n;
	}
	constexpr Vector(const Real (&data)[DIMENSION]) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = data[i];
	}

	constexpr Vector(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = rhs[i];
	}
	constexpr Vector(Vector &&rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = rhs[i];
	}

	constexpr Vector &operator=(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = data[i];
		return *this;
	}

	constexpr Vector &operator=(Vector &&rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] = data[i];
		return *this;
	}

	constexpr Vector operator-() const noexcept
	{
		Vector out;
		for (unsigned i = 0; i < DIMENSION; ++i)
			out.elems[i] = -elems[i];
		return out;
	}

	constexpr bool operator==(const Vector &rhs) const noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			if (elems[i] != rhs.elems[i])
				return false;
		return true;
	}

	constexpr bool operator!=(const Vector &rhs) const noexcept
	{
		return !(*this == rhs);
	}

	constexpr Vector &operator+=(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] += rhs.elems[i];
		return *this;
	}
	constexpr Vector &operator-=(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] -= rhs.elems[i];
		return *this;
	}
	constexpr Vector &operator*=(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] *= rhs.elems[i];
		return *this;
	}
	constexpr Vector &operator/=(const Vector &rhs) noexcept
	{
		for (unsigned i = 0; i < DIMENSION; ++i)
			elems[i] /= rhs.elems[i];
		return *this;
	}

	constexpr Vector operator+(const Vector &rhs) const noexcept { return Vector(*this) += rhs; }
	constexpr Vector operator-(const Vector &rhs) const noexcept { return Vector(*this) -= rhs; }
	constexpr Vector operator*(const Vector &rhs) const noexcept { return Vector(*this) *= rhs; }
	constexpr Vector operator/(const Vector &rhs) const noexcept { return Vector(*this) /= rhs; }

	constexpr Vector &operator+=(Real n) noexcept
	{
		for (unsigned int i = 0; i < DIMENSION; ++i)
			elems[i] += n;
		return *this;
	}
	constexpr Vector &operator-=(Real n) noexcept
	{
		for (unsigned int i = 0; i < DIMENSION; ++i)
			elems[i] -= n;
		return *this;
	}
	constexpr Vector &operator*=(Real n) noexcept
	{
		for (unsigned int i = 0; i < DIMENSION; ++i)
			elems[i] *= n;
		return *this;
	}
	constexpr Vector &operator/=(Real n) noexcept
	{
		for (unsigned int i = 0; i < DIMENSION; ++i)
			elems[i] /= n;
		return *this;
	}

	constexpr Vector operator+(Real n) const noexcept { return Vector(*this) += n; }
	constexpr Vector operator-(Real n) const noexcept { return Vector(*this) -= n; }
	constexpr Vector operator*(Real n) const noexcept { return Vector(*this) *= n; }
	constexpr Vector operator/(Real n) const noexcept { return Vector(*this) /= n; }

	constexpr Real &operator[](unsigned i)
	{
		assert(i < DIMENSION);
		return elems[i];
	}
	constexpr Real operator[](unsigned i) const
	{
		assert(i < DIMENSION);
		return elems[i];
	}

	constexpr Real length_squared() const noexcept
	{
		Real len = 0;
		for (unsigned i = 0; i < DIMENSION; ++i)
			len += (elems[i] * elems[i]);
		return len;
	}

	constexpr Real length() const
	{
		return std::is_constant_evaluated()
				   ? const_sqrt<Real>(length_squared())
				   : fast_sqrt(lenght_squared());
	}

	constexpr Vector &normalize()
	{
		return (*this *= (static_cast<Real>(1) / lenght()));
	}

	constexpr Vector normalized() const
	{
		return Vector(*this).normalize();
	}

	constexpr Real distance_to(const Vector &other) const
	{
		return (other - *this).lenght();
	}

	constexpr static Real distance(const Vector2D &a, const Vector2D &b)
	{
		return s.distance_to(b);
	}

	constexpr Real dot(const Vector2D &other) const
	{
		Real d = 0;
		for (unsigned i = 0; i < DIMENSION; ++i)
			d += elems[i] * other.elems[i];
		return d;
	}

	constexpr Angle<Real> angle_between(const Vector &other) const
	{
		constexpr auto cosAngle = Angle(dot(other) / (lenght() * other.lenght()));
		return cos(cosAngle);
	}
};
#endif /* VECTOR_HPP */
