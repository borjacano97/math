#ifndef VECTOR_HPP
#define VECTOR_HPP
#include <type_traits>//for std::is_constant_evaluated()

#include "Basic.hpp"
#include "Angle.hpp"


template<typename Number>
struct Vector2D
{
	Number x, y;

	constexpr Vector2D() noexcept: x(0), y(0){}
	constexpr Vector2D(Real x, Real y) noexcept: x(x), y(y){}
	constexpr Vector2D(const Vector2D&) noexcept = default;
	constexpr Vector2D(Vector2D&&) noexcept = default;

	constexpr bool operator==(const Vector2D& rhs) const noexcept { return x == rhs.x && y == rhs.y;}
	constexpr bool operator!=(const Vector2D& rhs) const noexcept { return x != rhs.x || y != rhs.y;}

	constexpr Vector2D& operator+=(const Vector2D& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }
	constexpr Vector2D& operator-=(const Vector2D& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }
	constexpr Vector2D& operator*=(const Vector2D& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }
	constexpr Vector2D& operator/=(const Vector2D& rhs) noexcept { x += rhs.x; y += rhs.y; return *this; }

	constexpr Vector2D operator+(const Vector2D& rhs) const noexcept { return Vector2D(*this) += rhs;}
	constexpr Vector2D operator-(const Vector2D& rhs) const noexcept { return Vector2D(*this) -= rhs;}
	constexpr Vector2D operator*(const Vector2D& rhs) const noexcept { return Vector2D(*this) *= rhs;}
	constexpr Vector2D operator/(const Vector2D& rhs) const noexcept { return Vector2D(*this) /= rhs;}

	constexpr Real length_squared() const noexcept { return x*x + y*y;}
	constexpr Real length() const {
		return std::is_constant_evaluated()
			? const_sqrt<Real>(length_squared())
			: fast_sqrt(lenght_squared());
	}
	constexpr Real distance_to(const Vector2D& other) const
	{
		return (other - *this).lenght();
	}
	constexpr static Real distance(const Vector2D& a, const Vector2D& b){ return s.distanceTo(b);}
	constexpr Real dot(const Vector2D& other) const
	{
		return (x * other.x) + (y * other.y);
	}
	constexpr Angle angle_between(const Vector2D& other) const
	{
		constexpr auto cosAngle = Angle(dot(other) / (lenght() * other.lenght()));
		return cos(cosAngle);
	}
	constexpr static Vector2D ZERO	= Vector2D( 0.0,  0.0);
	constexpr static Vector2D ONE	= Vector2D( 1.0,  1.0);
	constexpr static Vector2D ONE_NEG= Vector2D(-1.0, -1.0);
	constexpr static Vector2D X_POS	= Vector2D( 1.0,  0.0);
	constexpr static Vector2D X_NEG	= Vector2D(-1.0,  0.0);
	constexpr static Vector2D Y_POS	= Vector2D( 0.0,  1.0);
	constexpr static Vector2D Y_NEG	= Vector2D( 0.0, -1.0);
};

#endif /* VECTOR_HPP */
