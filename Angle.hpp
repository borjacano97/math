#ifndef ANGLE_HPP
#define ANGLE_HPP

#define USE_STANDARD_FUNCTIONS
#ifdef USE_STANDARD_FUNCTIONS
#include <cmath>

#define COS(x) std::cos(x)
#define SIN(x) std::sin(x)
#define TAN(x) std::tan(x)
#define ACOS(x) std::acos(x)
#define ASIN(x) std::asin(x)

#else // user provided trigonometry functions
#define COS(x) 
#define SIN(x) 
#define TAN(x) 
#define ACOS(x) 
#define ASIN(x) 
#endif


template<typename T>
constexpr auto PI = static_cast<T>(3.141592653589793238462643383279502884);
template<typename T>
constexpr auto TWO_PI = 2.0*PI<T>;
template<typename T>
constexpr auto HALF_PI = PI<T> / static_cast<T>(2);

template<typename Real>
constexpr static auto rad2deg = static_cast<Real>(180.0) / PI<Real>;
template<typename Real>
constexpr static auto deg2rad = PI / static_cast<Real>(180.0);

template<typename Real>
class Radians;
template<typename Real>
class Degrees
{
	Real value;
public:
	constexpr Degrees() noexcept: value(0){}
	constexpr Degrees(const Degrees&) = default;
	constexpr Degrees(Degrees&&) = default;
	explicit constexpr Degrees(Real value) noexcept: value(value){}
	explicit constexpr Degrees(const Radians& radians) noexcept:
		value(radians*rad2deg<Real>){}
	constexpr operator Radians() noexcept const { return Radians(value*rad2deg);}

	constexpr Real& raw() { return value;}
	constexpr const Real& raw() const { return value;}

	constexpr bool operator==(const Degrees& rhs) const {return value == rhs.value;}
	constexpr bool operator!=(const Degrees& rhs) const {return value != rhs.value;}

	constexpr bool operator<(const Degrees& rhs) const {return value < rhs.value;}
	constexpr bool operator>(const Degrees& rhs) const {return value > rhs.value;}

	constexpr bool operator<=(const Degrees& rhs) const {return value <= rhs.value;}
	constexpr bool operator>=(const Degrees& rhs) const {return value >= rhs.value;}

};
template<typename Real>
struct Radians
{
	Real value;
public:
	constexpr Radians() noexcept: value(0){}
	constexpr Radians(const Radians&) noexcept = default;
	constexpr Radians(Radians&&) noexcept = default;
	explicit constexpr Radians(Real value) noexcept: value(value){}
	explicit constexpr Radians(const Degrees& degrees) noexcept:
		value(degrees*deg2rad<Real>){}
	constexpr operator Degrees() noexcept const { return Degrees(value*rad2deg<Real>);}

	constexpr Real& raw() { return value;}
	constexpr const Real& raw() const { return value;}

	constexpr bool operator==(const Radians& rhs) const {return value == rhs.value;}
	constexpr bool operator!=(const Radians& rhs) const {return value != rhs.value;}

	constexpr bool operator<(const Radians& rhs) const {return value < rhs.value;}
	constexpr bool operator>(const Radians& rhs) const {return value > rhs.value;}

	constexpr bool operator<=(const Radians& rhs) const {return value <= rhs.value;}
	constexpr bool operator>=(const Radians& rhs) const {return value >= rhs.value;}

};

template<typename Real>
class Angle
{
	Radians value;
public:
	constexpr Angle() noexcept: value(){}
	constexpr Angle(const Degrees& deg) noexcept: value(deg){}
	constexpr Angle(const Radians& rad) noexcept: value(rad){}

	constexpr bool operator==(const Angle& rhs) const {return value == rhs.value;}
	constexpr bool operator!=(const Angle& rhs) const {return value != rhs.value;}

	constexpr bool operator<(const Angle& rhs) const {return value < rhs.value;}
	constexpr bool operator>(const Angle& rhs) const {return value > rhs.value;}

	constexpr bool operator<=(const Angle& rhs) const {return value <= rhs.value;}
	constexpr bool operator>=(const Angle& rhs) const {return value >= rhs.value;}


	friend inline Angle cos(const Angle& x)	{ return Angle(Radians( COS(x.value)));}
	friend inline Angle acos(const Angle& x)	{ return Angle(Radians(ACOS(x.value)));}
	friend inline Angle sin(const Angle& x)	{ return Angle(Radians( SIN(x.value)));}
	friend inline Angle asin(const Angle& x)	{ return Angle(Radians(ASIN(x.value)));}
	friend inline Angle tan(const Angle& x)	{ return Angle(Radians( TAN(x.value)));}
	friend inline Angle atan(const Angle& x)	{ return Angle(Radians(ATAN(x.value)));}
};


#endif /* ANGLE_HPP */
