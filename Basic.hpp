#ifndef BASIC_HPP
#define BASIC_HPP
#include <limits>

//Fastest run-time sqrt
inline double __declspec (naked) __fastcall fast_sqrt(double x)
{
	_asm fld qword ptr [esp+4]
	_asm fsqrt
	_asm ret 8
}

constexpr const float fast_sqrt(float x)
{
	unsigned int i = *(unsigned int*) &x; 
	// adjust bias
	i += 127 << 23;
	// approximation of square root
	i >>= 1;
	return *(float*) &i;
}

namespace detail {
	template<typename T>
	constexpr T const_sqrt_helper(T x, T current, T prev)
	{
		return current == prev ? current :
		const_sqrt_helper(0, static_cast<T>(0.5) * (current + x / current), current);
	}
}
template<typename T>
inline constexpr T const_sqrt(T x)
{
	return x >= static_cast<T>(0) && x < std::numeric_limits<T>::infinity()
		? detail::const_sqrt_helper(x, x, 0)
		: std::numeric_limits<T>::quiet_NaN();
}
#endif /* BASIC_HPP */
