/*
*/
///math functions using metaprogramming

#pragma once

namespace mMath {

	template <unsigned n>
	class factorial {
	public:
		enum { value = n* factorial<n - 1>::value, };
	};

	template <>
	class factorial < 1 > {
	public:
		enum { value = 1, };
	};

	template <unsigned n, unsigned m>
	class H {
	public:
		enum { value = factorial<n + m - 1>::value / (factorial<m>::value * factorial<n - 1>::value), };
	};

}