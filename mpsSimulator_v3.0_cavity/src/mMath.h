/*
*/
///math functions using metaprogramming

#pragma once

namespace mMath {

	template <unsigned n>	struct factorial	{ enum { value = n* factorial<n - 1>::value, }; };
	template <>				struct factorial<0> { enum { value = 1, }; };

	template <unsigned n, unsigned m>
	struct H { enum { value = factorial<n + m - 1>::value / (factorial<m>::value * factorial<n - 1>::value), };	};

}