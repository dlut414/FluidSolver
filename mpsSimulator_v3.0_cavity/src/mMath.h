/*
*/
///math functions using metaprogramming

#pragma once

namespace mMath {

	template <unsigned n>	struct Factorial	{ enum { value = n* Factorial<n - 1>::value, }; };
	template <>				struct Factorial<0> { enum { value = 1, }; };

	template <unsigned n, unsigned m>
	struct H { enum { value = Factorial<n + m - 1>::value / (Factorial<m>::value * Factorial<n - 1>::value), };	};

}