/*
*/
#pragma once
#include <vector>
#include <Eigen/Dense>

namespace SIM {

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

	template <unsigned i, unsigned j>
	class coefTable {
	public:
		enum {value = 1,};
	};

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A {
	public:
		enum { n = H<D, P>::value + Polynomial_A<R,D,P-1>::n, };
		typedef Eigen::Matrix<R, D, 1> vec;
		typedef Eigen::Matrix<R, n, 1> retVec;
		const retVec genPoly(const vec& v) const {
			retVec ret;
		}
	};
	template <typename R, unsigned D>
	class Polynomial_A < R, D, 0 > {
	public:
		enum { n = 0, };
	};

	template <typename R, unsigned D, unsigned P>
	class Polynomial_B {
	public:
		enum { n = H<D, P>::value + Polynomial_B<R, D, P - 1>::n, };
		typedef Eigen::Matrix<R, D, 1> vec;
	};
	template <typename R, unsigned D>
	class Polynomial_B < R, D, 0 > {
	public:
		enum { n = 1, };
	};

}