/*
*/
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "mMath.h"

namespace SIM {

	using namespace mMath;

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A {
	public:
		enum { n = H<D, P>::value + Polynomial_A<R, D, P-1>::n, };
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
		enum { n = H<D, P>::value + Polynomial_B<R, D, P-1>::n, };
		typedef Eigen::Matrix<R, D, 1> vec;
	};
	template <typename R, unsigned D>
	class Polynomial_B < R, D, 0 > {
	public:
		enum { n = 1, };
	};

}