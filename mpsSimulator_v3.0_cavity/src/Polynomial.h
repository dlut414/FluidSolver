/*
 * LICENSE
*/
///Polynomial.h
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "mMath.h"

namespace SIM {

	using namespace mMath;

	template <unsigned D, unsigned P, unsigned I, unsigned J>
	struct MultisetOfOrderP { enum { value = MultisetOfOrderP<D, P, I-1, J>::value-1, }; };

	template<unsigned D, unsigned P, 0, unsigned J>					struct MultisetOfOrderP { enum { value = 0, }; };
	template<unsigned D, unsigned P, 0, 0>							struct MultisetOfOrderP { enum { value = P, }; };
	template<unsigned D, unsigned P, H<D, P>::value, unsigned J>	struct MultisetOfOrderP { enum { value = 0, }; };
	template<unsigned D, unsigned P, H<D, P>::value, D>				struct MultisetOfOrderP { enum { value = P, }; };

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A {
	public:
		enum { n = H<D, P>::value + Polynomial_A<R, D, P-1>::n, };

		template <typename V, typename U>
		void Gen(const V& v, U& out) const {
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
	};
	template <typename R, unsigned D>
	class Polynomial_B < R, D, 0 > {
	public:
		enum { n = 1, };
	};

}