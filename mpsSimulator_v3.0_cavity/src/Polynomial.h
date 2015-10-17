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

	template <unsigned D, unsigned J, unsigned ROrder> struct Repetition { enum { value = H<D-j-1,m>, }; };

	template <unsigned D, unsigned I, unsigned J, unsigned ROrder, bool Clear>	struct CurRepetition						{ enum { value, }; };
	template <unsigned D, unsigned J, unsigned ROrder, bool Clear>				struct CurRepetition <D,0,J,ROrder,Clear>	{ enum { value = 1, }; };
	template <unsigned D, unsigned I, unsigned J, unsigned ROrder>				struct CurRepetition <D,I,J,ROrder,1> { enum { value = 1, }; };
	template <unsigned D, unsigned I, unsigned J, unsigned ROrder>				struct CurRepetition <D,I,J,ROrder,0> { enum { value = CurRepetition<D,I-1,J,ROrder,CurRepetition<D,I-1,J,ROrder,0>::value>Repetition<D,J,ROrder>::value>::value + 1, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool Rep> struct Multiset { enum { value, }; };

	template<unsigned D, unsigned P, 0, unsigned J, unsigned ROrder>	struct Multiset { enum { value = 0, }; };
	template<unsigned D, unsigned P, 0, 0, unsigned ROrder>			struct Multiset { enum { value = P, }; };

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