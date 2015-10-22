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

	template <unsigned D, unsigned J, unsigned ROrder> struct Repetition { enum { value = H<D-j-1,ROrder>, }; };

	template <unsigned D, unsigned I, unsigned J, unsigned ROrder, bool Clear>	struct CurRepetition						{ enum { value, }; };
	template <unsigned D, unsigned J, unsigned ROrder, bool Clear>				struct CurRepetition<D,0,J,ROrder,Clear>	{ enum { value = 1, }; };
	template <unsigned D, unsigned I, unsigned J, unsigned ROrder>				struct CurRepetition<D,I,J,ROrder,1>		{ enum { value = 1, }; };
	template <unsigned D, unsigned I, unsigned J, unsigned ROrder>				struct CurRepetition<D,I,J,ROrder,0>		{ enum { value = CurRepetition<D,I-1,J,ROrder,CurRepetition<D,I-1,J,ROrder,0>::value>Repetition<D,J,ROrder>::value>::value + 1, }; };
	template <unsigned D, unsigned I, unsigned J, unsigned ROrder>
	struct CurRepetitionIns { enum { value = CurRepetition<D,I,J,ROrder,CurRepetition<D,I,J,ROrder,0>::value>Repetition<D,J,ROrder>::value>, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool RepOver>	struct Multiset { enum { value, }; };
	template <unsigned D, unsigned P, unsigned J, unsigned ROrder, bool RepOver>				struct Multiset<D,P,0,J,ROrder,RepOver> { enum { value = 0, }; };
	template <unsigned D, unsigned P, unsigned ROrder, bool RepOver>							struct Multiset<D,P,0,0,ROrder,RepOver> { enum { value = P, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder>
	struct Multiset<D,P,I,J,ROrder,0> { enum { value = Multiset<D,P,I-1,J,ROrder,CurRepetitionIns<D,I-1,J,ROrder>::value>=Repetition<D,J,ROrder>::value&&CurRepetitionIns<D,I-1,J,ROrder>::value!=0>::value, }; };
	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder>
	struct Multiset<D,P,I,J,ROrder,1> { enum { value = Multiset<D,P,I-1,J,ROrder,CurRepetitionIns<D,I-1,J,ROrder>::value>=Repetition<D,J,ROrder>::value&&CurRepetitionIns<D,I-1,J,ROrder>::value!=0>::value - 1, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, int Left, bool RepOver>
	struct MultisetUp { enum { value, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool RepOver>
	struct MultisetUp<D,P,I,J,ROrder,0,RepOver> { enum { value = 0, }; };
	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool RepOver>
	struct MultisetUp<D,P,I,J,ROrder,1,RepOver> { enum { value = ROrder + 1, }; };

	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool Left, bool RepOver>
	struct MultisetUp<D,P,I,J,ROrder,Left,1,RepOver> { enum { value = 0, }; };
	template <unsigned D, unsigned P, unsigned I, unsigned J, unsigned ROrder, bool Left, bool RepOver>
	struct MultisetUp<D,P,I,J,ROrder,Left,0,RepOver> { enum { value = 0, }; };

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