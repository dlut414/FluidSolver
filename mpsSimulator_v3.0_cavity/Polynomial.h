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

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A_ {
	public:
		enum { value = H<D, P>::value + Polynomial_A_<R,D,P-1>::value, };
		template <unsigned N>	static __forceinline R Power(const R& x) { return x*Power<N - 1>(x); }
		template <>				static __forceinline R Power<0>(const R& x) { return 1; }
	};
	template <typename R, unsigned D>	class Polynomial_A_<R,D,0> { public: enum { value = 0, }; };

	template <typename R, unsigned D, unsigned P>	class Polynomial_A : public Polynomial_A_<R,D,P>{};

	template <typename R, unsigned P>				class Polynomial_A<R,1,P> : public Polynomial_A_<R,1,P> {
	private:
		template <unsigned P_ = 1>	static __forceinline void Gen_		(const R& x, R* const out) { out[P_-1] = Power<P_>(x); Gen_<P_+1>(x, out); }
		template <>					static __forceinline void Gen_<P>	(const R& x, R* const out) { out[P-1] = Power<P>(x); }
	public:
		static __forceinline void Gen(const R& x, R* const out) { Gen_(x, out); }
	};

	template <typename R, unsigned P>				class Polynomial_A<R,2,P> : public Polynomial_A_<R,2,P> {
	private:
		template <unsigned P_ = 1, unsigned PX = P_, unsigned I = 0>	static __forceinline void Gen_		(const R* const x, R* const out) { out[I] = Power<PX>(x[0])*Power<P_-PX>(x[1]); Gen_<P_+1>(x, out); }
		template <>										static __forceinline void Gen_<P>	(const R* const x, R* const out) { out[P-1] = Power<P>(x); }
	public:
		static __forceinline void Gen(const R* const x, R* const out) { Gen_(x, out); }
	};
	


	template <typename R, unsigned D, unsigned P>
	class Polynomial_B_ {
	public:
		enum { value = H<D, P>::value + Polynomial_B_<R, D, P-1>::value, };
	};
	template <typename R, unsigned D>
	class Polynomial_B_ < R, D, 0 > {
	public:
		enum { value = 1, };
	};

}