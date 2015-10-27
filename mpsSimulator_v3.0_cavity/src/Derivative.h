/*
* LICENSE
*/
///Derivative.h
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "mMath.h"
#include "Polynomial.h"

namespace mMath {

	template <typename R, unsigned D, unsigned P>
	class Derivative_A_ {
	public:
		template <unsigned N>	static __forceinline R Power(const R& x)	{ return x*Power<N - 1>(x); }
		template <>				static __forceinline R Power<0>(const R& x) { return 1; }
	};

	template <typename R, unsigned D, unsigned P>	class Derivative_A : public Derivative_A_<R,D,P>{};

	template <typename R, unsigned P>				class Derivative_A<R,1,P> : public Derivative_A_<R,1,P> {
	private:
		template <unsigned Der, unsigned P_ = 1, bool Zero = (Der>P_) > 
		struct Gen_					{
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = Der*Power<P_-Der>(x); Gen_<Der,P_+1,(Der>(P_+1))>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = Der*Power<P_-Der>(x/s)*Power<Der>(1/s); Gen_<Der,P_+1,(Der>(P_+1))>::Gen(s, x, out); }
		};
		template <unsigned Der, unsigned P_>
		struct Gen_<Der,P_,true>		{
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = 0; Gen_<Der,P_+1,(Der>(P_+1))>::Gen(x, out); }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = 0; Gen_<Der,P_+1,(Der>(P_+1))>::Gen(s, x, out); }
		};
		template <unsigned Der>
		struct Gen_<Der,P,true>		{
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = 0; }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = 0; }
		};
		template <unsigned Der>
		struct Gen_<Der,P,false>		{
			static __forceinline void Gen(const R& x, R* const out) { out[P_-1] = Der*Power<P_-Der>(x); }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = Der*Power<P_-Der>(x/s)*Power<Der>(1/s); }
		};
	public:
		template <unsigned Der> static __forceinline void Gen(const R& x, R* const out)				{ Gen_<Der>::Gen(x, out); }
		template <unsigned Der> static __forceinline void Gen(const R& s, const R& x, R* const out) { Gen_<Der>::Gen(s, x, out); }
	};

	template <typename R, unsigned P>				class Derivative_A<R,2,P> : public Derivative_A_<R,2,P> {
	private:
		template <unsigned P_ = 1, unsigned PX = P_, unsigned I = 0>
		struct Gen__			{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = Power<PX>(x[0])*Power<P_-PX>(x[1]); Gen__<P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = Power<PX>(x[0]/s)*Power<P_-PX>(x[1]/s); Gen__<P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned P_, unsigned I>						
		struct Gen__<P_,0,I>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = Power<P_>(x[1]); }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = Power<P_>(x[1]/s); }
		};
		template <unsigned P_ = 1, unsigned I = 0, bool Over = !(P-P_)>
		struct Gen_				{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<P_,P_,I>::Gen(x, out); Gen_<P_+1,I+H<2,P_>::value>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<P_,P_,I>::Gen(s, x, out); Gen_<P_+1,I+H<2,P_>::value>::Gen(s, x, out); }
		};
		template <unsigned P_, unsigned I>
		struct Gen_<P_,I,true>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<P_,P_,I>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<P_,P_,I>::Gen(s, x, out); }
		};
	public:
		static __forceinline void Gen(const R* const x, R* const out)				{ Gen_<>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen_<>::Gen(s, x, out); }
	};
	


	template <typename R, unsigned D, unsigned P>
	class Derivative_B_ {
	public:
		enum { value = H<D,P>::value + Derivative_B_<R,D,P-1>::value, };
	};
	template <typename R, unsigned D>
	class Derivative_B_ <R,D,0> {
	public:
		enum { value = 1, };
	};

}