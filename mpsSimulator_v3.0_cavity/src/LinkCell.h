/*
*/
#pragma once
#include <iostream>
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include "Header.h"

namespace SIM {

	template <typename R, unsigned D>
	class LinkCell {
		typedef Eigen::Matrix<R,D,1>	vecD;
		typedef Eigen::Matrix<int,D,1>	iVecD;
	public:
		template <typename T>
		LinkCell(const BBox<T,D>& b, const R& c) : box(b), cSize(c) { init(); }
		~LinkCell() { fina(); }

		template <unsigned D_ = 0, bool Over = (D_==(D-1))> struct Convert
		{ template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out) { out[D_] = static_cast<V>(in[D_]); Convert<D_+1>::Gen(in,out); } };
		template <unsigned D_>								struct Convert<D_,true>
		{ template <typename U, typename V> static __forceinline void Gen(const U* const in, V* const out) { out[D_] = static_cast<V>(in[D_]); } };

		template <unsigned D_ = D>	__forceinline const unsigned hash		(const iVecD& c)	const {}
		template <>					__forceinline const unsigned hash<1>	(const iVecD& c)	const { return unsigned((c.x + 100000) % cNum); }
		template <>					__forceinline const unsigned hash<2>	(const iVecD& c)	const { return unsigned((c.x + dv.x* c.y + 100000) % cNum); }
		template <>					__forceinline const unsigned hash<3>	(const iVecD& c)	const { return unsigned((c.x + dv.x* c.y + sheet* c.z + 100000) % cNum); }
		
		__forceinline const unsigned hash(const vecD& p) const { return hash( iCoord(p) ); }
		__forceinline const iVecD iCoord(const vecD& p) const { iVecD ret; Convert<>::Gen(p.data(), ret.data()); return ret; }

		void update(const std::vector<vecD>& pos) {
			for (unsigned i = 0; i < linkList.size(); i++) {
				linkList[i].clear();
			}
			for (unsigned i = 0; i < pos.size(); i++) {
				linkList[hash(pos[i])].push_back(i);
			}
		}

	public:
		BBox<R,D> box;
		std::vector< std::vector<unsigned> > linkList;

	private:
		void init() {
			vecD v = vecD(box.pMax - box.pMin);
			for (int d = 0; d < int(D); d++) { dv[d] = int(v[d] / cSize) + 1; dv[d] = dv[d] > 3 ? dv[d] : 3; }
			cNum = 1;
			for (int d = 0; d < int(D); d++) { cNum *= unsigned(dv[d]); }
			if(D==3) sheet = unsigned(dv[0] * dv[1]);
			else sheet = 1;

			if (cNum >= linkList.max_size()) {
				std::cout << " linkList overflow ! " << std::endl;
				exit(0);
			}
			linkList.clear();
			linkList = std::vector< std::vector<unsigned> >(cNum);
			std::cout << " cell num: " << cNum << std::endl;
		}
		void fina() { linkList.clear(); }

	private:
		iVecD dv;
		unsigned sheet;
		unsigned cNum;
		R cSize;
	};

}