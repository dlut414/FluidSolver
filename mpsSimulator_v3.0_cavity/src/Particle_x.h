/*
*/
#pragma once
#include "Header.h"
#include "Particle.h"
#include "Polynomial.h"
#include "Derivative.h"

#define UPWIND_VEL 0

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class Particle_x : public Particle<R,D,Particle_x<R,D,P>> {
		typedef mMath::Polynomial_A<R, D, P> PN;
		typedef mMath::Derivative_A<R, D, P> DR;
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1>		Vec;
		typedef Eigen::Matrix<R, D, D>		Mat;
		typedef Eigen::Matrix<R, PN::value, 1>	VecP;
		typedef Eigen::Matrix<R, PN::value, D>	MatPD;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
	public:
		Particle_x() : Particle() {}
		~Particle_x() {}
		
		__forceinline void poly(const Vec& in, VecP& out) const { PN::Gen(varrho, in.data(), out.data()); }

		const R func(const std::vector<R>& phi, const unsigned& p) const {
			return phi[p];
		}

		const Vec func(const std::vector<Vec>& u, const unsigned& p) const {
			return u[p];
		}

		const Vec grad(const std::vector<R>& phi, const unsigned& p) const {
			VecP  vv = VecP::Zero();
			const iVec c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec ne = c + iVec(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto	dr = pos[q] - pos[p];
							const auto	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto	w = w3(dr1);
							VecP npq;
							poly(dr, npq);
							vv += w * (phi[q] - phi[p]) * npq;
						}
					}
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		const Mat grad(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							VecP npq;
							poly(dr, npq);
							vv += (w* npq)* (u[q] - u[p]).transpose();
						}
					}
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		template <typename T, typename U, typename V>
		const T grad(const U& phi, const V& p) const {
		}

		const R div(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const iVec c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec ne = c + iVec(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto	dr = pos[q] - pos[p];
							const auto	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto	w = w3(dr1);
							VecP npq;
							poly(dr, npq);
							vv += (w* npq)* (u[q] - u[p]);
						}
					}
				}
			}
			const auto a = invMat[p] * vv;
			auto R ret = static_cast<R>(0);
			for (unsigned d = 0; d < D; d++) ret += pn_p_o.block<1,PN::value>(d, 0) * a.block<PN::value>(0, d);
			return ret;
		}

		const R lap(const std::vector<R>& phi, const unsigned& p) const {
			VecP vv = VecP::Zero();
			const iVec c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec ne = c + iVec(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto	dr = pos[q] - pos[p];
							const auto	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto	w = w3(dr1);
							VecP	npq;
							poly(dr, npq);
							vv += w * (phi[q] - phi[p])* npq;
						}
					}
				}
			}
			const VecP a = invMat[p] * vv;
			return (pn_lap_o*a);
		}

		const Vec lap(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const iVec c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec ne = c + iVec(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto	dr = pos[q] - pos[p];
							const auto	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto	w = w3(dr1);
							VecP	npq;
							poly(dr, npq);
							vv += (w * npq) * (u[q] - [p]);
						}
					}
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_lap_o*a).transpose();
		}

		template <typename T>
		const T rot(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							VecP npq;
							poly(dr, npq);
							vv += (w * npq) * (u[q] - u[p]);
						}
					}
				}
			}
			const auto a = invMat[p] * vv;
			const auto der = pn_p_o*a;
			switch (D) {
			case 1:
				return R(0);
				break;
			case 2:
				return R(der(0,1)-der(1,0));
				break;
			case 3:
				return R(0);
				break;
			default:
				return R(0);
			}
		}

		const R func(const std::vector<R>& phi, const Vec& p) const {
			auto rid = 0;
			auto flag = 0;
			auto rr = std::numeric_limits<R>::max();
			auto c = cell->iCoord(p);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							const auto dr = pos[q] - p;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							flag = 1;
							if (dr1 < rr) {
								rr = dr1;
								rid = q;
							}
						}
					}
				}
			}
			if (!flag) return Vec(0.);

			VecP vv = VecP::Zero();
			c = cell->iCoord(pos[rid]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(rid, q)) continue;
#endif
							const auto dr = pos[q] - pos[rid];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv += w * (phi[q] - phi[rid]) * npq;
						}
					}
				}
			}
			auto a = invMat[rid] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			auto gd = Vec(px.dot(a), 0., pz.dot(a));
			return u[rid] + (p - pos[rid])*gd;
		}

		const Vec func(const std::vector<Vec>& u, const Vec& p) const {
			auto rid = 0;
			auto flag = 0;
			auto rr = std::numeric_limits<R>::max();
			auto c = cell->iCoord(p);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							const auto dr = pos[q] - p;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							flag = 1;
							if (dr1 < rr) {
								rr = dr1;
								rid = q;
							}
						}
					}
				}
			}
			if (!flag) return Vec(0.);

			MatPD vv = MatPD::Zero();
			c = cell->iCoord(pos[rid]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(rid, q)) continue;
#endif
							const auto dr = pos[q] - pos[rid];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[rid].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[rid].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[rid].z) * npq;
						}
					}
				}
			}
			auto a = invMat[rid] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			auto gd = Mat(Vec(px.dot(a.block<5, 1>(0, 0)), 0., px.dot(a.block<5, 1>(0, 2))),
				Vec(0.),
				Vec(pz.dot(a.block<5, 1>(0, 0)), 0., pz.dot(a.block<5, 1>(0, 2))));
			return u[rid] + (p - pos[rid])*gd;
		}

		const Vec func_nobd2(const std::vector<Vec>& u, const Vec& p) const {
			return  Vec(0.);
		}

		const Vec func_mafl(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec(0.);
				auto ww = R(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) ww = 1.;
				ret[fp] = ret[fp] / ww;
			}
			return ret[3] - (dmove.mag() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
		}

		const Vec func_mafl_mmt(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec(0.);
				auto ww = R(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) continue;
				ret[fp] = ret[fp] / ww;
			}
			auto ret_mmt = ret[3] - (dmove.mag() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
			auto ret_min = ret[1];
			auto ret_max = ret[1];
			for (int fp = 1; fp <= 4; fp++) {
				if (ret[fp].mag2() < ret_min.mag2()) ret_min = ret[fp];
				if (ret[fp].mag2() > ret_max.mag2()) ret_max = ret[fp];
			}
			if (ret_mmt.mag2() < ret_min.mag2()) ret_mmt = ret_min* ret_mmt.norm();
			if (ret_mmt.mag2() > ret_max.mag2()) ret_mmt = ret_max* ret_mmt.norm();
			return ret_mmt;
		}

		const Vec func_mls_a(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {
			Matpp mm = Matpp::Zero();
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			Matpp inv = Matpp::Zero();
			if (abs(mm.determinant()) < eps_Mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_Mat) {
					inv = Matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<5, 1>(0, 0);
			const auto az = a.block<5, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			const auto dp = p_new - pos[p];
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az));
			return ret;
		}

		const Vec func_mls_a_upwind(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dp.norm();
#endif
			Matpp mm = Matpp::Zero();
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							if (dr*up < 0) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			Matpp inv = Matpp::Zero();
			if (abs(mm.determinant()) < eps_Mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_Mat) {
					inv = Matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<5, 1>(0, 0);
			const auto az = a.block<5, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az));
			return ret;
		}

		const Vec func_mls_a_downwind(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dp.norm();
#endif
			Matpp mm = Matpp::Zero();
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							if (dr*up > 0) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 5>(0, 0) += w * npq[0] * npq;
							mm.block<1, 5>(1, 0) += w * npq[1] * npq;
							mm.block<1, 5>(2, 0) += w * npq[2] * npq;
							mm.block<1, 5>(3, 0) += w * npq[3] * npq;
							mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							vv.block<5, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<5, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<5, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			Matpp inv = Matpp::Zero();
			if (abs(mm.determinant()) < eps_Mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_Mat) {
					inv = Matpp::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<5, 1>(0, 0);
			const auto az = a.block<5, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az));
			return ret;
		}

		const Vec func_mls_b(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}

		const Vec func_mls_b_upwind(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				Matpp mm = Matpp::Zero();
				const auto c = cell->iCoord(pos[p]);
				for (auto k = -1; k <= 1; k++) {
					for (auto j = -1; j <= 1; j++) {
						for (auto i = -1; i <= 1; i++) {
							const auto ne = c + iVec(i, j, k);
							const auto key = cell->hash(ne);
							for (auto m = 0; m < cell->linkList[key].size(); m++) {
								const auto q = cell->linkList[key][m];
#if BD_OPT
								if (bdOpt(p, q)) continue;
#endif
								const auto dr = pos[q] - pos[p];
								const auto dr1 = dr.mag();
								if (dr1 > r0) continue;
								const auto w = w3(dr1);
								const auto npq = poly(dr);
								mm.block<1, 5>(0, 0) += w * npq[0] * npq;
								mm.block<1, 5>(1, 0) += w * npq[1] * npq;
								mm.block<1, 5>(2, 0) += w * npq[2] * npq;
								mm.block<1, 5>(3, 0) += w * npq[3] * npq;
								mm.block<1, 5>(4, 0) += w * npq[4] * npq;
							}
						}
					}
				}
				invMat[p] = Matpp::Zero();
				if (abs(mm.determinant()) < eps_Mat) {
#if DEBUG
					std::cout << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_Mat) {
						invMat[p] = Matpp::Zero();
						continue;
					}
					invMat[p].block<2, 2>(0, 0) = mm_.inverse();
					continue;
				}
				invMat[p] = mm.inverse();
			}
		}

		void updateDiver() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(np); p++) {
				diver[p] = div(part->vel1, p);
			}
		}

		void init_x() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(MatPP());
			}

			varrho = 1.*dp;
			Vec zero = Vec::Zero();
			switch (D) {
			case 1:
				DR::Gen<1>(zero.data(), pn_p_o.data());
				DR::Gen<2>(zero.data(), pn_pp_o.data());
				DR::Gen<2>(zero.data(), pn_lap_o.data());
				break;
			case 2:
				DR::Gen<1,0>(zero.data(), pn_p_o.block<1,PN::value>(0, 0).data());
				DR::Gen<0,1>(zero.data(), pn_p_o.block<1,PN::value>(1, 0).data());
				DR::Gen<2,0>(zero.data(), pn_pp_o.block<1,PN::value>(0, 0).data());
				DR::Gen<1,1>(zero.data(), pn_pp_o.block<1,PN::value>(1, 0).data());
				DR::Gen<0,2>(zero.data(), pn_pp_o.block<1,PN::value>(2, 0).data());
				pn_lap_o = pn_pp_o.block<1,PN::value>(0, 0) + pn_pp_o.block<1,PN::value>(2, 0);
				break;
			case 3:
				//DR::Gen<1,0,0>(zero.data(), pn_p_o.block<1,PN::value>(0, 0).data());
				//DR::Gen<0,1,0>(zero.data(), pn_p_o.block<1,PN::value>(1, 0).data());
				//DR::Gen<0,0,1>(zero.data(), pn_p_o.block<1,PN::value>(2, 0).data());
				//DR::Gen<2,0,0>(zero.data(), pn_pp_o.block<1,PN::value>(0, 0).data());
				//DR::Gen<1,1,0>(zero.data(), pn_pp_o.block<1,PN::value>(1, 0).data());
				//DR::Gen<1,0,1>(zero.data(), pn_pp_o.block<1,PN::value>(2, 0).data());
				//DR::Gen<0,2,0>(zero.data(), pn_pp_o.block<1,PN::value>(3, 0).data());
				//DR::Gen<0,1,1>(zero.data(), pn_pp_o.block<1,PN::value>(4, 0).data());
				//DR::Gen<0,0,2>(zero.data(), pn_pp_o.block<1,PN::value>(5, 0).data());
				//pn_lap_o = pn_pp_o.block<1,PN::value>(0, 0) + pn_pp_o.block<1,PN::value>(3, 0) + pn_pp_o.block<1,PN::value>(5, 0);
				break;
			default:
				break;
			}
		}

	public:
		std::vector<MatPP> invMat;

		R varrho;
		Eigen::Matrix<R,D,PN::value>				pn_p_o;
		Eigen::Matrix<R,mMath::H<D,2>,PN::value>	pn_pp_o;
		Eigen::Matrix<R,1,PN::value>				pn_lap_o;
	};

}