/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//Interpolation.h
///defination of class Interpolation
#pragma once

namespace SIM {

	template <typename R, int D, int P>
	class Interpolation {
		typedef Particle_x<R,D,P> Ptc;
	public:
		Interpolation() {}
		~Interpolation() {}

//		template <int Stencils = 8, int Dimension = D>	struct interpolateWENO_ {};
//		template <int Stencils>							struct interpolateWENO_<Stencils, 2> {
//			template <typename U> 
//			static const U Gen(const Ptc* const part, const std::vector<U>& phi, const unsigned& p, const Vec& p_new) {
//				const auto dp = p_new - part->pos[p];
//				const auto up = dp.normalized();
//				static const auto alpha = 2.* M_PI / Stencils;
//				Ptc::Vec dir[Stencils];
//				for (auto i = 0; i < Stencils; i++) {
//					const auto theta = i* alpha;
//					const auto ct = cos(theta);
//					const auto st = sin(theta);
//					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
//				}
//				Ptc::MatPP mm[Stencils];
//				Ptc::VecP vv[Stencils];
//				for (auto i = 0; i < Stencils; i++) {
//					mm[i] = MatPP::Zero();
//					vv[i] = VecP::Zero();
//				}
//				const auto c = part->cell->iCoord(pos[p]);
//				for (auto i = 0; i < part->cell->blockSize::value; i++) {
//					const auto key = part->cell->hash(c, i);
//					for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
//						const auto q = part->cell->linkList[key][m];
//#if BD_OPT
//						if (bdOpt(p, q)) continue;
//#endif
//						const auto dr = part->pos[q] - part->pos[p];
//						const auto dr1 = dr.norm();
//						if (dr1 > part->r0) continue;
//						for (auto stc = 0; stc < Stencils; stc++) {
//							if (dr.dot(dir[stc]) < 0) continue;
//							const auto w = w3(dr1);
//							Ptc::VecP npq;
//							part->poly(dr, npq);
//							mm[stc] += (w* npq)* npq.transpose();
//							vv[stc] += (w* npq)* (phi[q] - phi[p]);
//						}
//					}
//				}
//				R oscillationIndicator[Stencils];
//				R stencilWeight[Stencils];
//				R stencilWeightNorm[Stencils];
//				Ptc::VecP polyCoef[Stencils];
//				for (auto i = 0; i < Stencils; i++) {
//					Ptc::MatPP inv = Ptc::MatPP::Zero();
//					if (abs(mm.determinant()) < part->eps_mat) {
//						auto mm_ = mm.block<2, 2>(0, 0);
//						if (abs(mm_.determinant()) < part->eps_mat) {
//							inv = Ptc::MatPP::Zero();
//						}
//						else inv.block<2, 2>(0, 0) = mm_.inverse();
//					}
//					else inv = mm[i].inverse();
//					polyCoef[i] = inv * vv[i];
//					const auto offset = PN::value - mMath::H<D, P>::value;
//					oscillationIndicator[i] = R(0.);
//					for (auto term = offset; term < Ptc::PN::value; term++) {
//						oscillationIndicator[i] += abs(polyCoef[i][term]);
//					}
//				}
//				static const R epsilon = 1.e-5;
//				static const int magnifier = 4;
//				for (auto i = 0; i < Stencils; i++) {
//					stencilWeight[i] = 1. / pow(epsilon + oscillationIndicator[i], magnifier);
//				}
//				R stencilWeightSum = R(0.);
//				for (auto i = 0; i < Stencils; i++) {
//					stencilWeightSum += stencilWeight[i];
//				}
//				for (auto i = 0; i < Stencils; i++) {
//					stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
//				}
//				auto combinedCoef = Ptc::VecP::Zero();
//				for (auto i = 0; i < Stencils; i++) {
//					combinedCoef += stencilWeightNorm[i] * polyCoef[i];
//				}
//				const auto gd = pn_p_o * combinedCoef;
//				const auto mgd = pn_pp_o * combinedCoef;
//				Ptc::Mat hes;
//				int counter = 0;
//				for (auto i = 0; i < D; i++) {
//					for (auto j = i; j < D; j++) {
//						hes(i, j) = mgd(counter++);
//					}
//				}
//				for (auto i = 0; i < D; i++) {
//					for (auto j = 0; j < i; j++) {
//						hes(i, j) = hes(j, i);
//					}
//				}
//				const auto dpt = dp.transpose();
//				auto ret = phi[p];
//				ret = ret + (dpt*gd) + 0.5*dpt * hes * dp;
//				return ret;
//			}
//		};
//		template <int Stencils>							struct interpolateWENO_<Stencils, 3> {
//			template <typename U> 
//			static const U Gen(const Ptc* const part, const std::vector<U>& phi, const unsigned& p, const Vec& p_new) {}
//		};
	};

}