/*
*/
#pragma once
#include "Header.h"
#include "Particle.h"
#include "Polynomial.h"
#include "Derivative.h"

namespace SIM {

	template <typename R, int D, int P>
	class Particle_x : public Particle<R,D,Particle_x<R,D,P>> {};

	template <typename R, int P>
	class Particle_x<R,2,P> : public Particle<R,2,Particle_x<R,2,P>> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef EiRun::Matrix<int,2,1> iVec;
		typedef EiRun::Matrix<R,2,1> Vec;
		typedef EiRun::Matrix<R,2,2> Mat;
		typedef EiRun::Matrix<R,PN::value,1> VecP;
		typedef EiRun::Matrix<R,PN::value,2> MatPD;
		typedef EiRun::Matrix<R,PN::value,PN::value> MatPP;
	public:
		Particle_x() : Particle() {}
		~Particle_x() {}

		__forceinline void poly(const R* in, R* out) const { PN::Run(varrho, in, out); }

		const R DerX(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_px_o * aa);
		}

		const R DerY(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_py_o * aa);
		}

		const R DerXX(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_pxx_o * aa);
		}

		const R DerYY(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_pyy_o * aa);
		}

		const Vec Grad(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_p_o * aa);
		}

		const R Div(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vvx += w * (phix[q] - phix[p]) * npq;
					vvy += w * (phiy[q] - phiy[p]) * npq;
				}
			}
			const VecP aax = invMat[p] * vvx;
			const VecP aay = invMat[p] * vvy;
			return (pn_px_o * aax + pn_py_o * aay);
		}

		const R Lap(const R* const phi, const int& p) const {
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const VecP aa = invMat[p] * vv;
			return (pn_lap_o * aa);
		}

		const Vec Lap(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vvx += w * (phix[q] - phix[p]) * npq;
					vvy += w * (phiy[q] - phiy[p]) * npq;
				}
			}
			const VecP aax = invMat[p] * vvx;
			const VecP aay = invMat[p] * vvy;
			Vec ret;
			ret[0] = pn_lap_o * aax;
			ret[1] = pn_lap_o * aay;
			return ret;
		}

		const R Rot(const R* const phix, const R* const phiy, const int& p) const {
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					vvx += w * (phix[q] - phix[p]) * npq;
					vvy += w * (phiy[q] - phiy[p]) * npq;
				}
			}
			const VecP aax = invMat[p] * vvx;
			const VecP aay = invMat[p] * vvy;
			return (pn_py_o * aax - pn_px_o * aay);
		}

		const R interpolateLSA(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					mm += (w * npq) * npq.transpose();
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					return phi[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aa = inv * vv;
			const R Px = pn_px_o* aa;
			const R Py = pn_py_o* aa;
			const R Pxx = pn_pxx_o* aa;
			const R Pxy = pn_pxy_o* aa;
			const R Pyy = pn_pyy_o* aa;
			return phi[p] + (dx*Px + dy*Py) + R(0.5)* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		const Vec interpolateLSA(const R* const phix, const R* const phiy, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					mm += (w * npq) * npq.transpose();
					vvx += w * (phix[q] - phix[p]) * npq;
					vvy += w * (phiy[q] - phiy[p]) * npq;
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					Vec retv;
					retv[0] = phix[p];
					retv[1] = phiy[p];
					return retv;
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aax = inv * vvx;
			const VecP aay = inv * vvy;
			const R Px[2] = { pn_px_o* aax, pn_px_o* aay };
			const R Py[2] = { pn_py_o* aax, pn_py_o* aay };
			const R Pxx[2] = { pn_pxx_o* aax, pn_pxx_o* aay };
			const R Pxy[2] = { pn_pxy_o* aax, pn_pxy_o* aay };
			const R Pyy[2] = { pn_pyy_o* aax, pn_pyy_o* aay };
			const Vec ret;
			ret[0] = phix[p] + (dx*Px[0] + dy*Py[0]) + R(0.5)* (dx*dx*Pxx[0] + R(2)*dx*dy*Pxy[0] + dy*dy*Pyy[0]);
			ret[1] = phix[p] + (dx*Px[1] + dy*Py[1]) + R(0.5)* (dx*dx*Pxx[1] + R(2)*dx*dy*Pxy[1] + dy*dy*Pyy[1]);
			return ret;
		}

		const R interpolateLSAU(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0 || (dx*dr[0] + dy*dr[1]) < R(0)) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					mm += (w * npq) * npq.transpose();
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					return phi[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aa = inv * vv;
			const R Px = pn_px_o* aa;
			const R Py = pn_py_o* aa;
			const R Pxx = pn_pxx_o* aa;
			const R Pxy = pn_pxy_o* aa;
			const R Pyy = pn_pyy_o* aa;
			return phi[p] + (dx*Px + dy*Py) + R(0.5)* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		const Vec interpolateLSAU(const R* const phix, const R* const phiy, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			MatPP mm = MatPP::Zero();
			VecP vvx = VecP::Zero();
			VecP vvy = VecP::Zero();
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0 || (dx*dr[0] + dy*dr[1]) < R(0)) continue;
					const R w = ww(dr1);
					VecP npq;
					poly(dr, npq.data());
					mm += (w * npq) * npq.transpose();
					vvx += w * (phix[q] - phix[p]) * npq;
					vvy += w * (phiy[q] - phiy[p]) * npq;
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					Vec retv;
					retv[0] = phix[p];
					retv[1] = phiy[p];
					return retv;
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const VecP aax = inv * vvx;
			const VecP aay = inv * vvy;
			const R Px[2] = { pn_px_o* aax, pn_px_o* aay };
			const R Py[2] = { pn_py_o* aax, pn_py_o* aay };
			const R Pxx[2] = { pn_pxx_o* aax, pn_pxx_o* aay };
			const R Pxy[2] = { pn_pxy_o* aax, pn_pxy_o* aay };
			const R Pyy[2] = { pn_pyy_o* aax, pn_pyy_o* aay };
			const Vec ret;
			ret[0] = phix[p] + (dx*Px[0] + dy*Py[0]) + R(0.5)* (dx*dx*Pxx[0] + R(2)*dx*dy*Pxy[0] + dy*dy*Pyy[0]);
			ret[1] = phix[p] + (dx*Px[1] + dy*Py[1]) + R(0.5)* (dx*dx*Pxx[1] + R(2)*dx*dy*Pxy[1] + dy*dy*Pyy[1]);
			return ret;
		}

		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY>
		const R interpolateWENO_A_(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			const R dd = sqrt(dx*dx + dy*dy);
			if (dd < eps) return phi[p];
			const R upx = dx / dd;
			const R upy = dy / dd;
			const R alpha = 2.* M_PI / StencilsX;
			Vec dir[StencilsX];
			Vec ctr[Stencils];
			for (int i = 0; i < StencilsX; i++) {
				const R theta = i* alpha;
				const R ct = cos(theta);
				const R st = sin(theta);
				dir[i] << ct*upx + st*upy, ct*upx - st*upy;
			}
			for (int j = 0; j < StencilsY; j++) {
				//const R dis = r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
				const R dis = r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
				for (int i = 0; i < StencilsX; i++) {
					const int stcId = i* StencilsY + j;
					ctr[stcId][0] = pos[0][p] + dis*dir[i][0];
					ctr[stcId][1] = pos[1][p] + dis*dir[i][1];
				}
			}
			MatPP mm[Stencils];
			VecP vv[Stencils];
			for (int i = 0; i < Stencils; i++) {
				mm[i] = MatPP::Zero();
				vv[i] = VecP::Zero();
			}
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					for (int stcId = 0; stcId < Stencils; stcId++) {
						const R disx = (pos[0][q] - ctr[stcId][0]);
						const R disy = (pos[1][q] - ctr[stcId][1]);
						const R dis = sqrt(disx*disx + disy*disy);
						const R w = ww(dis);
						VecP npq;
						part->poly(dr, npq.data());
						mm[stcId] += (w* npq)* npq.transpose();
						vv[stcId] += (w* npq)* (phi[q] - phi[p]);
					}
				}
			}
			R oscillationIndicator[Stencils];
			R stencilWeight[Stencils];
			R stencilWeightNorm[Stencils];
			VecP polyCoef[Stencils];
			for (int i = 0; i < Stencils; i++) {
				MatPP inv = MatPP::Zero();
				if (abs(mm[i].determinant()) < part->eps_mat) {
					auto mm_ = mm[i].block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < part->eps_mat) {
						inv = MatPP::Zero();
						inv(0, 0) = R(1) / eps;
						inv(1, 1) = R(1) / eps;
					}
					else inv.block<2, 2>(0, 0) = mm_.inverse();
				}
				else inv = mm[i].inverse();
				polyCoef[i] = inv * vv[i];
				const int offset = PN::value - mMath::H<D, P>::value;
				oscillationIndicator[i] = R(0.);
				for (int term = offset; term < PN::value; term++) {
					oscillationIndicator[i] += abs(polyCoef[i][term]);
				}
			}
			const R epsilon = 1.e-6;
			const int magnifier = 5;
			for (int i = 0; i < Stencils; i++) {
				stencilWeight[i] = R(1) / pow(epsilon + oscillationIndicator[i], magnifier);
			}
			R stencilWeightSum = R(0);
			for (int i = 0; i < Stencils; i++) {
				stencilWeightSum += stencilWeight[i];
			}
			for (int i = 0; i < Stencils; i++) {
				stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
			}
			VecP combinedCoef = VecP::Zero();
			for (int i = 0; i < Stencils; i++) {
				combinedCoef += stencilWeightNorm[i] * polyCoef[i];
			}
			const R Px = pn_px_o* combinedCoef;
			const R Py = pn_py_o* combinedCoef;
			const R Pxx = pn_pxx_o* combinedCoef;
			const R Pxy = pn_pxy_o* combinedCoef;
			const R Pyy = pn_pyy_o* combinedCoef;
			return phi[p] + (dx*Px + dy*Py) + R(0.5)* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY>
		const R interpolateWENO_B_(const R* const phi, const int& p, const R& px, const R& py) const {
			const R dx = px - pos[0][p];
			const R dy = py - pos[1][p];
			const R dd = sqrt(dx*dx + dy*dy);
			if (dd < eps) return phi[p];
			const R upx = dx / dd;
			const R upy = dy / dd;
			const R alpha = 2.* M_PI / StencilsX;
			Vec dir[StencilsX];
			Vec ctr[Stencils];
			for (int i = 0; i < StencilsX; i++) {
				const R theta = i* alpha;
				const R ct = cos(theta);
				const R st = sin(theta);
				dir[i] << ct*upx + st*upy, ct*upx - st*upy;
			}
			for (int j = 0; j < StencilsY; j++) {
				//const R dis = r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
				const R dis = r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
				for (int i = 0; i < StencilsX; i++) {
					const int stcId = i* StencilsY + j;
					ctr[stcId][0] = pos[0][p] + dis*dir[i][0];
					ctr[stcId][1] = pos[1][p] + dis*dir[i][1];
				}
			}
			MatPP mm[Stencils];
			VecP vv[Stencils];
			for (int i = 0; i < Stencils; i++) {
				mm[i] = MatPP::Zero();
				vv[i] = VecP::Zero();
			}
			const int cx = cell->pos2cell(pos[0][p]);
			const int cy = cell->pos2cell(pos[1][p]);
			for (int i = 0; i < cell->blockSize::value; i++) {
				const int key = cell->hash(cx, cy, i);
				for (int m = 0; m < cell->linkList[key].size(); m++) {
					const int q = cell->linkList[key][m];
					if (type[p] == BD2) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 > r0) continue;
					for (int stcId = 0; stcId < Stencils; stcId++) {
						const R disx = (pos[0][q] - ctr[stcId][0]);
						const R disy = (pos[1][q] - ctr[stcId][1]);
						const R dis = sqrt(disx*disx + disy*disy);
						const R w = ww(dis);
						VecP npq;
						part->poly(dr, npq.data());
						mm[stcId] += (w* npq)* npq.transpose();
						vv[stcId] += (w* npq)* (phi[q] - phi[p]);
					}
				}
			}
			R oscillationIndicator[Stencils];
			R stencilWeight[Stencils];
			R stencilWeightNorm[Stencils];
			VecP polyCoef[Stencils];
			for (int i = 0; i < Stencils; i++) {
				MatPP inv = MatPP::Zero();
				if (abs(mm[i].determinant()) < part->eps_mat) {
					return phi[p];
					auto mm_ = mm[i].block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < part->eps_mat) {
						inv = MatPP::Zero();
						inv(0, 0) = R(1) / eps;
						inv(1, 1) = R(1) / eps;
					}
					else inv.block<2, 2>(0, 0) = mm_.inverse();
				}
				else inv = mm[i].inverse();
				polyCoef[i] = inv * vv[i];
				oscillationIndicator[i] = R(0.);
			}

			for (auto i = 0; i < Stencils; i++) {
				const R A = polyCoef[i][0] * polyCoef[i][0];
				const R B = polyCoef[i][1] * polyCoef[i][1];
				const R C = polyCoef[i][2] * polyCoef[i][2];
				const R D = polyCoef[i][3] * polyCoef[i][3];
				const R E = polyCoef[i][4] * polyCoef[i][4];
				const R beta1 = (dp*dp*dp)*(A + B) + (dp*dp*dp*dp*dp / R(6))*(R(2)*C + D + R(2)*E);
				const R beta2 = (dp*dp*dp*dp*dp)*(R(4)*C + D + R(4)*E);
				//oscillationIndicator[i] = R(4)*beta1 - (R(1) / R(3))*beta2;
				oscillationIndicator[i] = beta1 + beta2;
			}
			const R epsilon = 1.e-6;
			const int magnifier = 5;
			for (int i = 0; i < Stencils; i++) {
				stencilWeight[i] = R(1) / pow(epsilon + oscillationIndicator[i], magnifier);
			}
			R stencilWeightSum = R(0);
			for (int i = 0; i < Stencils; i++) {
				stencilWeightSum += stencilWeight[i];
			}
			for (int i = 0; i < Stencils; i++) {
				stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
			}
			VecP combinedCoef = VecP::Zero();
			for (int i = 0; i < Stencils; i++) {
				combinedCoef += stencilWeightNorm[i] * polyCoef[i];
			}
			const R Px = pn_px_o* combinedCoef;
			const R Py = pn_py_o* combinedCoef;
			const R Pxx = pn_pxx_o* combinedCoef;
			const R Pxy = pn_pxy_o* combinedCoef;
			const R Pyy = pn_pyy_o* combinedCoef;
			return phi[p] + (dx*Px + dy*Py) + R(0.5)* (dx*dx*Pxx + R(2)*dx*dy*Pxy + dy*dy*Pyy);
		}

		__forceinline const R interpolateWENO(const R* const phi, const int& p, const R& px, const R& py) const {
			return interpolateWENO_B_<>(phi, p, px, py);
		}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				MatPP mm = MatPP::Zero();
				const int cx = cell->pos2cell(pos[0][p]);
				const int cy = cell->pos2cell(pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (type[p] == BD2) continue;
						const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > r0) continue;
						const R w = ww(dr1);
						VecP npq;
						poly(dr, npq.data());
						mm += (w* npq) * npq.transpose();
					}
				}

				if (type[p] == BD1) {
					const Vec& inner = bdnorm.at(p);
					MatPP nn = MatPP::Zero();
					nn.block<D, D>(0, 0) = w3(0)* inner * inner.transpose();
					MatPP mpn = mm + nn;
					MatPP& inv = invNeu.at(p);
					inv = MatPP::Zero();
					if (abs(mpn.determinant()) < eps_mat) {
#if DEBUG
						std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mpn.determinant() << std::endl;
#endif
						auto mpn_ = mpn.block<2, 2>(0, 0);
						if (abs(mpn_.determinant()) < eps_mat) inv = MatPP::Zero();
						else inv.block<2, 2>(0, 0) = mpn_.inverse();
					}
					else inv = mpn.inverse();
				}

				MatPP& invRef = invMat.at(p);
				invRef = MatPP::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) invRef = MatPP::Zero();
					else invRef.block<2, 2>(0, 0) = mm_.inverse();
				}
				else invRef = mm.inverse();
			}
		}

		void init_x() {
			invNeu.clear();
			invMat.clear();
			for (int p = 0; p < np; p++) {
				invMat.push_back(MatPP());
				if (type[p] == BD1) invNeu[p] = MatPP::Zero();
			}
			varrho = 1. / (1.*dp);
			Vec zero = Vec::Zero();
			DR::Run<1, 0>(varrho, zero.data(), pn_px_o.data());
			DR::Run<0, 1>(varrho, zero.data(), pn_py_o.data());
			DR::Run<2, 0>(varrho, zero.data(), pn_pxx_o.data());
			DR::Run<1, 1>(varrho, zero.data(), pn_pxy_o.data());
			DR::Run<0, 2>(varrho, zero.data(), pn_pyy_o.data());
			DR::Run<1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Run<0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Run<2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Run<1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Run<0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			pn_lap_o = pn_pxx_o + pn_pyy_o;
		}

	public:
		std::vector<MatPP> invMat;
		std::unordered_map<int, MatPP> invNeu;

		R varrho;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_px_o;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_py_o;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_pxx_o;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_pxy_o;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_pyy_o;
		EiRun::Matrix<R, 2, PN::value, EiRun::RowMajor>						pn_p_o;
		EiRun::Matrix<R, mMath::H<2, 2>::value, PN::value, EiRun::RowMajor>	pn_pp_o;
		EiRun::Matrix<R, 1, PN::value, EiRun::RowMajor>						pn_lap_o;
	};

	template <typename R, int P>
	class Particle_x<R,3,P> : public Particle<R,3,Particle_x<R,3,P>> {};

}