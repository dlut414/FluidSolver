/*
*/
#pragma once
#include <list>
#include "Header.h"
#include "Parameter.h"
#include "Particle.h"

namespace SIM {
	
	template <typename R, int D>
	class Shifter {};

	template <typename R>
	class Shifter<R,1> {};

	template <typename R>
	class Shifter<R,2> {
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,2,2> Mat;
	public:
		Shifter() {}
		~Shifter() {}

		template <typename Der>
		void SpringModel(Particle<R,2,Der>* const part, const Parameter<R,2>& para) const {
			std::vector<R> Dposx(part->np, R(0));
			std::vector<R> Dposy(part->np, R(0));
			std::vector<R> Du1x(part->np, R(0));
			std::vector<R> Du1y(part->np, R(0));
			std::vector<R> Du2x(part->np, R(0));
			std::vector<R> Du2y(part->np, R(0));
			std::vector<R> Du3x(part->np, R(0));
			std::vector<R> Du3y(part->np, R(0));
			std::vector<R> Dtemp(part->np, R(0));
			const R coef = para.umax* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				R Dpq[2] = { 0.0, 0.0 };
				const auto& cell = part->cell;
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
						const R coeff = w / dr1;
						Dpq[0] -= coeff * dr[0];
						Dpq[1] -= coeff * dr[1];
					}
				}
				Dposx[p] = coef* Dpq[0];
				Dposy[p] = coef* Dpq[1];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				dp[p] = part->pos_m1[p];
				tmp1[p] = part->derived().func_lsA_upwind(part->vel1, p, dp[p]);
				tmp2[p] = part->derived().func_lsA_upwind(part->vel2, p, dp[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = dp[p];
				part->vel1[p] = tmp1[p];
				part->vel2[p] = tmp2[p];
			}
		}

		template <typename T, typename U>
		void StaticUpwindModel(T* const part, U& phi) const {
			U tmp(part->np);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				tmp[p] = part->derived().func_lsA_upwind(phi, p, part->pos_m1[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				part->pos[p] = part->pos_m1[p];
				phi[p] = tmp[p];
			}
		}

		template <typename T, typename U>
		void StaticWENOModel(T* const part, U& phi) const {
			U tmp(part->np);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				tmp[p] = part->derived().interpolateWENO(part->phi, p, part->pos_m1[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				part->pos[p] = part->pos_m1[p];
				phi[p] = tmp[p];
			}
		}

		template <typename R, class Der>
		void shiftNearest(Particle<R, 2, Der>* const part, const Parameter<R, 2>& para) const {
			std::vector<Vec> dp(part->pos);
			std::vector<Vec> tmp(part->np, Vec(0));
			for (int loop = 0; loop < 1; loop++) {
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < int(part->np); p++) {
					if (part->type[p] != FLUID) continue;
					const auto inf = std::numeric_limits<R>::max();
					std::list<int> knn(8, 0);
					std::list<R> knd(8, inf);
					//auto gc = Vec(0., 0., 0.);
					const auto c = part->cell->iCoord(part->pos[p]);
					for (auto k = -1; k <= 1; k++) {
						for (auto j = -1; j <= 1; j++) {
							for (auto i = -1; i <= 1; i++) {
								const auto ne = c + iVec3(i, j, k);
								const auto key = part->cell->hash(ne);
								for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
									const auto q = part->cell->linkList[key][m];
									if (q == p) continue;
									const auto dr = part->pos[q] - part->pos[p];
									const auto dr2 = dr.mag2();
									if (dr2 > part->r0*part->r0) continue;
									std::list<int>::iterator it1 = knn.begin();
									std::list<R>::iterator it2 = knd.begin();
									for (; it1 != knn.end(); ++it1, ++it2) {
										if (dr2 < *it2) {
											knn.insert(it1, q);
											knn.pop_back();
											knd.insert(it2, dr2);
											knd.pop_back();
											break;
										}
									}
								}
							}
						}
					}
					auto dpq = Vec(0., 0., 0.);
					for (std::list<int>::iterator it = knn.begin(); it != knn.end(); ++it) {
						const auto dr = part->pos[*it] - part->pos[p];
						const auto dr1 = dr.mag();
						dpq -= w_spring(dr1, part->dp*1.5)* dr.norm();
					}
					//if (part->isFs(p)) {
					//	gc = gc.norm();
					//	dpq = dpq - 1.*(dpq*gc)*gc;
					//}
					part->pos[p] += para.umax*para.dt* dpq;
				}
			}
			/*explicit Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const auto tmp = part->pos[p];
				part->pos[p] = dp[p];
				dp[p] = tmp;
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				if (part->isFs(p)) continue;
				tmp[p] = part->derived().func_mls_a_upwind(part->vel2, p, dp[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = dp[p];
				part->vel2[p] = tmp[p];
			}
		}

		template <typename R, , class Der>
		void shiftPnd(Particle<R, 2, Der>* part, Parameter<R, 2>& para) const {
			std::vector<Vec> dn(part->np, Vec(0.));
			for (int p = 0; p < part->pnd.size(); p++) {
				dn[p] = part->grad_suzuki_pnd(part->pnd, p);
			}
			R nn = 0.;
			for (int p = 0; p < dn.size(); p++) {
				nn += dn[p] * dn[p];
			}
			for (int p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				Vec gc = Vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const int key = part->cell->hash(ne);
							for (int m = 0; m < part->cell->linkList[key].size(); m++) {
								const int q = part->cell->linkList[key][m];
								if (q == p) continue;
								const Vec	dr = part->pos[p] - part->pos[q];
								const R  dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								gc += ( /*part->w1(dr1)*/1. / dr1)* dr;
							}
						}
					}
				}
				const R corr = part->nbd[p] > 0 ? R(part->pn0) / R(part->pn[p]) : 1.;
				Vec dp = (part->n0 - part->pnd[p] * corr) / nn * dn[p];
				if ((part->type[p] == FLUID) && part->isFs(p)) {
					gc = gc.norm();
					//dp = dp - (dp*gc)*gc;
				}
				dp = 20.*dp;
				part->vel1[p] += part->grad_suzuki(part->vel1, p) * dp;
				part->pos[p] += dp;
			}
		}

		template <typename R, class Der>
		void shiftPbf(Particle<R, 2, Der>* part, Parameter<R, 2>& para) const {
			std::vector<Vec> dp(part->pos.size());
			for (int p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				Vec dpq = Vec(0., 0., 0.);
				Vec gc = Vec(0., 0., 0.);
				R lam = 0.;
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const int key = part->cell->hash(ne);
							for (int m = 0; m < part->cell->linkList[key].size(); m++) {
								const int q = part->cell->linkList[key][m];
								if (q == p) continue;
								const Vec	dr = part->pos[q] - part->pos[p];
								const R	dr2 = dr.mag2();
								const R  dr1 = sqrt(dr2);
								if (dr1 > part->r0) continue;
								const Vec cp = (-part->r0 / (dr2 * dr1)) * dr;
								lam += cp*cp;
								dpq -= cp;

								gc -= part->w1(dr1) * dr;
							}
						}
					}
				}
				lam += dpq* dpq;
				lam = (part->n0 - part->pnd[p]) / (lam + 1.e-3);
				if ((part->type[p] == FLUID) && (part->pn[p] < part->pn0 * para.beta)) {
					gc = gc.norm();
					dpq = dpq - (dpq*gc)*gc;
				}
				dp[p] = 0.001*lam * dpq;
			}
			for (int p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += part->grad(part->vel1, p) * dp[p];
			}
			for (int p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += dp[p];
				part->vel1[p] = part->vel2[p];
			}
		}

		template <typename R, class Der>
		void shiftXu(Particle<R, 2, Der>* const part, const Parameter<R, 2>& para) const {
			std::vector<Vec> dp(part->pos.size());
			for (int p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				R pa = 0.;
				int count = 0;
				Vec dpq = Vec(0., 0., 0.);
				Vec gc = Vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const int key = part->cell->hash(ne);
							for (int m = 0; m < part->cell->linkList[key].size(); m++) {
								const int q = part->cell->linkList[key][m];
								if (q == p) continue;
								const Vec	dr = part->pos[q] - part->pos[p];
								const R	dr2 = dr.mag2();
								const R  dr1 = sqrt(dr2);
								if (dr1 > part->r0) continue;
								pa += dr1;
								count++;
								dpq -= (1. / dr2) * (dr / dr1);
								gc -= /*w2(dr1, part->r0)**/ (dr / dr1);
							}
						}
					}
				}
				count = (count == 0) ? 1 : count;
				pa = pa / count;

				if (part->isFs(p)) {
					gc = gc.norm();
					dpq = dpq - 1.*(dpq*gc)*gc;
				}
				dp[p] = 0.05* para.umax * para.dt *pa*pa * dpq;
				//dp[p] = 0.05* para.cfl*para.dp *pa*pa * dpq;
			}
			for (int p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel2[p] += dp[p] * (part->derived().grad(part->vel1, p));
				//if ((part->vel2[p]).mag() > 5.) part->vel2[p] = 5. * (part->vel2[p]).norm();
			}
			for (int p = 0; p < part->pos.size(); p++) {
				if ((part->type[p] != FLUID)) continue;
				//part->vel1[p] = part->vel2[p];
				part->pos[p] += dp[p];
			}
		}

	private:
		inline const R w_spline(const R& r, const R& r0) const {
			/*cubic spline*/
			const R q = r / r0;
			if (q <= 0.5) return (1. - 6 * q*q + 6 * q*q*q);
			else if (q <= 1.) return (2.*pow((1. - q), 3));
			else return 0.;
			/*inverse r2*/
			//return pow(r0 / r, 3);
		}
		inline const R w2(const R& r, const R& r0) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				//return r0 / r - 1.;
				return pow((1 - r / r0), 2);
			}
		}
		inline const R w3(const R& r, const R& r0) const {
			/*parabolic*/
			const R q = r / r0;
			if (q < 1.) return 2.5*(1. - q*q);
			else return 0.;
		}
		inline const R w_spring(const R&r, const R& r0) const {
			if (r >= r0) return 0.;
			return 1. - r / r0;
		}
	};

	template <typename R>
	class Shifter<R,3> {};

}