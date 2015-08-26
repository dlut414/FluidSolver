/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_2d_h.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Mls_2d_h : public Simulator < real, dim, Mls_2d_h<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Matrix<real, 6, 1> vecp;
		typedef Eigen::Triplet<real> tpl;
	public:
		Mls_2d_h() {}
		~Mls_2d_h() {}

		void step() {
			calInvMat();
			shift();
			calCell(); //update cells
			//insertRand();
			calInvMat();
			convect1();
			//calPnd();
			//makeFs();
			makeDirichlet();
			makeMat();
			makeSource();
			//bvpSource();
			solvMat();
			convect2();
			calInvMat();
			//makeFs();
			calForVis();
			check();
			//calInvMat(); //for sensor
			//profileOut();
			//sensorOut();
			////bvpAvgError();
			//bvpMaxError();
			//lapMaxError();
			//gradMaxError();
			//surfCol();
			//pthOrderVelSpatialFilter();
		}

		void convect1() {
			/*standard*/
			/*
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] += para.dt * part->vel2[p];
			}
			*/
			//-------------------------------------------------------------------------//
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
		}

		void convect2() {
			/*standard*/
			/*
			std::vector<vec> dash(part->np);
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			dash[p] = -para.dt / para.rho * part->grad_mls_poly2d_d(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel1[p] = part->vel2[p] = part->vel2[p] + dash[p];
			part->pos[p] += para.dt * dash[p];
			}
			*/
			//------------------------------------------------------------------------//
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -(para.dt / para.rho) * part->grad(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
				part->vel1[p] = part->vel2[p];
			}
			//------------------------------------------------------------------------//
			/*fourth-order Ronge Kuta*/
			/*
			std::vector<vec> x0(part->np);
			std::vector<vec> k2(part->np);
			std::vector<vec> k3(part->np);
			std::vector<vec> k4(part->np);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				x0[p] = part->pos[p];
				part->pos[p] += 0.5* para.dt* part->vel1[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				k2[p] = part->vel2[p] + 0.5* para.dt* (-1. / para.rho)* part->grad(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = x0[p] + 0.5* para.dt* k2[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				k3[p] = part->vel2[p] + 0.5* para.dt* (-1. / para.rho)* part->grad(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = x0[p] + para.dt* k3[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				k4[p] = part->vel2[p] + para.dt* (-1. / para.rho)* part->grad(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] = x0[p] + (para.dt / 6.)* (part->vel1[p] + 2.*k2[p] + 2.*k3[p] + k4[p]);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -(para.dt / para.rho) * part->grad(part->pres, p);
				part->vel1[p] = part->vel2[p];
			}
			*/
			//------------------------------------------------------------------------//
		}

		void makeMat() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(tpl(p, p, 1.));
					coef.push_back(tpl(p, part->bbMap.at(p), -1.));
					continue;
				}
				if (part->isFs(p)) {
					coef.push_back(tpl(p, p, 1.));
					continue;
				}
				real pqsum = 0.;
				//std::vector<unsigned> used;
				const auto mm = part->invMat[p];
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				//if (abs(mm.determinant()) > part->eps_mat) {
					for (int k = -1; k <= 1; k++) {
						for (int j = -1; j <= 1; j++) {
							for (int i = -1; i <= 1; i++) {
								const iVec3 ne = c + iVec3(i, j, k);
								const unsigned key = part->cell->hash(ne);
								//for (unsigned us = 0; us < used.size(); us++) { if (key == used[us]) std::cout << "used!!!!!!!!!!" << std::endl; }
								//used.push_back(key);
								for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
									const unsigned q = part->cell->linkList[key][m];
#if BD_OPT
									if (part->bdOpt(p, q)) continue;
#endif
									const auto	dr = part->pos[q] - part->pos[p];
									const auto	dr1 = dr.mag();
									if (dr1 > part->r0) continue;
									const auto w = part->w3(dr1);
									const auto npq = part->poly(dr);
									const auto a = mm * (w* npq);
									const auto lp = part->poly_lap_0;
									real pq = lp.dot(a);
									if (part->type[p] == BD1) {
										//pq += a.dot(part->poly2d_d_px(vec(0.)) + part->poly2d_d_pz(vec(0.)));
									}
									coef.push_back(tpl(p, q, pq));
									pqsum += abs(pq);
								}
							}
						}
					}
					if (pqsum < para.eps) coef.push_back(tpl(p, p, 1.));
				//}
//				else {
//					auto pp = 0.;
//					for (int k = -1; k <= 1; k++) {
//						for (int j = -1; j <= 1; j++) {
//							for (int i = -1; i <= 1; i++) {
//								const iVec3 ne = c + iVec3(i, j, k);
//								const unsigned key = part->cell->hash(ne);
//								for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
//									const unsigned q = part->cell->linkList[key][m];
//#if BD_OPT
//									if (part->bdOpt(p, q)) continue;
//#endif
//									const auto	dr = part->pos[q] - part->pos[p];
//									const auto	dr1 = dr.mag();
//									if (dr1 > part->r0) continue;
//									const auto w = part->w3(dr1);
//									pp -= w;
//									coef.push_back(tpl(p, q, w));
//								}
//							}
//						}
//					}
//					coef.push_back(tpl(p, p, pp));
//				}
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeSource() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = 0.;
					continue;
				}
				const real coef = para.rho / para.dt;
				mSol->b[p] = coef * part->div(part->vel2, p);
			}
		}

		void init_() {
			part = new Particle_2d_h<real, dim>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->b2b();
			part->b2norm();
			//part->updateTeam();
			part->init2d_x();
		}

	public:
		Particle_2d_h<real, dim>*  part;

	private:
		std::vector<tpl> coef;
	};

}