/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//FractionalStep_KM_DC_A.h
///defination of class FractionalStep_KM_DC_A (Kim & Moin - using deferred correction)
///using Deferred Correction (DC) in the source term of pressure poisson equation
#pragma once
#include "Simulator.h"
#include "Particle_x.h"
#include "Shifter.h"

namespace SIM {
	
	template <typename R, int D, int P>
	class FractionalStep_KM_DC_A : public Simulator<R,D,FractionalStep_KM_DC_A<R,D,P>> {};

	template <typename R, int P>
	class FractionalStep_KM_DC_A<R,1,P> : public Simulator<R,1,FractionalStep_KM_DC_A<R,1,P>> {};

	template <typename R, int P>
	class FractionalStep_KM_DC_A<R,2,P> : public Simulator<R,2,FractionalStep_KM_DC_A<R,2,P>> {
		typedef mMath::Polynomial_A<R,2,P> PN;
		typedef mMath::Derivative_A<R,2,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<R,2,1> Vec;
		typedef Eigen::Matrix<R,PN::value,PN::value> MatPP;
		typedef Eigen::Triplet<R> Tpl;
	public:
		FractionalStep_KM_DC_A() {}
		~FractionalStep_KM_DC_A() {}

		void init_() {
			part = new Particle_x<R,2,P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k);
			part->buildCell();
			part->makeBdc();
			part->b2b();
			part->b2normal();
			part->b2neumann();
			part->b2dirichlet();
			part->init_x();
			sen = new Sensor<R,2,Particle_x<R,2,P>>(part);
			*sen << "Sensor.in";
			DGP_old.resize(part->np);
		}

		void step() {
			calInvMat();

			visTerm_i_q2r0();
			presTerm_i_q2();
			temperatureTerm_i_q1();

			syncPos();
			updateVelocity_q2();
			updatePosition_s2();

			calCell();
			calInvMat();
			calForVis();
			check();

			Redistribute();

			sync();
		}

		void Redistribute() {
			shi.StaticUpwindModel(part);
			//shi.StaticWENOModel(part);
			//shi.SpringUpwindModel(part, para);
		}

		void visTerm_i_q2r1() {
			makeLhs_v_q2();
			makeRhs_v_q2r1();
			solvMat_v();
		}

		void visTerm_i_q1r0() {
			makeLhs_v_q1();
			makeRhs_v_q1r0();
			solvMat_v();
		}

		void visTerm_i_q2r0() {
			makeLhs_v_q2();
			makeRhs_v_q2r0();
			solvMat_v();
		}

		void presTerm_i_q2() {
			for (int i = 0; i < 5; i++) {
				makeLhs_p();
				makeRhs_p_q2();
				solvMat_phi();
			}
			//for (int p = 0; p < part->np; p++) {
			//	//std::cout << (R(3.0) / (R(2.0)* para.dt))* part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p) - part->Lap(part->phi.data(), p) - mSol->x[part->np] << std::endl;
			//	//std::cout << (mSol->a.row(p).dot( mSol->x)) - mSol->b[p] << std::endl;
			//}
		}

		void presTerm_i_q1() {
			makeLhs_p();
			makeRhs_p_q1();
			solvMat_phi();
		}

		void temperatureTerm_i_q1() {
			makeLhs_t();
			makeRhs_t_q1();
			solvMat_t();
		}

		void updateVelocity_q1() {
			const R coefL = para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					part->pres[p] = part->phi[p];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec du = -coefL * part->Grad(part->phi.data(), p);
					part->vel_p1[0][p] += du[0];
					part->vel_p1[1][p] += du[1];
				}
			}
		}

		void updateVelocity_q2() {
			const R coefL = (2.0* para.dt) / (3.0);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					part->pres[p] = part->phi[p];
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					const Vec du = -coefL * part->Grad(part->phi.data(), p);
					part->vel_p1[0][p] += du[0];
					part->vel_p1[1][p] += du[1];
				}
				else if (part->type[p] == BD1) {
					part->vel_p1[0][p] = part->vel[0][p];
					part->vel_p1[1][p] = part->vel[1][p];
				}
			}
		}

		void updatePosition_s1() {
			const R coefL = R(0.5)* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += coefL * (part->vel[0][p] + part->vel_p1[0][p]);
					part->pos[1][p] += coefL * (part->vel[1][p] + part->vel_p1[1][p]);
				}
			}
		}

		void updatePosition_s2() {
			const R coefL = 0.5* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == FLUID) {
					part->pos[0][p] += coefL * (3.0* part->vel[0][p] - 1.0* part->vel_m1[0][p]);
					part->pos[1][p] += coefL * (3.0* part->vel[1][p] - 1.0* part->vel_m1[1][p]);
				}
			}
		}

	public:
		Particle_x<R,2,P>* part;
		Sensor<R,2,Particle_x<R,2,P>>* sen;

	private:
		void makeLhs_v_q2() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					coef.push_back(Tpl(2 * p, 2 * p, R(1)));
					coef.push_back(Tpl(2 * p + 1, 2 * p + 1, R(1)));
					continue;
				}
				R pp = R(0);
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = mm * (w* npq);
						const R pq = -(part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(2 * p, 2 * q, pq));
						coef.push_back(Tpl(2 * p + 1, 2 * q + 1, pq));
					}
				}
				pp += 3.0 / (2.0 * para.dt * para.Pr);
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void makeLhs_v_q1() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					coef.push_back(Tpl(2 * p, 2 * p, 1.0));
					coef.push_back(Tpl(2 * p + 1, 2 * p + 1, 1.0));
					continue;
				}
				R pp = R(0);
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = mm * (w* npq);
						const R pq = -(part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(2 * p, 2 * q, pq));
						coef.push_back(Tpl(2 * p + 1, 2 * q + 1, pq));
					}
				}
				pp += 1.0 / (para.dt * para.Pr);
				coef.push_back(Tpl(2 * p, 2 * p, pp));
				coef.push_back(Tpl(2 * p + 1, 2 * p + 1, pp));
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_v_q2r1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = part->vel[0][p];
					mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const Vec Gp = part->Grad(part->pres, p);
				const R coefL = 1.0 / (2.0* para.dt * para.Pr);
				const R rhsx = coefL* (4.0* part->vel[0][p] - part->vel_m1[0][p]) - (1.0 / para.Pr)* Gp[0];
				const R rhsy = coefL* (4.0* part->vel[1][p] - part->vel_m1[1][p]) - (1.0 / para.Pr)* Gp[1] + para.Ra* part->temp[p];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}

		void makeRhs_v_q1r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					mSol->rhs[2 * p + 0] = part->vel[0][p];
					mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const R coefL = (R(1.0) / para.dt * para.Pr);
				const R rhsx = coefL* part->vel[0][p];
				const R rhsy = coefL* part->vel[0][p] + para.Ra* part->temp[p];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}

		void makeRhs_v_q2r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					const R coef_local = R(2) / R(3);
					const Vec lap_local = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
					const R volumeForce_local = para.Ra* para.Pr* part->temp[p];
					mSol->rhs[2 * p + 0] = coef_local* (para.dt*(para.Pr* lap_local[0]) + R(2)*part->vel[0][p] - R(0.5)*part->vel_m1[0][p]);
					mSol->rhs[2 * p + 1] = coef_local* (para.dt*(para.Pr* lap_local[1] + volumeForce_local) + R(2)*part->vel[1][p] - R(0.5)*part->vel_m1[1][p]);
					//mSol->rhs[2 * p + 0] = part->vel[0][p];
					//mSol->rhs[2 * p + 1] = part->vel[1][p];
					continue;
				}
				const R coefL = R(1) / (R(2)* para.dt * para.Pr);
				const R rhsx = coefL* (R(4)* part->vel[0][p] - part->vel_m1[0][p]);
				const R rhsy = coefL* (R(4)* part->vel[1][p] - part->vel_m1[1][p]) + para.Ra* part->temp[p];
				mSol->rhs[2 * p + 0] = rhsx;
				mSol->rhs[2 * p + 1] = rhsy;
			}
		}


		void makeLhs_p() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, R(1)));
					continue;
				}
				R pqsum = R(0);
				R pp = R(0);
				MatPP* mm;
				if (IS(part->bdc[p], P_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = (*mm) * (w* npq);
						const R pq = part->pn_lap_o* aa;
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
						pqsum += abs(pq);
					}
				}
				coef.push_back(Tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, R(1.0)));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_p_q2() {
			const R coefL = R(3) / (R(2)* para.dt);
			part->DivGrad(part->phi.data(), DGP_old.data());
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				const R div_local = part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
				const R LP_old_local = part->Lap(part->phi.data(), p);
				mSol->b[p] = coefL * div_local +(LP_old_local - DGP_old[p]);
				if (IS(part->bdc[p], P_NEUMANN)) {
					Vec& normal = part->bdnorm.at(p);
					VecP inner = VecP::Zero();
					inner.block<2,1>(0, 0) = normal;
					const Vec lap_ustar_local = part->Lap(part->vel[0].data(), part->vel[1].data(), p);
					const R neumannX = para.Pr* lap_ustar_local[0];
					const R neumannY = para.Pr* lap_ustar_local[1] + para.Ra* para.Pr* part->temp[p];
					const R neumann = neumannX* normal[0] + neumannY* normal[1];
					const VecP aa = part->invNeu.at(p)* inner;
					//const R cst = part->p_neumann.at(p)*part->ww(R(0))* (R(1) / part->varrho) * (part->pn_lap_o.dot(aa));
					const R cst = neumann *part->ww(R(0))* (R(1) / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void makeRhs_p_q1() {
			const R coefL = R(1.0) / para.dt;
			part->DivGrad(part->phi.data(), DGP_old.data());
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.0;
					continue;
				}
				const R div_local = part->Div(part->vel_p1[0].data(), part->vel_p1[1].data(), p);
				const R LP_old_local = part->Lap(part->phi.data(), p);
				mSol->b[p] = coefL * div_local + (LP_old_local - DGP_old[p]);
				if (IS(part->bdc[p], P_NEUMANN)) {
					VecP inner = VecP::Zero();
					inner.block<2, 1>(0, 0) = part->bdnorm.at(p);
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = part->p_neumann.at(p)*part->ww(0.0)* (1.0 / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}

		void makeLhs_t() {
			coef.clear();
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, 1.0));
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					coef.push_back(Tpl(p, p, 1.0));
					continue;
				}
				R pqsum = R(0.0);
				R pp = R(0.0);
				MatPP* mm;
				if (IS(part->bdc[p], T_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
				const R coefL = -R(0.5)* para.dt;
				const auto& cell = part->cell;
				const int cx = cell->pos2cell(part->pos[0][p]);
				const int cy = cell->pos2cell(part->pos[1][p]);
				for (int i = 0; i < cell->blockSize::value; i++) {
					const int key = cell->hash(cx, cy, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (part->type[q] == BD2) continue;
						const R dr[2] = { part->pos[0][q] - part->pos[0][p], part->pos[1][q] - part->pos[1][p] };
						const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
						if (dr1 > part->r0) continue;
						const R w = part->ww(dr1);
						VecP npq;
						part->poly(dr, npq.data());
						const VecP aa = (*mm) * (w* npq);
						const R pq = coefL* (part->pn_lap_o.dot(aa));
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
					}
				}
				pp += 1.0;
				coef.push_back(Tpl(p, p, pp));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_t_q1() {
			const R coefL = 0.5* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.0;
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					mSol->b[p] = part->t_dirichlet.at(p);
					continue;
				}
				mSol->b[p] = part->temp[p] + coefL* part->Lap(part->temp.data(), p);
				if (IS(part->bdc[p], T_NEUMANN)) {
					VecP inner = VecP::Zero();
					inner.block<2,1>(0, 0) = part->bdnorm.at(p);
					const VecP aa = part->invNeu.at(p)* inner;
					const R cst = part->p_neumann.at(p)*part->ww(0.0)* (1.0 / part->varrho) * (part->pn_lap_o.dot(aa));
					mSol->b[p] -= cst;
				}
			}
		}


		void sync() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->vel_m1[0][p] = part->vel[0][p];
				part->vel_m1[1][p] = part->vel[1][p];
				part->vel[0][p] = part->vel_p1[0][p];
				part->vel[1][p] = part->vel_p1[1][p];
			}
		}
		void syncPos() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->pos_m1[0][p] = part->pos[0][p];
				part->pos_m1[1][p] = part->pos[1][p];
			}
		}

	private:
		Shifter<R,2> shi;
		std::vector<Tpl> coef;
		std::vector<R> DGP_old;
	};

	template <typename R, int P>
	class FractionalStep_KM_DC_A<R,3,P> : public Simulator<R,3,FractionalStep_KM_DC_A<R,3,P>>  {};
}