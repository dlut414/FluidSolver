/*
*/
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include "Header.h"

#define AUGMENT (1)
#define AG AUGMENT

namespace SIM {

	template <typename R, int D>
	class MatSolver {
		typedef Eigen::Triplet<R> Tpl;
		typedef Eigen::Matrix<R, Eigen::Dynamic, 1> dVec;
		typedef Eigen::SparseMatrix<R, Eigen::RowMajor> sMat;
#if AUGMENT
		typedef Eigen::IncompleteLUT<R> preconditioner;
#else
		typedef Eigen::DiagonalPreconditioner<R> preconditioner;
#endif
	public:
		MatSolver(const int& _n, const R& e)
			: n(_n), Dn(D*_n), eps(e), 
			a(_n + AG, _n + AG), x(_n + AG), b(_n + AG), 
			au(D*_n, D*_n), u(D*_n), rhs(D*_n), 
			vr(_n+AG), vr_hat(_n+AG), vv(_n+AG), vp(_n+AG), vt(_n+AG), vh(_n+AG), vs(_n+AG),
			Dvr(D*_n), Dvr_hat(D*_n), Dvv(D*_n), Dvp(D*_n), Dvt(D*_n), Dvh(D*_n), Dvs(D*_n) {
			init();
		}
		~MatSolver() {}

		void biCg() {
			//solverBiCg.compute(a);
			////x = solver1.solveWithGuess(b, 0.5*x);
			//x = solverBiCg.solve(b);
			//std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			//std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
			R rho;
			R alpha;
			R beta;
			R omega;
			R residual;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < n; p++) {
				const R ppInv = R(1.0) / a(p, p);
				for (int q = 0; q < n; q++) {
					a(p, q) *= ppInv;
				}
				b[p] *= ppInv;
				x[p] = R(0.0);
			}
			vr = b - a*x;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < n; p++) {
				vr_hat[p] = vr[p];
				vv[p] = vp[p] = R(0.0);
			}
			rho = alpha = omega = R(1.0);
			int loop = 0;
			for (loop = 0; loop < maxIter; loop++) {
				const R rho_m1 = rho;
				rho = vr_hat.dot(vr);
				beta = (rho / rho_m1)*(alpha / omega);
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < n; p++) {
					vp[p] = vr[p] + beta*(vp[p] - omega* vv[p]);
				}
				vv = a* vp;
				alpha = rho / (vr_hat.dot(vv));
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < n; p++) {
					vh[p] = x[p] + alpha*vp[p];
					vs[p] = vr[p] - alpha* vv[p];
				}
				vt = a*vs;
				omega = (vt.dot(vs)) / (vt.dot(vt));
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < n; p++) {
					x[p] = vh[p] + omega*vs[p];
				}
				vr = b - a*x;
				residual = vr.dot(vr);
				if (residual < eps*eps) break;
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < n; p++) {
					vr[p] = vs[p] - omega*vt[p];
				}
			}

			std::cout << " iterations ----------> " << loop+1 << std::endl;
			std::cout << " error ---------------> " << sqrt(residual) << std::endl;
		}
		void biCg_v() {
			//solverBiCg.compute(au);
			////u = solver2.solveWithGuess(rhs, 0.5*u);
			//u = solverBiCg.solve(rhs);
			//std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			//std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
			R rho;
			R alpha;
			R beta;
			R omega;
			R residual;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < Dn; p++) {
				const R ppInv = R(1.0) / a(p, p);
				for (int q = 0; q < Dn; q++) {
					au(p, q) *= ppInv;
				}
				rhs[p] *= ppInv;
				u[p] = R(0.0);
			}
			vr = b - a*x;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < Dn; p++) {
				Dvr_hat[p] = Dvr[p];
				Dvv[p] = Dvp[p] = R(0.0);
			}
			rho = alpha = omega = R(1.0);
			int loop = 0;
			for (loop = 0; loop < maxIter; loop++) {
				const R rho_m1 = rho;
				rho = Dvr_hat.dot(Dvr);
				beta = (rho / rho_m1)*(alpha / omega);
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < Dn; p++) {
					Dvp[p] = Dvr[p] + beta*(Dvp[p] - omega* Dvv[p]);
				}
				Dvv = au* Dvp;
				alpha = rho / (Dvr_hat.dot(Dvv));
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < Dn; p++) {
					Dvh[p] = u[p] + alpha*Dvp[p];
					Dvs[p] = Dvr[p] - alpha* Dvv[p];
				}
				Dvt = au*Dvs;
				omega = (Dvt.dot(Dvs)) / (Dvt.dot(Dvt));
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < Dn; p++) {
					u[p] = Dvh[p] + omega*Dvs[p];
				}
				Dvr = rhs - au*u;
				residual = Dvr.dot(Dvr);
				if (residual < eps*eps) break;
#if OMP
#pragma omp parallel for
#endif
				for (int p = 0; p < Dn; p++) {
					Dvr[p] = Dvs[p] - omega*Dvt[p];
				}
			}

			std::cout << " iterations ----------> " << loop + 1 << std::endl;
			std::cout << " error ---------------> " << sqrt(residual) << std::endl;
		}
		void qr() {
			a.makeCompressed();
			solverQR.compute(a);
			x = solverQR.solve(b);
			std::cout << " rank ----------------> " << solverQR.rank() << std::endl;
			if (solverQR.info() != Eigen::Success) {
				std::cout << " info ----------------> " << solverQR.lastErrorMessage() << std::endl;
			}
		}
		void lnBiCg() {
			sMat at(n, n);
			at = a.transpose();
			a = (a*at);
			solverBiCg.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCg.solve(b);
			x = (at*x);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void ccBiCg_augment(const std::vector<enum pType>& type) {
			sMat d(n + AG, n + AG);
			std::vector<Tpl> coef;
			for (int p = 0; p<n; p++) {
				if (type[p] == BD2) continue;
				coef.push_back(Tpl(n, p, 1.));
				coef.push_back(Tpl(p, n, 1.));
			}
			d.setFromTriplets(coef.begin(), coef.end());
			a = a + d;
			b[n] = 0.;
			solverBiCg.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCg.solve(b);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void ccCgs_horibata() {
			sMat at(n, n);
			at = a.transpose();
			dVec e(n), zero(n);
			solverBiCg.compute(at);
			e = solverBiCg.solve(0.*zero);
			const auto e2 = e.dot(e);
			if (e2 > eps) {
				b = b - (b.dot(e) / e2)*e;
			}
			const dVec r0 = b - a*x;
			auto r1 = r0;
			auto r2 = r1;
			auto p = r0;
			auto e1 = b;
			auto res = r1.dot(r1);
			auto h = b;
			int k;
			for (k = 0; k < maxIter; k++) {
				res = r1.dot(r1);
				if (res < eps) break;
				const auto alpha = res / r0.dot(a*p);
				h = e1 - alpha*(a*p);
				x += alpha*(e1 + h);
				r2 = r1 - alpha*(e1 + h);
				const auto beta = r0.dot(r2) / r0.dot(r1);
				e1 = r2 + beta* h;
				p = e1 + beta*(h + beta*p);
				r1 = r2;
			}
			std::cout << " iterations ----------> " << k << std::endl;
			std::cout << " error ---------------> " << res << std::endl;
		}
		void lsqr() {}

	public:
		int n;
		int Dn;
		int maxIter;
		R eps;
		sMat a, au;
		dVec x, b, u, rhs;
		Eigen::BiCGSTAB< sMat, preconditioner > solverBiCg;
		Eigen::SparseQR< sMat, Eigen::NaturalOrdering<int> > solverQR;

	private:
		void init() {
			maxIter = 1000;
			for (int i = 0; i < n; i++) {
				x[i] = b[i] = 0.;
			}
			for (int i = 0; i < D*n; i++) {
				u[i] = rhs[i] = 0.;
			}
#if AUGMENT
			solverBiCg.preconditioner().setDroptol(eps);
			solverBiCg.preconditioner().setFillfactor(1);
#endif
			solverBiCg.setMaxIterations(maxIter);
			solverBiCg.setTolerance(eps);
			solverQR.setPivotThreshold(1.0 / n);
		}
		void fina() {}

	private:
		dVec vr, vr_hat, vv, vp, vt, vh, vs;
		dVec Dvr, Dvr_hat, Dvv, Dvp, Dvt, Dvh, Dvs;
	};

}