/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
//Simulator.h
///defination of class Simulator
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "Header.h"
#include "Parameter.h"
#include "Particle.h"
#include "MatSolver.h"
#include "Sensor.h"

namespace SIM {

	template <typename R, int D, typename Derived>
	class Simulator {};

	template <typename R, typename Derived>
	class Simulator<R,1,Derived> {};

	template <typename R, typename Derived>
	class Simulator<R,2,Derived> {
		typedef Eigen::Matrix<R, 2, 1> Vec;
		typedef Eigen::Matrix<R, 2, 2> Mat;
		typedef Eigen::Triplet<R> Tpl;
	public:
		Simulator() { timeStep = 0; }
		~Simulator() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void operator >> (const std::string& str) const {
			saveData(str);
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " No file Para. found ! " << std::endl;
			file >> para.k >> para.Pr >> para.Ra >> para.cfl >> para.dtMax >> para.tt >> para.eps >> para.alpha >> para.beta;
			std::cout << " Effective radius (times of dp)   : " << para.k << std::endl;
			std::cout << " Prandtl number                   : " << para.Pr << std::endl;
			std::cout << " Rayleigh number                  : " << para.Ra << std::endl;
			std::cout << " CFL number                       : " << para.cfl << std::endl;
			std::cout << " Maximum time step (1)            : " << para.dtMax << std::endl;
			std::cout << " Total time (1)                   : " << para.tt << std::endl;
			std::cout << " EPS                              : " << para.eps << std::endl;
			std::cout << " Arbitrary parameter Alpha        : " << para.alpha << std::endl;
			std::cout << " Arbitrary parameter Beta         : " << para.beta << std::endl;
			std::cout << " Reading Para. done " << std::endl;
			file.close();
		}

		void init() {
			*this << "Para.txt";
			derived().init_();
			mSol = new MatSolver<R,2>(int(derived().part->np), para.eps);
			std::cout << " Particle number : " << derived().part->np << std::endl;
			sen << "Sensor.in";
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			timeStep = int(derived().part->ct / para.dt);
		}

		void mainLoop() {
			auto* const part = derived().part;
			while (part->ct <= para.tt) {
				std::cout << " step ----------------------------------> " << timeStep << std::endl;
				R tmp = cfl();
				para.dt = tmp < para.dtMax ? tmp : para.dtMax;
				part->updateCell();
				derived().step();
				part->ct += para.dt;	timeStep++;
				std::cout << " time --------> " << part->ct << std::endl;
				std::cout << " dt ----------> " << para.dt << std::endl;
			}
			saveData();
		}

		R stepGL() {
			auto* const part = derived().part;
			if (part->ct > para.tt) {
				saveData();
			}
			std::cout << " step ----------------------------------> " << timeStep << std::endl;
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			part->updateCell();
			derived().step();
			part->ct += para.dt;	timeStep++;
			std::cout << " time --------> " << part->ct << std::endl;
			std::cout << " dt ----------> " << para.dt << std::endl;
			return part->ct;
		}

		void sensorOut() {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			sen.writeVect(derived().part);
			sen >> convert.str();
		}

		void profileOut(const R& rt) {
			static std::string pf = "profile";
			sen.writeScal(derived().part);
			sen.profile(rt, pf);
		}

		void saveData() const {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().part) >> ("./out/" + convert.str() + ".out");
		}
		void saveData(const std::string& str) const {
			*(derived().part) >> ("./out/" + str + ".out");
		}

		__forceinline const R* positionX() const {
			return derived().part->pos[0].data();
		}
		__forceinline const R* positionY() const {
			return derived().part->pos[1].data();
		}
		__forceinline const R* scalar() const {
			return derived().part->phi.data();
		}
		__forceinline const int* type() const {
			return (int*)(derived().part->type.data());
		}

	public:
		Parameter<R, 2> para;
		MatSolver<R, 2>* mSol;
		Sensor<R, 2> sen;

	protected:
		void step() {}
		void convect() {}
		void visTerm_e() {}
		void visTerm_i() {}
		void presTerm_e() {}
		void presTerm_i() {}
		void makeDirchlet_v() {}

		void makeNeumann_p() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != BD1) continue;
				part->neumann[p] = R(0);
			}
		}

		void solvMat_p() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->pres[p] = mSol->x[p];
				if (part->pres[p] < -1.e5) part->pres[p] = -1.e5;
				if (part->pres[p] > 1.e5) part->pres[p] = 1.e5;
			}
		}

		void solvMat_phi() {
			auto* const part = derived().part;
			mSol->ccBiCg_augment(part->type);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->phi[p] = mSol->x[p];
				if (part->phi[p] < -1.e5) part->phi[p] = -1.e5;
				if (part->phi[p] > 1.e5) part->phi[p] = 1.e5;
			}
		}

		void solvMat_v() {
			mSol->biCg_v();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[0][p] = mSol->u[D*p];
				part->vel2[1][p] = mSol->u[D*p+1];
			}
		}

		R cfl() {
			R umax = 0.;
			const auto* const part = derived().part;
			for (int p = 0; p < part->np; p++) {
				const R ux = part->vel1[0][p];
				const R uy = part->vel1[1][p];
				const R tmp = sqrt(ux*ux + uy*uy);
				if (tmp > umax) umax = tmp;
			}
			para.umax = umax;
			return para.cfl * part->dp / umax;
		}

		void calCell() {
			derived().part->updateCell();
		}

		void calInvMat() {
			derived().part->updateInvMat();
		}

		void makeFs() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) part->fs[p] = part->_isFs(p);
		}

		void calForVis() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < part->np; p++) {
				part->vort[p] = part->rot(part->vel1, p);
			}
		}

		void check() const {
			const auto* const part = derived().part;
			R velMax = std::numeric_limits<R>::min();
			R phiMax = std::numeric_limits<R>::min();
			R divMax = std::numeric_limits<R>::min();
			int idv = 0, idp = 0, idd = 0;
			for (int p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) continue;
				const R ux = part->vel1[0][p];
				const R uy = part->vel1[1][p];
				const R vel = sqrt(ux*ux + uy*uy);
				const R phi = part->phi[p];
				const R div = part->div(part->vel1, p);
				if (vel > velMax) {
					velMax = vel;
					idv = p;
				}
				if (abs(phi) > abs(phiMax)) {
					phiMax = phi;
					idp = p;
				}
				if (abs(div) > abs(divMax)) {
					divMax = div;
					idd = p;
				}
			}
			std::cout << " max vel: " << velMax << " --- id: " << idv << std::endl;
			std::cout << " max phi: " << phiMax << " --- id: " << idp << std::endl;
			std::cout << " max Div: " << divMax << " --- id: " << idd << std::endl;
		}

		void insertRand() {
			auto* const part = derived().part;
			R coef = 0.25;
			std::default_random_engine gen;
			std::normal_distribution<R> dis(0., 0.5);
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const R dr = coef* part->dp* dis(gen);
				const R theta = 2.* M_PI * (R(rand()) / RAND_MAX);
				const R dx = cos(theta)*dr;
				const R dy = sin(theta)*dr;
				part->pos[0][p] += dx;
				part->pos[1][p] += dy;
				part->pos_m1[0][p] += dx;
				part->pos_m1[1][p] += dy;
			}
		}

	protected:
		int timeStep;
	};

	template <typename R, typename Derived>
	class Simulator<R,3,Derived> {};
}