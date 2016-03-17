/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "LinkCell.h"
#include "Header.h"
#include "Base.h"
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif

namespace SIM {

	template <typename R, int D, typename Derived>
	class Particle : public Base<R,D> {};

	template <typename R, typename Derived>
	class Particle<R,1,Derived> : public Base<R,1> {};

	template <typename R, typename Derived>
	class Particle<R,2,Derived> : public Base<R,2> {
	public:
		typedef Eigen::Matrix<int,2,1> iVec;
		typedef Eigen::Matrix<R,2,1> Vec;
	public:
		Particle() {}
		~Particle() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void clean() {
			type.clear();
			pos[0].clear(); pos[1].clear(); pos_m1[0].clear(); pos_m1[1].clear();
			vel[0].clear(); vel[1].clear(); vel_p1[0].clear(); vel_p1[1].clear(); vel_m1[0].clear(); vel_m1[1].clear();
			pres.clear();
			phi.clear(); vort.clear();
		}
		void operator >> (const std::string str) const {
			std::ofstream file(str, std::ofstream::out);
			file << std::scientific << std::setprecision(6) << ct << std::endl;
			file << std::scientific << std::setprecision(6) << dp << std::endl;
			file << np << " " << bd1 << " " << bd2 << std::endl;
			for (int p = 0; p < np; p++) {
				file << std::scientific << std::setprecision(6);
				file << type[p] << " " << pos[0][p] << " " << pos[1][p] << " " << vel[0][p] << " " << vel[1][p] << " " << temp[p] << std::endl;
			}
			std::cout << " Writing Geo. done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			int n;	 int t;		Vec p;		Vec	v;	R tp;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File Geo. not found ! " << std::endl;
			file >> ct >> dp >> np >> bd1 >> bd2;
			n = np;
			while (n-- > 0) {
				file >> t;
				file >> p[0] >> p[1];
				file >> v[0] >> v[1];
				file >> tp;
				addPart(pType(t), p, v, tp);
			}
			file.close();
			std::cout << " Reading Geo. done " << std::endl;
		}

		void addPart(const pType& t, const Vec& p, const Vec& v, const R& tp) {
			type.push_back(t);
			pos[0].push_back(p[0]);	pos[1].push_back(p[1]);
			pos_m1[0].push_back(p[0]); pos_m1[1].push_back(p[1]);
			vel[0].push_back(v[0]); vel[1].push_back(v[1]);
			vel_p1[0].push_back(v[0]); vel_p1[1].push_back(v[1]);
			vel_m1[0].push_back(v[0]); vel_m1[1].push_back(v[1]);
			temp.push_back(tp); pres.push_back(R(0)); phi.push_back(R(0)); vort.push_back(R(0));
			bdc.push_back(0);
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (int p = 0; p < np; p++) {
				Vec pp;
				pp[0] = pos[0][p];
				pp[1] = pos[1][p];
				b += pp;
			}
			b.Expand(0.1);
			cell = new LinkCell<R,2>(b, r0);
			updateCell();
		}
		void updateCell() {
			cell->update(pos[0].data(), pos[1].data(), np);
		}

		void b2b() {
			bbMap.clear();
			for (int p = 0; p < np; p++) {
				if (type[p] != BD2) continue;
				R tmpdr = std::numeric_limits<R>::max();
				int tmpbb = 0;
				for (int q = 0; q < np; q++) {
					if (q == p || type[q] != BD1) continue;
					const R dr[2] = { pos[0][q] - pos[0][p], pos[1][q] - pos[1][p] };
					const R dr1 = sqrt(dr[0] * dr[0] + dr[1] * dr[1]);
					if (dr1 < tmpdr) {
						tmpdr = dr1;
						tmpbb = q;
					}
				}
				bbMap[p] = tmpbb;
			}
		}

		void b2norm() {
			bdnorm.clear();
			for (const auto& p : bbMap) {
				Vec n;
				n[0] = pos[0][p.second] - pos[0][p.first];
				n[1] = pos[1][p.second] - pos[1][p.first];
				bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				Vec n;
				n[0] = pos[0][p.second] - pos[0][p.first];
				n[1] = pos[1][p.second] - pos[1][p.first];
				const Vec tmp = bdnorm.at(p.second);
				if (tmp.norm() < n.norm()) bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				bdnorm[p.second] = bdnorm.at(p.second).normalized();
			}
		}

		void b2neumann() {
			p_neumann.clear();
			t_neumann.clear();
			for (int p = 0; p < np; p++) {
				if (IS(bdc[p], P_NEUMANN)) p_neumann[p] = R(0);
				if (IS(bdc[p], T_NEUMANN)) t_neumann[p] = R(0);
			}
		}
		void b2dirichlet() {
			p_dirichlet.clear();
			t_dirichlet.clear();
			for (int p = 0; p < np; p++) {
				if (IS(bdc[p], P_DIRICHLET)) p_dirichlet[p] = R(0);
				if (IS(bdc[p], T_DIRICHLET0)) t_dirichlet[p] = R(0);
				if (IS(bdc[p], T_DIRICHLET1)) t_dirichlet[p] = R(1);
			}
		}

	public:
		R ct;
		int np, bd1, bd2;
		std::vector<R> pos[2];
		std::vector<R> pos_m1[2];
		std::vector<R> vel[2];
		std::vector<R> vel_p1[2];
		std::vector<R> vel_m1[2];

		std::vector<R> temp;
		std::vector<R> pnd;
		std::vector<R> pres;
		std::vector<pType> type;
		std::vector<int> bdc;
		std::vector<R> phi;
		std::vector<R> vort;
		std::unordered_map<int, Vec> bdnorm;
		std::unordered_map<int, R> p_dirichlet;
		std::unordered_map<int, R> t_dirichlet;
		std::unordered_map<int, R> p_neumann;
		std::unordered_map<int, R> t_neumann;
		std::unordered_map<int, int> bbMap;

		LinkCell<R,2>* cell;
	};

	template <typename R, typename Derived>
	class Particle<R,3,Derived> : public Base<R,3> {};

}