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

namespace SIM {

	template <typename R, int D, typename Derived>
	class Particle : public Base<R, D> {
	public:
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1> vec;
		typedef Eigen::Matrix<R, 1, 3> vec13;
		typedef Eigen::Matrix<R, 5, 1> vec5;
		typedef Eigen::Matrix<R, 6, 1> vec6;
		typedef Eigen::Matrix<R, 7, 1> vec7;
		typedef Eigen::Matrix<R, 8, 1> vec8;
		typedef Eigen::Matrix<R, 9, 1> vec9;
		typedef Eigen::Matrix<R, 5, 5> mat55;
		typedef Eigen::Matrix<R, 6, 6> mat66;
		typedef Eigen::Matrix<R, 7, 7> mat77;
		typedef Eigen::Matrix<R, 8, 8> mat88;
		typedef Eigen::Matrix<R, 9, 9> mat99;
		typedef Eigen::Matrix<R, 5, 3> mat53;
		typedef Eigen::Matrix<R, 6, 3> mat63;
		typedef Eigen::Matrix<R, 7, 3> mat73;
		typedef Eigen::Matrix<R, 8, 3> mat83;
		typedef Eigen::Matrix<R, 9, 3> mat93;
		typedef Eigen::Matrix<R, D, D> mat;
	public:
		Particle() {}
		~Particle() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void clean() {
			type.clear();
			for (int d = 0; d < D; d++) {
				pos[d].clear(); pos_m1[d].clear();
				vel[d].clear(); vel_p1[d].clear(); vel_m1[d].clear();
			}
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
				file << type[p] << " ";
				for (int d = 0; d < D; d++) {
					file << pos[d][p] << " ";
				}
				for (int d = 0; d < D; d++) {
					file << vel[d][p] << " ";
				}
				file << temp[p] << std::endl;
			}
			std::cout << " Writing Geo. done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			int n;	 int t;		vec p;		vec	v;	R tp;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File Geo. not found ! " << std::endl;
			file >> ct >> dp >> np >> bd1 >> bd2;
			n = np;
			while (n-- > 0) {
				file >> t;
				for (int d = 0; d < D; d++) file >> p[d];
				for (int d = 0; d < D; d++) file >> v[d];
				file >> tp;
				addPart(pType(t), p, v, tp);
			}
			file.close();
			std::cout << " Reading Geo. done " << std::endl;
		}

		void addPart(const pType& t, const vec& p, const vec& v) {
			type.push_back(t);
			for (int d = 0; d < D; d++) {
				pos[d].push_back(p[d]); pos_m1[d].push_back(p[d]);
				vel[d].push_back(v[d]); vel_p1.push_back(v[d]); vel_m1[d].push_back(v[d]);
			}
			pres.push_back(R(0));
			phi.push_back(R(0)); vort.push_back(R(0));
			norm.push_back(vec::Zero());
			bdc.push_back(0);
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (int p = 0; p < pos.size(); p++) {
				b += pos[p];
			}
			b.Expand(0.1);
			cell = new LinkCell<R, D>(b, r0);
			updateCell();
		}
		inline void updateCell() {
			cell->update(pos);
		}

		void b2b() {
			bbMap.clear();
			for (int p = 0; p < pos.size(); p++) {
				if (type[p] != BD2) continue;
				auto tmpdr = std::numeric_limits<R>::max();
				int tmpbb = 0;
				const auto c = cell->iCoord(pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (int m = 0; m < cell->linkList[key].size(); m++) {
						const int q = cell->linkList[key][m];
						if (q == p || type[q] != BD1) continue;
						const auto dr1 = (pos[q] - pos[p]).norm();
						if (dr1 < tmpdr) {
							tmpdr = dr1;
							tmpbb = q;
						}
					}
				}
				bbMap[p] = tmpbb;
			}
		}

		void b2norm() {
			bdnorm.clear();
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				const auto tmp = bdnorm.at(p.second);
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
		std::vector<R> pos[D];
		std::vector<R> pos_m1[D];
		std::vector<R> vel[D];
		std::vector<R> vel_p1[D];
		std::vector<R> vel_m1[D];

		std::vector<R> temp;
		std::vector<R> pnd;
		std::vector<R> pres;
		std::vector<pType> type;
		std::vector<int> bdc;
		std::vector<R> phi;
		std::vector<R> vort;
		std::unordered_map<int, vec> bdnorm;
		std::unordered_map<int, R> p_dirichlet;
		std::unordered_map<int, R> t_dirichlet;
		std::unordered_map<int, R> p_neumann;
		std::unordered_map<int, R> t_neumann;
		std::unordered_map<int, int> bbMap;

		LinkCell<R, D>* cell;

	private:
	};

}