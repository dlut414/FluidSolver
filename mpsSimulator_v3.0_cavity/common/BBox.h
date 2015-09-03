/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#ifndef BBOX_H
#define BBOX_H

#include <limits>
#include <algorithm>
#include <Eigen/Dense>

template <typename R, int D = 2>
class BBox {
	typedef Eigen::Matrix<R, D, 1> vec;
public:
    BBox() {
        R _infinity = std::numeric_limits<R>::infinity();
		for (int i = 0; i < D; i++) {
			pMin[i] = _infinity;
			pMax[i] = -_infinity;
		}
    }

	template <typename T>
	BBox(const BBox<T, D>& b) {
		for (int i = 0; i < D; i++) {
			pMin[i] = static_cast<R>(b.pMin[i]);
			pMax[i] = static_cast<R>(b.pMax[i]);
		}
	}

    BBox(const vec& p) {
		for (int i = 0; i < D; i++) {
			pMin[i] = p[i];
			pMax[i] = p[i];
		}
	}

    BBox(const vec& p1, const vec& p2) {
		for (int i = 0; i < D; i++) {
			pMin[i] = std::min(p1[i], p2[i]);
			pMax[i] = std::max(p1[i], p2[i]);
		}
    }

    void operator+=(const vec& p) {
		for (int i = 0; i < D; i++) {
			pMin[i] = std::min(pMin[i], p[i]);
			pMax[i] = std::max(pMax[i], p[i]);
		}
    }

    BBox Union(const BBox& b1, const BBox &b2) const {
        BBox ret = b1;
        ret.pMin.x = std::min(b1.pMin.x, b2.pMin.x);
        ret.pMin.y = std::min(b1.pMin.y, b2.pMin.y);
        ret.pMin.z = std::min(b1.pMin.z, b2.pMin.z);
        ret.pMax.x = std::max(b1.pMax.x, b2.pMax.x);
        ret.pMax.y = std::max(b1.pMax.y, b2.pMax.y);
        ret.pMax.z = std::max(b1.pMax.z, b2.pMax.z);
        return ret;
    }

    bool bOverlaps(const BBox &b) const {
        bool x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
        bool y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
        bool z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
        return (x && y && z);
    }

    bool bInside(const vec &p) const {
        return (p.x >= pMin.x && p.x <=pMax.x &&
                p.y >= pMin.y && p.y <=pMax.y &&
                p.z >= pMin.z && p.z <=pMax.z);
    }

    vec pCenter() const {
        return vec( (pMax.x+pMin.x), (pMax.y+pMin.y), (pMax.z+pMin.z) ) * 0.5f;
    }

    void Expand(const real& s) {
        vec v = s * (pMax - pMin);
        pMax += v;
        pMin -= v;
    }

public:
    vec pMin, pMax;

};

#endif // BBOX_H
