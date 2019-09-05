#pragma once
#ifndef RungeKutta_H
#define RungeKutta_H

#include <cuda_runtime.h>

namespace iki { namespace math { namespace device { 
	template <typename T, typename F_t>
	struct Runge4th {
		__device__ T operator()(size_t x_cnt, T y) const {
			T x = x0 + x_cnt * dx;
			T k1 = F(x, y);
			T k2 = F(x + dx / T(2.), y + k1 * dx / T(2.));
			T k3 = F(x + dx / T(2.), y + k2 * dx / T(2.));
			T k4 = F(x + dx, y + k2 * dx);
			return dx / T(6.) * (k1 + T(2.) * k2 + T(2.) * k3 + k4);
		}

		__device__ Runge4th(T x0, T dx, F_t F): x0(x0), dx(dx), F(F) { }

	private:
		T const x0, dx;
		F_t F;
	};

	template <typename T, typename F_t>
	__device__ Runge4th<T, F_t> make_runge4th(T x0, T dx, F_t F) { return { x0,dx,F }; }
} /* device */ } /* math */ } /* iki */

#endif /* RungeKutta_H */
