#pragma once
#ifndef ZFunc_H
#define ZFunc_H

#include "SimpleTable.h"

#include <stddef.h>

#include <cuda_runtime.h>
	
namespace iki {	namespace whfi { namespace device {
	template <typename T>
	struct ZFuncODE {
		__device__ T operator()(T x, T y) const { return -x * y - T(1.); }
	};

	template <typename T>
	struct ZFunc {
		__device__ T operator()(T arg) {
			T farg = fabs(arg);
			auto idx = size_t(farg / zfunc_table.space.axes[0].step);
			if ((idx + 1u) < collapsed_size(&table.bounds)) {
				return (arg > T(0.) ? T(1.) : T(-1.)) * ((zfunc_table.data[idx + 1u] - zfunc_table.data[idx]) / zfunc_table.space.axes[0].step * (farg - zfunc_table.space.axes[0].step * idx) + zfunc_table.data[idx]);
			}
			else { //asymptotic
				T over = T(1.) / arg, square = over * over;
				return -over * (T(1.) + square + T(3.) * square * square);
			}
		}

	private:
		UniformSimpleTable<T, 1u, 1u> zfunc_table;
	};
} /* cuda */ } /* math */ } /* iki */

#endif