#pragma once
#ifndef ZFunc_CH
#define ZFunc_CH

#include "SimpleTable.h"
#include <stddef.h>

namespace iki {	namespace whfi {
	template <typename T>
	struct ZFuncODE {
		T operator()(T x, T y) const { return -x * y - T(1.); }
	};

	template <typename T>
	struct ZFunc {
		T operator()(T arg) {
			T farg = fabs(arg);
			auto idx = size_t(farg / zfunc_table.space.axes[0].step);
			if ((idx + 1u) < collapsed_size(&zfunc_table.bounds)) {
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
} /* math */ } /* iki */

#endif