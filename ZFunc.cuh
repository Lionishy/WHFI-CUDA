#pragma once
#ifndef ZFunc_H
#define ZFunc_H

#include <stddef.h>

#include <cuda.h>
#include <cuda_runtime.h>
	
namespace iki {	namespace math { namespace cuda{
	template <typename T>
	struct ZFuncTable {
		size_t size = 0u;
		T step = T(1);
		T *values = NULL;
	};

	template <typename T>
	__device__ T zfunc(ZFuncTable<T> const *table_ptr, T arg) {
	}
} /* cuda */ } /* math */ } /* iki */

#endif