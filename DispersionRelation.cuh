#pragma once
#ifndef DispersionRelation_H
#define DispersionRelation_H

#include "PhysicalParameters.h"

#include <cuda.h>
#include <cuda_runtime.h>

namespace iki {	namespace whfi { namespace cuda {
	template <typename T>
	__device__ T dispersion_realation(PhysicalParamenters<T> params, T omega, T k) {
	}

} /* cuda */ } /* whfi */ } /* iki */

#endif