#pragma once
#ifndef DispersionRelation_CUH
#define DispersionRelation_CUH

#include "PhysicalParameters.h"
#include <cuda.h>
#include <cuda_runtime.h>

namespace iki {	namespace whfi { namespace device {
	template <typename T, typename ZFunc_t>
	struct DispersionRelation {
		__device__ T operator()(T omega, T k) {
			return T(1. / 1836.) + k * k
				- p.nc * (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * Z((omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c)
				- p.nh * (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * Z((omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h);
		}

		__device__ DispersionRelation(ZFunc_t Z, PhysicalParamenters<T> p): Z(Z), p(p) { }

	private:
		ZFunc_t Z;
		PhysicalParamenters<T> p;
	};

	template <typename T, typename ZFunc_t>
	struct DispersionRelationDerivative {
		__device__ T operator()(T omega, T k) {
			T arg_c = (omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c;
			T arg_h = (omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h;
			T Zc = Z(arg_c), Zh = Z(arg_h);
			return p.nc / (k * p.betta_root_c) * (-Zc + (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * (Zc * arg_c + T(1.)))
				+ p.nh / (k * p.betta_root_h) * (-Zh + (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * (Zh * arg_h + T(1.)));
		}

		__device__ DispersionRelationDerivative(ZFunc_t Z, PhysicalParamenters<T> p) : Z(Z), p(p) { }

	private:
		ZFunc_t Z;
		PhysicalParamenters<T> p;
	};

} /* cuda */ } /* whfi */ } /* iki */

#endif