#pragma once
#ifndef StepSolver_CUH
#define StepSolver_CUH

#include "DeviceOptional.h"
#include <cuda_runtime.h>

namespace iki { namespace math { namespace device {
	template <typename T, typename F_t>
	struct StepSolver {
		__device__ Optional<T> operator()(T min, T max, T step) const {
			size_t counter = 0u; T arg_curr = min, arg_next = min + step, f_curr = eqn(arg_curr), f_next = eqn(arg_next);
			while (arg_next < max) {
				if (f_curr * f_next < 0.) { Optional<T> opt; opt.is_present = true; opt.value.data = (arg_curr + arg_next) * T(0.5); return opt; }
				arg_curr = arg_next; f_curr = f_next;
				arg_next = min + ++counter * step; f_next = eqn(arg_next);
			}
			Optional<T> opt; opt.is_present = false; opt.value.stub = 0;
			return opt;
		}

		__device__ StepSolver(F_t equation) : eqn(equation) { }

	private:
		F_t eqn;
	};
} /* device */ } /* math */ } /* iki */

#endif