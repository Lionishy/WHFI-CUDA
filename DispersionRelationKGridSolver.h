#pragma once

#include "DispersionRelation.h"
#include "StepSolver.h"

#include <vector>
#include <utility>
#include <optional>

namespace iki { namespace whfi {
	template <typename T, typename ZFunc_t>
	std::vector<std::pair<T,std::optional<T>>> DRKGridSolve(T k_min, T k_max, T k_step, DispersionRelation<T,ZFunc_t> const &dr) {
		auto step_solver_kernell = [&dr] (T k)->std::pair<T, std::optional<T>> { 
			return {
				k
				, StepSolver([&dr,k] (T omega) { return dr(omega,k); })(T(0.),T(1.),T(1.e-4)
			};  
		};
	}
} /* whfi */ } /* iki */