#pragma once
#ifndef StepSolver_H
#define StepSolver_H

#include <optional>

namespace iki { namespace math {
	template <typename T, typename F_t>
	struct StepSolver {
		std::optional<T> operator()(T min, T max, T step) const {
			size_t counter = 0u; T arg_curr = min, arg_next = min + step, f_curr = eqn(arg_curr), f_next = eqn(arg_next);
			while (arg_next < max) {
				if (f_curr * f_next < 0.) { return std::make_optional((arg_curr + arg_next) * T(0.5)); }
				arg_curr = arg_next; f_curr = f_next;
				arg_next = min + ++counter * step; f_next = eqn(arg_next);
			}
			return std::nullopt;
		}

		StepSolver(F_t equation): eqn(equation) { }

	private:
		F_t eqn;
	};

	template <typename T, typename F_t>
	auto make_step_solver(F_t eqn) { return StepSolver<T, F_t>(eqn); }
} /* math */ } /* iki */

#endif