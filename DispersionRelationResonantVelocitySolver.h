#pragma once
#ifndef DispersionRelationResonantVelocitySolver_H
#define DispersionRelationResonantVelocitySolver_H

#include "PhysicalParameters.h"
#include "StepSolver.h"

#include <optional>
#include <utility>
#include <cmath>

namespace iki { namespace whfi {
	template <typename T, typename ZFunc_t>
	std::optional<std::pair<T, T>> DRResonantVelocitySolve(T resonant_v, PhysicalParamenters<T> const &params, ZFunc_t Z) {
		auto Zc = Z(resonant_v - params.bulk_to_term_c)
			, Zh = Z(resonant_v * std::sqrt(params.TcTh_ratio) - params.bulk_to_term_h);
		auto dispersion_relation = [&params,resonant_v,Zc,Zh] (T omega) {
			return T(1. / 1836.) 
				+ (omega - 1) * (omega - 1) / (resonant_v * params.betta_root_c) / (resonant_v * params.betta_root_c)
				- params.nc * ((omega * resonant_v) / (omega - T(1.)) - params.bulk_to_term_c) * Zc
				- params.nh * ((omega * resonant_v * std::sqrt(params.TcTh_ratio)) / (omega - T(1.)) - params.bulk_to_term_h) * Zh;
		};

		auto omega_opt = math::make_step_solver<T>(dispersion_relation)(T(0.), T(1.), T(1.e-4));
		if (!omega_opt) return std::nullopt;
		return std::optional<std::pair<T,T>>({ (*omega_opt - T(1.)) / (resonant_v * params.betta_root_c) ,*omega_opt });
	}
} /* whfi */ } /* iki */

#endif /* DispersionRelationResonantVelocitySolver_H */