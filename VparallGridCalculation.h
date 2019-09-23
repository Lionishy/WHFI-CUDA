#pragma once
#ifndef VparallGridCalculation_H
#define VparallGridCalculation_H

#include "SimpleTable.h"
#include "PhysicalParameters.h"
#include "DispersionRelation.h"
#include "DispersionRelationResonantVelocitySolver.h"

#include <optional>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <cmath>
#include <future>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <iterator>

namespace iki { namespace whfi { 
	template <typename T, typename ZFunc_t>
	struct VparallGridCalculator {
		UniformSimpleTable<T, 1u, 6u> operator()(UniformAxis<T> const &v_parall_axis, std::vector<T> & v_parall_data) {
			auto v_parall_kernell = [v_parall_axis, this] (size_t counter)->std::array<T, 6u> {
				T v_resonant = v_parall_axis.begin + v_parall_axis.step * counter;
				auto k_omega_opt = rvsolver(v_resonant);
				if (!k_omega_opt) {
					std::stringstream error_text;
					error_text << "Can't find root for the resonant velocity " << v_resonant;
					throw std::runtime_error(error_text.str());
				}

				std::array<T, 6u> result_array;
				result_array[0] = k_omega_opt->first;
				result_array[1] = k_omega_opt->second;
				result_array[2] = omega_derive(k_omega_opt->second, k_omega_opt->first);
				result_array[3] = k_derive(k_omega_opt->second, k_omega_opt->first);
				result_array[4] = -result_array[3] / result_array[2];
				result_array[5] = T(1. / std::fabs(v_resonant - result_array[4]/params.betta_root_c));
				return result_array;
			};

			size_t counter = 0u; T k_prev = { 0 };
			while (true) {
				auto result_array = v_parall_kernell(counter);
				if (result_array[0] < k_prev) break;
				std::copy(std::begin(result_array), std::end(result_array), std::back_inserter(v_parall_data));
				++counter;
				k_prev = result_array[0];
			}

			UniformSimpleTable<T, 1u, 6u> v_parall_table;
			v_parall_table.space.axes[0] = v_parall_axis;
			v_parall_table.bounds.components[0] = counter;
			v_parall_table.data = v_parall_data.data();
			return v_parall_table;
		}

		VparallGridCalculator(ZFunc_t Z, PhysicalParamenters<T> params): params(params), rvsolver(Z,params), omega_derive(Z,params), k_derive(Z,params) { }
	private:
		PhysicalParamenters<T> params;
		ResonantVelocitySolver<T, ZFunc_t> rvsolver;
		DispersionRelationOmegaDerivative<T, ZFunc_t> omega_derive;
		DispersionRelationKDerivative<T, ZFunc_t> k_derive;
	};
} /* whfi */ } /* iki */

#endif