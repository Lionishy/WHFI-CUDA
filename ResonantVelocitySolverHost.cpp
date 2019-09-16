
#include "SimpleTable.h"
#include "SimpleTableIO.h"
#include "ZFunc.h"
#include "ZFuncIO.h"
#include "PhysicalParameters.h"
#include "DispersionRelationResonantVelocitySolver.h"

#include <vector>
#include <fstream>
#include <iomanip>
#include <tuple>
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>
#include <algorithm>

void ResonantVelocitySolverHost() {
	using namespace std;
	using namespace iki;

	/* Load ZFunc table */
	whfi::ZFunc<float> zfunc;
	vector<float> zfunc_data;
	{
		std::ifstream binary_is("./fZFunc.tbl", ios::binary);
		binary_is.exceptions(ios::badbit | ios::failbit);
		whfi::ZFuncImport(binary_is, zfunc.zfunc_table, zfunc_data);
	}

	/* Set physical parameters */
	whfi::PhysicalParamenters<float> params = whfi::init_parameters(0.85f, 1.f / 0.85f, 0.25f, -9.f);
	
	try {
		vector<float> v_resonant_data;
		UniformSimpleTable<float, 1u, 2u> v_resonant_table;

		float k_prev = 0.f, v_resonant; float v_begin = -12.5f, v_step = 2.5e-2f; unsigned counter = 0u;
		while (true) {
			v_resonant = v_begin + counter * v_step;
			auto k_omega_opt = DRResonantVelocitySolve(v_resonant, params, zfunc);
			if (!k_omega_opt) {
				std::stringstream error_text;
				error_text << "Can't find root for the resonant velocity " << v_resonant;
				throw std::runtime_error(error_text.str());
			}
			if (k_omega_opt->first < k_prev) break;
			k_prev = k_omega_opt->first;
			v_resonant_data.push_back(k_omega_opt->second);
			v_resonant_data.push_back(k_omega_opt->first);
			++counter;
		}
		cout << v_resonant_data.size() << " " << counter << endl;
		reverse(v_resonant_data.begin(), v_resonant_data.end());
		v_resonant_table.space.axes[0].begin = v_resonant;
		v_resonant_table.space.axes[0].step = -v_step;
		v_resonant_table.bounds.components[0] = v_resonant_data.size() / 2u;
		v_resonant_table.data = v_resonant_data.data();

		{
			std::ofstream ascii_os("./fResonantVelocity.txt");
			ascii_os << setprecision(7) << fixed;
			ascii_os << v_resonant_table;
		}
	}
	catch (std::runtime_error const &ex) {
		cerr << ex.what() << endl;
	}
}