
#include "SimpleTable.h"
#include "SimpleTableIO.h"
#include "ZFunc.h"
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

void ResonantVelocitySolverHost() {
	using namespace std;
	using namespace iki;

	/* Load ZFunc table */
	whfi::ZFunc<float> zfunc;
	vector<float> zfunc_table_data;
	{
		zfunc.zfunc_table.space.axes[0] = { 0.f, 1.e-4f };
		zfunc.zfunc_table.bounds.components[0] = 100'001u;
		zfunc_table_data.resize(zfunc.zfunc_table.bounds.components[0]);
		zfunc.zfunc_table.data = zfunc_table_data.data();
		{
			ifstream binary_is;
			binary_is.exceptions(ios::badbit | ios::failbit);
			binary_is.open("./fZfunc.tbl", ios::binary);
			read_binary(binary_is, zfunc.zfunc_table);
		}
	}

	/* Set physical parameters */
	whfi::PhysicalParamenters<float> params = whfi::init_parameters(0.85f, 1.f / 0.85f, 0.25f, -9.f);
	
	std::vector<std::pair<float, float>> k_omega;
	try {
		float k_prev = 0.f;
		unsigned counter = 0u; float v_resonant_start = -12.5f, v_resonant_step = 2.5e-2f;
		for (; counter < 512u; ++counter) {
			float v_resonant = v_resonant_start + counter * v_resonant_step;
			auto k_omega_opt = DRResonantVelocitySolve(v_resonant, params, zfunc);
			if (!k_omega_opt) {
				std::stringstream error_text;
				error_text << "Can't find root for the resonant velocity " << v_resonant;
				throw std::runtime_error(error_text.str());
			}
			if (k_omega_opt->first < k_prev) break;
			k_prev = k_omega_opt->first;
			k_omega.push_back({ k_omega_opt->first, k_omega_opt->second });
		}
	}
	catch (std::runtime_error const &ex) {
		cerr << ex.what() << endl;
	}

	{
		ofstream ascii_out("./fKOmegaResonant.txt");
		ascii_out << setprecision(7) << fixed;
		for (auto &[k, omega] : k_omega)
			ascii_out << k << " " << omega << endl;
	}
	
	
	
}