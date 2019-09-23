#include "VparallGridCalculation.h"
#include "ZFunc.h"
#include "ZFuncIO.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void DiffusionCoefficientsKernellHost() {
	using namespace std;
	using namespace iki;
	
	/* Load ZFunc table */
	whfi::ZFunc<float> zfunc;
	vector<float> zfunc_data;
	{
		std::ifstream binary_is("./data/fZFunc.tbl", ios::binary);
		binary_is.exceptions(ios::badbit | ios::failbit);
		whfi::ZFuncImport(binary_is, zfunc.zfunc_table, zfunc_data);
	}

	/* Set physical parameters */
	whfi::PhysicalParamenters<float> params = whfi::init_parameters(0.85f, 1.f / 0.85f, 0.25f, -11.f);

	/* Set Grid Calculator */
	auto vparall_grid_calculator = whfi::VparallGridCalculator(zfunc, params);

	/* Vparall Axis */
	UniformAxis<float> vparall_axis = { -12.5f, 2.4e-2f };

	/* Result table */
	std::vector<float> vparall_data;
	auto vparall_table = vparall_grid_calculator(vparall_axis, vparall_data);

	{
		stringstream filename_ss;
		filename_ss << setprecision(2) << fixed;
		filename_ss << "./data/fVparallGrid" << "-" << params.nc << "-" << params.betta_c << "-" << params.TcTh_ratio << "-" << (params.bulk_to_alfven_c < 0 ? "m" : "") << std::fabs(params.bulk_to_alfven_c) << ".txt";
		ofstream ascii_os(filename_ss.str());
		ascii_os << setprecision(7) << fixed << vparall_table;
	}
}