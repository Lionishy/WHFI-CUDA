#include "VparallGridCalculation.h"
#include "ZFunc.h"
#include "ZFuncIO.h"

#include <iostream>
#include <fstream>
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
	whfi::PhysicalParamenters<float> params = whfi::init_parameters(0.85f, 1.f / 0.85f, 0.25f, -9.f);

	/* Set Grid Calculator */
	auto vparall_grid_calculator = whfi::VparallGridCalculator(zfunc, params);

	/* Vparall Axis */
	UniformAxis<float> vparall_axis = { -12.5f, 2.4e-2f };

	/* Result table */
	std::vector<float> vparall_data;
	auto vparall_table = vparall_grid_calculator(vparall_axis, vparall_data);

	{
		ofstream ascii_os("./data/fVparallGrid.txt");
		ascii_os << setprecision(7) << fixed << vparall_table;
	}
}