
#include "SimpleTable.h"
#include "SimpleTableIO.h"
#include "ZFunc.h"
#include "PhysicalParameters.h"
#include "DispersionRelation.h"

#include <vector>
#include <fstream>
#include <exception>
#include <stdexcept>

void KGridSolverHost() {
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
	whfi::PhysicalParamenters<float> params = whfi::init_parameters(0.85f,1.f/0.85f,0.25f,-9.f);
	auto dispersion_relation = whfi::DispersionRelation(zfunc,params);


}