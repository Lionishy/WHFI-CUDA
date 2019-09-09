
#include "SimpleTable.h"
#include "SimpleTableIO.h"
#include "ZFunc.h"
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
		zfunc_table.space.axes[0] = { 0.f, 1.e-4f };
		zfunc_table.bounds.components[0] = 100'001u;
		zfunc_table_data.resize(zfunc_table.bounds.components[0]);
		zfunc_table.data = zfunc_table_data.data();
		{
			ifstream binary_is;
			binary_is.exceptions(ios::badbit | ios::failbit);
			binary_is.open("./fZfunc.tbl", ios::binary);
			read_binary(binary_is, zfunc_table);
		}
	}

	
}