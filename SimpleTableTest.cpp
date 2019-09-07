#include "SimpleTable.h"
#include "SimpleTableIO.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

void simple_table_zfunc_test() {
	using namespace std;
	using namespace iki;

	vector<float> zfunc_table_data;
	UniformSimpleTable<float, 1u, 1u> zfunc_table;
	{
		cin >> zfunc_table.space;
		zfunc_table.bounds.components[0] = 100'001;
		zfunc_table_data.resize(zfunc_table.bounds.components[0]);
		zfunc_table.data = zfunc_table_data.data();
	}

	{
		ifstream ascii_is("./fZfunc.txt");
		ascii_is >> zfunc_table;
	}

	{
		ofstream ascii_os("./fZfunc2.txt");
		ascii_os << setprecision(7) << fixed;
		ascii_os << zfunc_table;
	}
}