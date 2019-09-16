#pragma once
#ifndef ZFuncIO_H
#define ZFuncIO_H

#include "SimpleTable.h"
#include "SimpleTableIO.h"
#include "ZFunc.h"
#include "RungeKutta.h"
#include "KahanSumSequence.h"

#include <vector>
#include <iostream>

namespace iki {	namespace whfi {
	template <typename T>
	std::istream& ZFuncImport(std::istream &binary_is, UniformSimpleTable<T, 1u, 1u> & zfunc_table, std::vector<T> &zfunc_data) {
		read_binary(binary_is, zfunc_table.space);
		read_binary(binary_is, zfunc_table.bounds);

		zfunc_data.resize(zfunc_table.bounds.components[0]);
		zfunc_table.data = zfunc_data.data();
		read_binary(binary_is, zfunc_table);
		return binary_is;
	}

	template <typename T>
	std::ostream& ZFuncExport(std::ostream &binary_os, T step, unsigned size, unsigned loop_size = 1u) {
		std::vector<float> zfunc_data;
		iki::UniformSimpleTable<float, 1u, 1u> zfunc_table;
		{
			zfunc_table.bounds.components[0] = size;
			zfunc_table.space.axes[0].begin = 0.f;
			zfunc_table.space.axes[0].step = step;
			zfunc_data.resize(zfunc_table.bounds.components[0]);
			zfunc_table.data = zfunc_data.data();
			zfunc_table.data[0] = 0.;
		}
		iki::math::kahan_tabulator_sequence(
			iki::math::make_runge4th<float>(0.f, step/loop_size, iki::whfi::ZFuncODE<float>())
			, size - 1u, zfunc_table.data + 1u, loop_size
		);

		{
			write_binary(binary_os, zfunc_table.space);
			write_binary(binary_os, zfunc_table.bounds);
			write_binary(binary_os, zfunc_table);
		}

		return binary_os;
	}
} /* whfi */ } /* iki */

#endif