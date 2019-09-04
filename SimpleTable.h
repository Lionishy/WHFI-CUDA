#pragma once
#ifndef SimpleTable_H
#define SimpleTable_H

#include <cstddef>

namespace iki {
	template <size_t Dim>
	struct Index {
		size_t components[Dim];
	};

	template <typename T, size_t Dim>
	struct Argument {
		T components[Dim];
	};

	template <typename T>
	struct UniformAxis {
		size_t size = 0u; T begin = T(0), step = T(1);
	};

	template <typename T, size_t Dim>
	struct UniformSpace {
		UniformAxis<T> axes[Dim];
	};

	template <typename T, size_t Dim, size_t Scale>
	struct SimpleTable {
		T *data;
		size_t bounds[Dim];
	};
} /* iki */

#endif
