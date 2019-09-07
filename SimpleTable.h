#pragma once
#ifndef SimpleTable_H
#define SimpleTable_H

#include <cstddef>

namespace iki {
	/* Index and Bounds */
	template <size_t Dim>
	struct Index {
		size_t components[Dim];
	};

	template <size_t Dim>
	struct Bounds {
		size_t components[Dim];
	};

	template <size_t Dim>
	size_t collapsed_size(Bounds<Dim> const *bounds) {
		size_t size = 1u;
		for (size_t bound_idx = 0u; bound_idx != Dim; ++bound_idx)
			size *= bounds->components[bound_idx];
		return size;
	}

	template <size_t Dim>
	Index<Dim> collapsed_bounds(Bounds<Dim> const *bounds) {
		Index<Dim> collapsed;
		collapsed.components[0] = 1u;
		for (size_t collapsed_idx = 1u; collapsed_idx != Dim; ++collapsed_idx)
			collapsed[collapsed_idx] = collapsed[collapsed_idx - 1u] * bounds->components[collapsed_idx - 1u];
		return collapsed;
	}

	template <size_t Dim>
	size_t collpased_index(Index<Dim> const *idx, Index<Dim> const *collapsed_bounds) {
		size_t collapsed = 0u;
		for (size_t idx_idx = 0u; idx_idx != Dim; ++idx_idx)
			collapsed += idx[idx_idx] * collapsed_bounds[idx_idx];
		return collapsed;
	}

	template <size_t Dim>
	Index<Dim> first_index() {
		Index<Dim> expanded_index;
		for (size_t idx_idx = 0u; idx_idx != Dim; ++idx_idx)
			expanded_index.components[idx_idx] = 0u;
		return expanded_index;
	}

	template <size_t Dim>
	Index<Dim> last_index(Bounds<Dim> const *bounds) {
		Index<Dim> expanded_index;
		for (size_t idx_idx = 0u; idx_idx != Dim; ++idx_idx)
			expanded_index.components[idx_idx] = bounds->components[idx_idx] == 0u ? 0u : bounds->components[idx_idx] - 1u;
		return expanded_index;
	}

	template <size_t Dim>
	void next_index(Index<Dim> *index, Bounds<Dim> const *bounds) {
		bool carry_flat = true;
		for (size_t idx_idx = 0u; idx_idx != Dim; ++idx_idx)
			if (carry_flag) {
				index->components[idx_idx] += 1u;
				carry_flag = false;
				if (index->components[idx_idx] == bounds->components[idx_idx]) {
					index->components[idx_idx] = 0u;
					carry_flag = true;
				}
			}
	}

	/* Argument */
	template <typename T, size_t Dim>
	struct Argument {
		T components[Dim];
	};

	/* Uniform Space */
	template <typename T>
	struct UniformAxis {
		T begin = T(0), step = T(1);
	};

	template <typename T, size_t Dim>
	struct UniformSpace {
		UniformAxis<T> axes[Dim];
	};

	template <typename T, size_t Dim>
	Argument<T, Dim> uniform_argument(Index<Dim> const *index, UniformSpace<T,Dim> const *space) {
		Argument<T, Dim> expanded_argument;
		for (size_t arg_idx = 0u; arg_idx != Dim; ++arg_idx)
			expanded_argument[arg_idx] = space->axes[arg_idx].begin + space->axes[arg_idx].step * index->components[arg_idx];
		return expanded_argument;
	}

	/* Table */
	template <typename T, size_t Dim, size_t Scale>
	struct SimpleTable {
		T *data;
		Bounds<Dim> bounds;
	};

	template <typename T, size_t Dim, size_t Scale>
	struct UniformSimpleTable {
		T *data;
		Bounds<Dim> bounds;
		UniformSpace<T, Dim> space;
	};

	

	
} /* iki */

#endif
