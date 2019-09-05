#pragma once
#ifndef KahanSumSequence_H
#define KahanSumSequence_H

#include <cuda_runtime.h>

#include <cstddef>

namespace iki { namespace math { namespace device { 
	template <typename T, typename Seq>
	__device__ T kahan_summation_sequence(Seq seq, size_t seq_size) {
		T s = T(0), c = T(0);
		for (size_t counter = 0u; counter != seq_size; ++counter) {
			T y, t;
			y = seq(counter,s) - c;
			t = s + y;
			c = (t - s) - y;
			s = t;
		}
		return s;
	}

	template <typename T, typename Seq>
	__device__ void kahan_tabulator_sequence(Seq seq, size_t seq_size, T *begin, size_t loop_size = 1u) {
		T s = T(0), c = T(0);
		for (size_t seq_counter = 0u; seq_counter != seq_size; ) {
			for (size_t loop_count = 0u; loop_count != loop_size; ++loop_count, ++seq_counter) {
				T y, t;
				y = seq(seq_counter, s) - c;
				t = s + y;
				c = (t - s) - y;
				s = t;
			}
			*begin++ = s;
		}
	}
} /* device */ } /* math */ } /* iki */

#endif /* KahanSumSequence_H */
