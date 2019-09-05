
#include "KahanSumSequence.cuh"

#include <cuda_runtime.h>

template <typename T>
struct GeometricProgression {
	__device__ T operator()(size_t counter, T _) {
		T tmp = current; current *= q; return tmp;
	}

	__device__ GeometricProgression(T q, T q0): q(q), q0(q0), current(q0) { }
	__device__ void reset() { current = q0; }
private:
	T const q;
	T const q0;
	T current;
};

template <typename T>
struct HarmonicSeries {
	__device__ T operator()(size_t counter, T _) {
		return counter == 0 ? T(1.) : T(1. / counter);
	}

	__device__ HarmonicSeries() { }
	__device__ void reset() { }

};

template <typename T>
struct ZFuncSeries {
	__device__ T operator()(size_t counter, T s) {
		return - ((counter * step) * s + T(1.))*step;
	}

	__device__ ZFuncSeries(T step): step(step) { }

private:
	T const step;
};

template <typename T>
__global__ void test_summation_kernell_geometric_progression(T *res) {
	*res = iki::math::device::kahan_summation_sequence<T, GeometricProgression<T>>(GeometricProgression<T>(T(0.5),T(0.5)), 1'000u);
}

template <typename T>
__global__ void test_summation_kernell_harmonic_series(T *res) {
	*res = iki::math::device::kahan_summation_sequence<T, HarmonicSeries<T>>(HarmonicSeries<T>(), 100'000'000u);
}


template <typename T>
__global__ void test_tabulator_kernell_zfunc(T *res, size_t size) {
	*res = T(0.);
	iki::math::device::kahan_tabulator_sequence<T, ZFuncSeries<T>>(ZFuncSeries<T>(1.e-2f), size, res+1);
}

#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

void run_geometic_progression_test() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto End;
	}

	{
		float *res_dev_ptr, res_host;
		if (cudaSuccess != cudaMalloc((void **)& res_dev_ptr, sizeof(float))) {
			cout << "Can't alocate device memory: " << sizeof(float) << " bytes" << endl;
			goto End;
		}
		test_summation_kernell_geometric_progression<float><<<1,1>>>(res_dev_ptr);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			cout << "geometric progression kernell launch failed: " << cudaGetErrorString(cudaStatus) << endl;
		}
		else {
			cudaMemcpy(&res_host, res_dev_ptr, sizeof(float), cudaMemcpyDeviceToHost);
			cout << res_host << endl;
		}
		cudaFree(res_dev_ptr);
	}

End:
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		cout << "Error in cuda device reset" << endl;
	}
}

void run_harmonic_series_test() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto End;
	}

	{
		float *res_dev_ptr, res_host;
		if (cudaSuccess != cudaMalloc((void **)& res_dev_ptr, sizeof(float))) {
			cout << "Can't alocate device memory: " << sizeof(float) << " bytes" << endl;
			goto End;
		}
		test_summation_kernell_harmonic_series<float><<<1,1>>>(res_dev_ptr);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			cout << "harmonic series kernell launch failed: " << cudaGetErrorString(cudaStatus) << endl;
		}
		else {
			cudaMemcpy(&res_host, res_dev_ptr, sizeof(float), cudaMemcpyDeviceToHost);
			cout << res_host << endl;
		}
		cudaFree(res_dev_ptr);
	}

End:
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		cout << "Error in cuda device reset" << endl;
	}
}

void run_zfunc_tabulator_test() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto End;
	}

	{
		float *res_dev_ptr; size_t size = 1000u;
		if (cudaSuccess != cudaMalloc((void **)& res_dev_ptr, size * sizeof(float))) {
			cout << "Can't alocate device memory: " << size*sizeof(float) << " bytes" << endl;
			goto End;
		}
		test_tabulator_kernell_zfunc<float><<<1,1>>>(res_dev_ptr,size);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			cout << "harmonic series kernell launch failed: " << cudaGetErrorString(cudaStatus) << endl;
		}
		else {
			std::vector<float> res_host(size);
			cudaMemcpy(res_host.data(), res_dev_ptr, size * sizeof(float), cudaMemcpyDeviceToHost);
			{
				ofstream ascii_os("./fZFunc.txt");
				size_t counter = 0u; float const step = 1.e-2f;
				for_each(begin(res_host), end(res_host), [&ascii_os, &counter, step] (auto x) { ascii_os << step * (counter++) << ' ' << x << '\n'; });
			}
		}
		cudaFree(res_dev_ptr);
	}

End:
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		cout << "Error in reseting cuda device" << endl;
	}
}