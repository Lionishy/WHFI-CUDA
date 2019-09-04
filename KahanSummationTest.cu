
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
__global__ void test_kernell(T *res) {
	*res = iki::math::device::kahan_summation_sequence<T, GeometricProgression<T>>(GeometricProgression<T>(T(0.5),T(0.5)), 10);
}

#include <iostream>

float run_kernell() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto Error;
	}
	{
		float *res_dev_ptr; float res_host = -1.;
		cudaStatus = cudaMalloc((void **)&res_dev_ptr, sizeof(float));

		test_kernell<float> <<<1, 1>>> (res_dev_ptr);
		cudaMemcpy(&res_host, res_dev_ptr, sizeof(float), cudaMemcpyDeviceToHost);
		cout << "Res: " << res_host << endl;

		cudaFree(res_dev_ptr);
	}

	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto Error;
	}

Error:
	return -1.f;
}