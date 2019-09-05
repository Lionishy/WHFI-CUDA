#include "RungeKutta.cuh"
#include "KahanSumSequence.cuh"
#include "ZFunc.cuh"

#include <cuda_runtime.h>

template <typename T>
__global__ void test_summation_kernell_runge_kutta_zfunc(T *res, size_t size, T step, size_t loop_count = 1u) {
	using namespace iki::math::device;
	*res = T(1.);
	kahan_tabulator_sequence(make_runge4th(T(0.),step,iki::whfi::device::ZFuncODE<T>()), size-1, res + 1,loop_count);
}

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

void run_zfunc_runge_kutta_test() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto End;
	}

	{
		float *res_dev_ptr; size_t size = 100'001u;
		if (cudaSuccess != cudaMalloc((void **)& res_dev_ptr, size * sizeof(float))) {
			cout << "Can't alocate device memory: " << size * sizeof(float) << " bytes" << endl;
			goto End;
		}
		test_summation_kernell_runge_kutta_zfunc<float><<<1, 1>>>(res_dev_ptr, size, 1.e-6f, 100u);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			cout << "harmonic series kernell launch failed: " << cudaGetErrorString(cudaStatus) << endl;
		}
		else {
			std::vector<float> res_host(size);
			if (cudaSuccess != cudaMemcpy(res_host.data(), res_dev_ptr, size*sizeof(float), cudaMemcpyDeviceToHost)) {
				cout << "Memory copy device->host failed!" << endl;
			} 
			else {
				ofstream ascii_os("./fZFunc.txt");
				size_t counter = 0u; float const step = 1.e-4f;
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