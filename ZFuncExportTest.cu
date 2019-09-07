#include "RungeKutta.cuh"
#include "KahanSumSequence.cuh"
#include "ZFunc.cuh"
#include "SimpleTable.h"
#include "SimpleTableIO.h"

#include <cuda_runtime.h>

template <typename T>
__global__ void test_zfunc_table_calculation(iki::UniformSimpleTable<T,1u,1u> res_table, size_t loop_count = 1u) {
	using namespace iki::math::device;

	*res_table.data = T(0.);
	kahan_tabulator_sequence(make_runge4th(T(0.), res_table.space.axes[0].step/loop_count, iki::whfi::device::ZFuncODE<T>()), res_table.bounds.components[0]-1u, res_table.data + 1u, loop_count);
}

#include <vector>
#include <fstream>
#include <iomanip>

void zfunc_export_test() {
	using namespace std;

	cudaError_t cudaStatus;
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		cout << "Error in starting cuda device" << endl;
		goto End;
	}

	{
		iki::UniformSimpleTable<float, 1u, 1u> ZFuncTable;
		ZFuncTable.space.axes[0] = { 0.f, 1.e-4f };
		ZFuncTable.bounds.components[0] = 100'001u;

		vector<float> hst_table_data(ZFuncTable.bounds.components[0]);
		float *dev_table_data;
		if (cudaSuccess != cudaMalloc((void **)& dev_table_data, ZFuncTable.bounds.components[0] * sizeof(float))) {
			cout << "Can't alocate device memory: " << ZFuncTable.bounds.components[0] * sizeof(float) << " bytes" << endl;
			goto End;
		}

		ZFuncTable.data = dev_table_data;
		test_zfunc_table_calculation<float> <<<1, 1>>> (ZFuncTable, 100u);
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			cout << "ZFuncTable calculation kernell launch failed: " << cudaGetErrorString(cudaStatus) << endl;
		}
		else {
			if (cudaSuccess != cudaMemcpy(hst_table_data.data(), dev_table_data, ZFuncTable.bounds.components[0] * sizeof(float), cudaMemcpyDeviceToHost)) {
				cout << "Memory copy device->host failed!" << endl;
			}
			else {
				ZFuncTable.data = hst_table_data.data();

				ofstream ascii_os("./fZFunc.txt");
				ascii_os.precision(8); ascii_os.setf(ios::fixed);
				ascii_os << ZFuncTable << endl;
			}
		}
		cudaFree(dev_table_data);
	}

End:
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		cout << "Error in reseting cuda device" << endl;
	}
}