#include <iostream>
#include <exception>

void DiffusionCoefficientsKernellHost();
int main() {
	using namespace std;
	try {
		DiffusionCoefficientsKernellHost();
	}
	catch (exception &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}