#include <iostream>
#include <exception>

void ResonantVelocitySolverHost();
int main() {
	using namespace std;
	try {
		ResonantVelocitySolverHost();
	}
	catch (exception &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}