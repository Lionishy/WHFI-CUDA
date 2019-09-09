#include <iostream>
#include <exception>

void KGridSolverHost();
int main() {
	using namespace std;
	try {
		KGridSolverHost();
	}
	catch (exception &ex) {
		cerr << ex.what() << endl;
	}
	return 0;
}