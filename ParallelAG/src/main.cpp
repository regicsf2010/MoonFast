#include <iostream>
#include <chrono>
using namespace std;
#include "auxiliaries/Configuration.h"
#include "ag/MoonFast.h"
#include "ag/Population.h"

int main(void) {
	MoonFast mf(Ackley::ID);

	auto start = chrono::steady_clock::now();
	Population *p = mf.run();
	auto finish = chrono::steady_clock::now();

	//double elapsed = chrono::duration_cast<chrono::milliseconds>(finish - start).count();
	double elapsed = chrono::duration_cast<chrono::duration<double>>(finish - start).count();

	cout << "Time: " << elapsed << " (s)" << endl;

	delete p;
	cout << "--END MOONFAST--" << endl;
	return 0;
}





