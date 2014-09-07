// complex vector test

#include <complex>
#include <iostream>
#include <vector>
using namespace std;

/*
typedef vector<complex<double>> dcomplex;

int main()
{
	vector<dcomplex> v;
	//dcomplex v;
	int rnum = 1; // # of rows
	int cnum = 3; // # of columns

	v.resize(rnum,dcomplex(cnum));

	//v[0][0] = complex<double>(1.0,2.0);
	//v[0][1] = complex<double>(3.0,4.0);
	//v[0][2] = complex<double>(5.0,6.0);
	v[0] = {complex<double>(1,2),complex<double>(3,4),complex<double>(5,6)};

	for (int n=0; n<cnum; n++) {
		cout << v[0][n] << endl;
		cout << "real:" << real(v[0][n]) << " imag:" << imag(v[0][n]) << endl;
	}

	return 0;
}
*/


int main()
{
	typedef complex<double> cd;
	vector<cd> v;
	v.resize(3);

	//v[0] = complex<double>(1.0,2.0);
	//v[0] = cd(1.0,2.0);
	v = {cd(1.0,2.0),cd(3.0,4.0),cd(5.0,6.0)};

	for (int n=0; n<v.size(); n++) {
		cout << v[n] << endl;
		cout << "real:" << real(v[n]) << " imag:" << imag(v[n]) << endl;
	}

	return 0;
}
