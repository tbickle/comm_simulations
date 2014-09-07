#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iterator>
using namespace std;
typedef vector<int> vi;
typedef vector<double> vd;
typedef complex<double> cd;
typedef vector<complex<double>> vcd;

#ifndef COMM_H
#define COMM_H

vcd sym_table(int,int);
double rms(vcd);
vd randn(int,double,double);
vi bitstream(vd);
vi symbol_map(vi,int,int);
vcd crandn(int K,double mean,double std);
vi int2bin(vi,int);
void display(vi,vcd,vi,vi,vi,vcd,vcd,vcd,vi,vi,vi,vi,vi,vi);

#endif
