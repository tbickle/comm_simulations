// normal_distribution

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
using namespace std;

// function prototype(s)
vector<double> randn(int,double,double);
vector<int> sign(vector<double>);
vector<int> bitstream(vector<double>);

int main(void)
{
	int K = 10; // # of tests
	vector<double> p(K);
	vector<int> sgn(K);
	vector<int> b(K);

	// generate vector of gaussian random numbers
	p = randn(K,0.0,1.0);
	sgn = sign(p);
	b = bitstream(p);

	//for (int i=0; i<K; i++) {cout<<p[i]<<endl;}	
	//for (int i=0; i<K; i++) {cout<<sgn[i]<<endl;}
	for (int i=0; i<K; i++) {cout<<b[i]<<endl;}	

	return 0;
}

// FUNCTIONS
vector<double> randn(int K,double mean,double std)
{
	// input: gaussian mean and standard deviation
	// output: a vector of gaussian distributed random numbers

	vector<double> p(K);
	random_device rd;
	mt19937 e2(rd());
	normal_distribution<double> dist(mean,std);

	for (int i=0; i<K; i++) {p[i]=dist(e2);}

	return p;
}

vector<int> sign(vector<double> in)
{
	vector<int> out(in.size());
	for (int n=0; n<in.size(); n++) {out[n]=copysign(1,in[n]);}
	return out;
}

vector<int> bitstream(vector<double> in)
{
	vector<int> out(in.size());
	for (int n=0; n<in.size(); n++) {out[n]=(copysign(1,in[n])+1)/2;}
	return out;
}

