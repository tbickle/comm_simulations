#include <iostream>
#include <random>
#include <math.h>
using namespace std;

// global variable(s)
const int K = 1000000; // # of tests & memory limitation

// function prototype(s)
void randn(float*, float, float);
float qfunc(float, float*);
float BER(float, int*, float, float);

int main(void)
{

// select transmission and channel parameters
	float Es = 3;
	float Eb = Es;
	float N0 = 0.5;    
	float mean = 0.0; // gaussian mean


// generate vector of gaussian random numbers
	float rnd[K] = {}; // vector of gaussian distributed random numbers
	randn(rnd, mean, 1.0); // generate Gaussian vector
	//for (int i=0; i<K; i++) {cout << rnd[i] << endl;} // print rand


// generate approximate qfunction: Q(x)
	//float Q = 0.0;
	//float x = 0.0;
	//Q = qfunc(x, n);
	//cout << Q << endl;


// generate cvs file to print graph
	// ???


// generate random bit stream
	int b[K] = {}; // random bit stream (vector)
	for (int i=0; i<K; i++) {
		if (rnd[i] > 0) {b[i] = 1;}
		else {b[i] = 0;}
	}
	//for (int i=0; i<K; i++) {cout << b[i] << endl;} // print b


// run BER trials

	// generate BER 
	float ber = 0.0;
	ber = BER(Es, b, mean, N0);
	cout << ber << endl;




// generate actual curve
	//p = qfunc(sqrt(2*Eb/N0));
    
// plot results
	//figure(2);
	//semilogy(Eb/N0,BER,Eb/N0,p);

return 0;
}


// **************************************
// ************* FUNCTIONS **************
// **************************************

void randn(float *p, float mean, float std)
{
	// input: gaussian mean and standard deviation
	// output: a vector of gaussian distributed random numbers

	random_device rd;
	mt19937 e2(rd());
	normal_distribution<float> dist(mean,std);
	float number = 0.0;

	for (int i=0; i<K; i++) {
		number = dist(e2);
		*p = number;
		p++;
	}
}


float qfunc(float x, float *rnd)
{
	// An single estimation generator of the Q(.) function.
	// theory: Q(x) ~ (# of realizations exceeding x)/K
	// x = single input value
	// rnd = is the randomly generated vector of N.Guassian realizations

	int count = 0;
	for(int i=0; i<K; i++) {
		if (*rnd > x) {count++;}	
		rnd++;
	}

	return count/(K*1.0); // divide by float to return a float
}

float BER(float Es, int *b, float mean, float N0)
{
	// Generates Bit Error Rate
	// Es = symbol energy
	// b = transmitted bit stream
	// n = channel noise

	float s[K] = {};
	for (int i=0; i<K; i++) {
	// i'm declaring Es = 1 -> b = 1
	// and Es = -1 -> b = 0
		if (b[i] == 1) {s[i] = sqrt(Es);}
		else {s[i] = -sqrt(Es);}
	}	

	// generate noise
	float var = N0/2.0; // variance
	float std = sqrt(var); // gaussian standard deviation
	float n[K] = {}; // vector of gaussian distributed random numbers
	randn(n, mean, std); // generate Gaussian vector of variance N0/2

	// generate observation equation
	float z[K] = {};
	for (int i=0; i<K; i++) {z[i] = s[i] + n[i];}

	// MAP detection to determine b_est
	float b_est[K] = {};
	for (int i=0; i<K; i++) {if (z[i] > 0) {b_est[i] = 1;}}

	// count differences b/w b & b_est to determine BER
	int temp = 0;
	for (int i=0; i<K; i++) {if (b[i] != b_est[i]) {temp++;}}

	// test outputs
	//for (int i=0; i<K; i++) {
		//cout << s[i] << endl; // print s
		//cout << n[i] << endl; // print n
		//cout << z[i] << endl; // print z
		//cout << b_est[i] << endl; // print b_est
	//}

	return temp/(K*1.0);
}

