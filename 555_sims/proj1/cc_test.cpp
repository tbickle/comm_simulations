#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "../include/base.h"
#include "../include/fsm.h"
#include "../include/viterbi.h"
using namespace std;
using namespace gr;
using namespace trellis;
typedef vector<int> vi;

int main(void)
{
	// define convolutional code and FSM
	int kv=1, nv=3;		// kv=k(viterbi)=#_of_blocks???, nv=n(viterbi)=#_of_output_bits->(k,n)_cc=(kv,nv)_cc
	vi G = {7,7,5};		// generator matrix (g1,g2,g3)_10 
	fsm cc(kv,nv,G);	// fsm instance of (k,n) conv. code: (kv=k,nv=n,G=generator)

	vi b = {1,1,0,0,0,1,1,0,0,0}; 	// original data bit stream
	vi coded = VA_encode(cc,b,nv);	// encode data with viterbi algorithm

	// modulate->channel->demodulate

	vi b_est = VA_decode(cc,coded,nv); // decode

	// display
	cout<<"    b:"; for(int i=0;i<b.size();i++) cout<<b[i]; cout<<endl;		// test
	cout<<"coded:"; for(int i=0;i<coded.size();i++) cout<<coded[i]; cout<<endl;	// test
	cout<<"b_est:"; for(int i=0;i<b_est.size();i++)  cout<<b_est[i]; cout<<endl;	// test

	return 0;
}

