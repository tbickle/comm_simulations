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

//vi VA_decode(fsm,vi,int);

int main(void)
{
	// define convolutional code
	int kv=1, nv=3;		// kv = k (viterbi) = # of blocks???, nv = n (viterbi) = # of output bits -> (k,n) cc = (kv,nv) cc
	vi G = {7,7,5};		// generator matrix (g1,g2,g3)_10 
	fsm cc(kv,nv,G);	// fsm instance of (k,n) conv. code: (kv=k,nv=n,G=generator)

	// original data
	vi b = {1,1,0,0,0,1}; 	// data bit stream
	// viterbi algorithm encoded data
	vi coded = VA_encode(cc,b,nv); // encode

	// modulate
	// channel
	// demodulate

	vi b_est = VA_decode(cc,coded,nv); // decode

	// display
	cout<<"    b:"; for(int i=0;i<b.size();i++) cout<<b[i]; cout<<endl;		// test
	cout<<"coded:"; for(int i=0;i<coded.size();i++) cout<<coded[i]; cout<<endl;	// test
	cout<<"b_est:"; for(int i=0;i<b_est.size();i++)  cout<<b_est[i]; cout<<endl;

	return 0;
}
/*
vi VA_decode(fsm cc,vi data,int n)
{
	// notes:
	// collect nv bits at a time of the rx'd signal
	// generate outputs that are a fcn of the current state, 0 bit input, and 1 bit input
	// calculate distance b/w 0-bit output, and 1-bit output and pick the minimum
	// store input data as the estimated data stream
	// generate the next state as fcn of the current state and the chosen input data
	// repeat

	int state, in, next, w_iter=0, dist0=0, dist1=0, b_iter=0;
	vi word(n), ns=cc.NS(), os=cc.OS(), b_est(data.size()/n);

	state=0; // assume initial state is 0
	for(int i=0;i<data.size();i++) {
		word[w_iter]=data[i];
		if(w_iter==(n-1)) {
			vi test0(n); dec2base(os[state*cc.I()+0],2,test0);
			vi test1(n); dec2base(os[state*cc.I()+1],2,test1);
			for(int h=0;h<word.size();h++) { // calc. distance b/w test0-word, & test1-word
				test0[h]=abs(word[h]-test0[h]);
				test1[h]=abs(word[h]-test1[h]);
			}
			dist0=accumulate(&test0[0],&test0[n-1],test0[n-1]);
			dist1=accumulate(&test1[0],&test1[n-1],test1[n-1]);
			if(dist0<=dist1) in=0; else in=1;
			b_est[b_iter]=in; // min. dist. est.
			b_iter++;
			state=ns[state*cc.I()+in];
			w_iter=0;
		} //endif
		else {w_iter++;}
	}
	return b_est;
}
*/
