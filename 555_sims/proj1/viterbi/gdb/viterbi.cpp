#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "base.h"
#include "fsm.h"
#include "viterbi.h"
using namespace gr;
using namespace trellis;

vector<int> VA_encode(fsm cc,vector<int> data,int n)
{
	int state, in, next;
	vector<int> out(data.size());
	vector<int> outb(data.size()*n);
	vector<int> ns = cc.NS();
	vector<int> os = cc.OS();
	state=0; // inital state is 0
	for(int i=0;i<data.size();i++) {
		in=data[i]; out[i]=os[state*cc.I()+in]; next=ns[state*cc.I()+in];
		//printf("state:%d, in:%d -> out:%d, next_state:%d\n",state,in,out[i],next);
		state=next;
	}
	outb = int2bin(out,n);
	return outb;
}

vector<int> VA_decode(fsm cc,vector<int> data,int n)
{
	// notes:
	// assign initial distance weights to all states: state0=0, all else=inf
	// collect nv bits at time t to for a rx'd word
	// for cc.S()
		// go to each state and calculate distances associated w/ all paths going into each state
		// for every state determine which of the cc.S() survivor paths it belongs to
		// append the state to one of the cc.S() survivor paths
	// build up a list of cc.S() number of survivor paths until assessment is finished.
		// this list will contain the survivor path's estimated bit stream and associated distance
	// at end time T choose the bit stream with the lowest accumulated associated survivor path distance


	// old:
	// generate outputs that are a fcn of the current state, 0 bit input, and 1 bit input
	// calculate distance b/w 0-bit output, and 1-bit output and pick the minimum
	// store input data as the estimated data stream
	// generate the next state as fcn of the current state and the chosen input data
	// repeat

	int state, in, next, w_iter=0, dist0=0, dist1=0, b_iter=0;
	vector<int> word(n), ns=cc.NS(), os=cc.OS(), b_est(data.size()/n);

	state=0; // assume initial state is 0
	for(int i=0;i<data.size();i++) {
		word[w_iter]=data[i];
		if(w_iter==(n-1)) {
			vector<int> test0(n); dec2base(os[state*cc.I()+0],2,test0);
			vector<int> test1(n); dec2base(os[state*cc.I()+1],2,test1);
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

vector<int> int2bin(vector<int> in,int k)
{
	// Note: For the purposes of converting in the comms simulation,
	//	 I don't need a binary number exceeds k (i.e. log2(M)).
	//	 Therefore, the needed vector is will be
	//	 K = k*ns = log2(M)*(#_of_symbols).
	// Note: Input = vector<int> of size ns, Output = vector<int> or vector<bool> of size K = k*ns

	int h,rem,dec,ns=in.size();
	vector<int> out(ns*k);
	for (int n=0;n<ns;n++) {
		h=(n+1)*k-1;
		dec=in[n];
		do {rem=dec%2; out[h]=rem; dec=dec/2; h--;}
		while(dec>0); //end_dowhile
	} //end_for
	return out;
}
