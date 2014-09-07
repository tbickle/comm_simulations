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
int dist(int,vi);

int main(void)
{
	cout<<endl;
	// define convolutional code and FSM
	int kv=1, nv=3;		// kv=k(viterbi)=#_of_blocks???, nv=n(viterbi)=#_of_output_bits->(k,n)_cc=(kv,nv)_cc
	vi G = {7,7,5};		// generator matrix (g1,g2,g3)_10 
	fsm cc(kv,nv,G);	// fsm instance of (k,n) conv. code: (kv=k,nv=n,G=generator)
	vector<int> ns = cc.NS();
	vector<int> os = cc.OS();

	//vi b = {1,1,0,0,0,1,1,0,0,0}; 	// original data bit stream
	//vi b_coded = VA_encode(cc,b,nv);	// encode data with viterbi algorithm

	//vi b = {1,1,0,0,0}; 	// original data bit stream
	//vi b_coded = {1,0,1,0,0,1,1,0,1,1,1,1,0,0,0};

	vi b={1}; vi b_coded={1,0,1};

	// modulate->channel->demodulate
	//
	//vi b_est = VA_decode(cc,b_coded,nv); // decode



	// *** Initialize ***
	// cc.S() x b.size()+1 matrix:
		// the 1st element in each of the cc.S() rows will be the accumulated sp distance information
		// the rest of the columns in each row will contain the sp's b_est data stream
		// there will be 1 of these row vectors for each cc.S() survivor path
	// initialize zero state w/ dist 0, and the rest w/ a large #
	// check initializations
	vector<vi> sp;
	int rnum = cc.S(); // # of rows
	int cnum = b.size()+1; // # of columns
	sp.resize(rnum,vi(cnum));
	int hi = 50; // arbitrarily high #
	sp[0][0]=0; // initialize zero state to dist. 0
	for (int r=1;r<rnum;r++) sp[r][0]=hi; // initialize rest of states to a high number
	//for (int r=0;r<rnum;r++) {for (int c=0;c<cnum;c++) {cout<<sp[r][c]; if(!c) cout<<"\t";} cout<<endl;} // test
	// *** end_INIT ***
	vector<vi> pS = cc.PS(); // pS[curr_state][trans#]
	vector<vi> pI = cc.PI(); // pI[curr_state][prev_input]
	//for (int h=0;h<pS.size();h++) {for(int i=0;i<pS[0].size();i++) cout<<pS[h][i]; cout<<endl;} // test
	//for (int h=0;h<pI.size();h++) {for(int i=0;i<pI[0].size();i++) cout<<pI[h][i]; cout<<endl;} // test
	cout<<endl; // test
/*
	for (int h=0;h<cc.S();h++) { // current state
		cout<<"curr_state:"<<h<<endl;
		for(int i=0;i<cc.I();i++) { // transition #
			cout<<"trans#"<<i<<":"<<" PS:"<<pS[h][i]<<" -> "<<os[pS[h][i]*cc.I()+pI[h][i]]<<"/"<<pI[h][i]<<endl;
		} cout<<endl;
	} // test
*/
	
	// *** ACS (add-compare-select) ***
	//for (int g=1;g<=b.size();g++) { // for discrete time t (up to T)
	//int state=0;
	vi state(cc.S());
	int w_iter=0, b_iter=1;
	for(int g=0;g<b_coded.size();g++) { // for discrete time t (up to T)
		vi word(nv);
		word[w_iter]=b_coded[g];
		if(w_iter==(nv-1)) {
			for (int h=0;h<cc.S();h++) { // for each survivor path
				int tran, d, in;
				int tran0=os[pS[h][0]*cc.I()+pI[h][0]];
				int tran1=os[pS[h][1]*cc.I()+pI[h][1]];
				int d0 = dist(tran0,word); // find the dist. b/w transition0 output and rx'd word
				int d1 = dist(tran1,word); // find the dist. b/w transition1 output and rx'd word
				if(d0<=d1) {tran=0;d=d0;in=pI[h][0];} // determine minimum
				else {tran=1;d=d1;in=pI[h][1];}

				// log accumulated distance, and prev_input
				sp[h][0]=sp[h][0]+d; // accum dist
				sp[h][b_iter]=in; // log estimated data
				//state[h]=ns[state*cc.I()+in];
			} // end_for
			b_iter++;
			w_iter=0;
		} // end_if
		else {w_iter++;}
	} // end_for
	// *** end_ACS ***
	for (int r=0;r<rnum;r++) {for (int c=0;c<cnum;c++) {cout<<sp[r][c]; if(!c) cout<<"\t";} cout<<endl;} // test

	// *** Traceback ***
	// xxx
	// *** end_TB ***

	// display
	cout<<endl;
	cout<<"      b:"; for(int i=0;i<b.size();i++) cout<<b[i]; cout<<endl;		// test
	cout<<"b_coded:"; for(int i=0;i<b_coded.size();i++) cout<<b_coded[i]; cout<<endl;	// test
	//cout<<"b_est:"; for(int i=0;i<b_est.size();i++)  cout<<b_est[i]; cout<<endl;	// test
	cout<<endl;

	return 0;
}

int dist(int t,vi word)
{
	int distance=0;
	int n = word.size();
	vi dist(n);
	vi trans(n); dec2base(t,2,trans);
	for(int h=0;h<n;h++) { // calc. distance
		dist[h]=abs(word[h]-trans[h]);
	}
	distance=accumulate(&dist[0],&dist[n-1],dist[n-1]);
	return distance;
}



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

// ***** PSEUDOCODE *****
// initialize: s_0=0, s_else=inf <- min. dist. assigned to the zero state (assumption)
// for t=1:T
	// for cc.S() <- # of states
		// for cc.I() <- # of input transitions
			// determine distances of cc.I() transitions into cc.S() states
		// end_for
		// choose path to current state w/ lowest accumulated distance (survivor path)
		// append appropriate infomation to the chosen survivor path
	// end_for
// end_for
// sp_fin = min(sp_1:sp_cc.S()) <- choose the minimum distance survivor path
// b_est = input_code(sp_fin) <- choose the estimated bit stream associated w/ the chosen survivor path
// output b_est
// **********************

