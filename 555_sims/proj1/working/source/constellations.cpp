// Digital Communications Channel & Demodulation/Decoding Simulator
// by Matthew Chase
// tbickle.matt@gmail.com


// NOTES
// Compile: g++ -o t constellations.cpp -std=c++11
//    Call: ./t (mod scheme #:0=pam,1=psk,2=qam,3=orth) (M-ary #:2,4,8,16)
// Derivation of N0:
//	define: Es=1, k=log2(M), EbN0dB=0:#
//	find: N0
//	recognize Eb/N0 = Es/(k*N0)
//	EbN0dB = 10log10(Es/(k*N0)
//	10^(EbN0dB/10) = Es/(k*N0), therefore
//	N0 = 1/k*10^(-EbN0dB/10)
//int MOD; cout << "Enter Modulation Scheme (0=PAM,1=PSK,2=QAM,3=Orth): "; cin >> MOD;
//int M; cout << "Enter # of Symbols (2,4,8,16): "; cin >> M;

#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iterator>
#include "../include/comm.h"
#include "../include/base.h"
#include "../include/fsm.h"
#include "../include/viterbi.h"

using namespace std;
using namespace gr;
using namespace trellis;

typedef vector<int> vi;
typedef vector<double> vd;
typedef complex<double> cd;
typedef vector<complex<double>> vcd;

// global vars
const bool CSV=0; // simulation outputs in .csv format (0=No,1=Yes)
//double RMS=1;

int main(int argc, char *argv[])
{
	if(!CSV) cout<<endl;
	// INITIALIZATION
	int MOD=atoi(argv[1]), M=atoi(argv[2]);
	const int ORD=5, MAX=20;

	// ***** conv. code *****
	bool CC=1; // conv. code active (0=no,1=yes)
	int kv=1, nv=3;
	vi G={7,7,5};
	//vi G={5,7};
	fsm cc(kv,nv,G); // define FSM
	// **********************

	int k=log2(M), K=k*pow(10,ORD)*nv, /*K=k*pow(10,ORD),*/ ns=K/k, err=0, ns2=0, dchk=0, K2=0;
	double N0=0;

	// vectors
	vi EbN0dB(MAX), b_uncoded(K/nv), destu(K/nv), bs(K), sym(ns), sym_est(ns), sym_err(ns), b_est(K), b_err(K), derru(K/nv);
	vcd sym_tbl(M), s(ns), n(ns), r(ns);
	// final result vectors
	vd SER(MAX), BER(MAX), BERcc(MAX);


	// SIMULATION
	for(int z=0; z<MAX; z++) {EbN0dB[z] = z;} // orig
	sym_tbl = sym_table(MOD,M);

	for(int z=0; z<MAX; z++) {

		//b_uncoded = bitstream(randn(K/nv,0,1)); // generate bit stream

		if(!CC&&nv==1) bs=b_uncoded;
		else bs=VA_encode(cc,b_uncoded,nv);

		sym = symbol_map(bs,K,k); // assign index
		N0 = 1/(1.0*k)*pow(10,EbN0dB[z]*-1/(1.0*10))*1; // determine N0 (see note)
		n = crandn(K,0,sqrt(N0/2)); // generate noise

		for(int y=0; y<ns; y++) {
			s[y] = sym_tbl[sym[y]]; // assign symbol to tx
			r[y] = s[y] + n[y]; // rx signal (post-channel)
		}

		// ML decision [1]: min._dist=min(abs(sym_tbl-r[x]))
		err=0; ns2=0; K2=0; fill(sym_err.begin(),sym_err.end(),0);
		for (int y=0; y<ns; y++) {
			vd temp(M);
			for (int x=0; x<M; x++) {temp[x] = abs(sym_tbl[x]-r[y]);}
			sym_est[y] = distance(temp.begin(),min_element(temp.begin(),temp.end())); // outputs index of min. value
			sym_err[y] = bool(sym[y]-sym_est[y]);

			//terminate early
			ns2++;
			if(sym_err[y]) err++;
			//if(err>=100) break;
		} // end_for

		// symbol error determination [1]
		SER[z] = err/(1.0*ns2);

		// bit error determination [1]
		K2 = ns2*k;
		b_est = int2bin(sym_est,k);
		for (int y=0; y<K2; y++) {b_err[y] = bool(bs[y]-b_est[y]);}
		BER[z] = accumulate(&b_err[0],&b_err[K2-1],b_err[K2-1])/(1.0*K2);

		if(CC&&nv!=1) destu = VA_decode(cc,b_est,nv); // decode
		for (int y=0;y<destu.size(); y++) {derru[y] = bool(b_uncoded[y]-destu[y]);}
		int temp8=0;
		for (int i=0;i<derru.size();i++) {if(derru[i]) temp8++;}
		BERcc[z] = temp8/(1.0*K2);

		// specific Eb/N0 metrics
		if(!CSV) printf("Eb/N0 (dB):%2u   err:%3u   ns2:%*u   SER:%1.2e   BER:%1.2e   BERcc:%1.2e\n",EbN0dB[z],err,ORD+1,ns2,SER[z],BER[z],BERcc[z]);

		// If Es/N0 is high, then terminate early.
		//if(err||ns2!=ns) dchk=0; // dchk is a double-check checking for 2 zeros
		//if(!err&&ns2==ns&&dchk) break;
		//if(!err&&ns2==ns) dchk++;

	} // end_for

	if(!CSV) {if(K<=nv*k*100) display(EbN0dB,sym_tbl,b_uncoded,bs,sym,s,n,r,sym_est,sym_err,destu,b_est,derru,b_err);}
	else { // output in .csv format
		  for(int z=0;z<MAX;z++) {cout<<SER[z]<<",";} cout<<endl;
		  for(int z=0;z<MAX;z++) {cout<<BER[z]<<",";} cout<<endl;
		  for(int z=0;z<MAX;z++) {cout<<BERcc[z]<<",";} cout<<endl;
	}

	return 0;
} // end_main()


// REFERENCES
// [1] MPSK, by K. Bell (11/22/99)
