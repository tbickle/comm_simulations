#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <iterator>
#include "../include/comm.h"
using namespace std;
typedef vector<int> vi;
typedef vector<double> vd;
typedef complex<double> cd;
typedef vector<complex<double>> vcd;

// FUNCTIONS
vcd sym_table(int mod_scheme,int M)
{
	vcd v(M);
	double pi=M_PI;

	if (mod_scheme==0) {
		//if(!CSV) {cout << "********** " << M << "-ary PAM" << " **********" << endl;}
		switch (M) {
			case 2:  v={cd( -1,0), cd(  1,0)}; break;
			case 4:  v={cd( -3,0), cd( -1,0), cd(  1,0), cd( 3,0)}; break;
			case 8:  v={cd( -7,0), cd( -5,0), cd( -3,0), cd(-1,0),
				    cd(  1,0), cd(  3,0), cd(  5,0), cd( 7,0)}; break;
			case 16: v={cd(-15,0), cd(-13,0), cd(-11,0), cd(-9,0),
				    cd( -7,0), cd( -5,0), cd( -3,0), cd(-1,0),
				    cd(  1,0), cd(  3,0), cd(  5,0), cd( 7,0),
				    cd(  9,0), cd( 11,0), cd( 13,0), cd(15,0)}; break;
			default: v={}; 	break;
		} // end_switch
	}
	else if (mod_scheme==1) {
		//if(!CSV) {cout << "********** " << M << "-ary PSK" << " **********" << endl;}
		switch (M) {
			case 2:  v={cd(-1,0), cd(1,0)}; break;
			case 4:  v={cd(1,0), cd(0,1), cd(-1,0), cd(0,-1)}; break;
			case 8:  v={cd( 1,0), cd( 1/sqrt(2), 1/sqrt(2)), cd(0, 1), cd(-1/sqrt(2), 1/sqrt(2)),
				    cd(-1,0), cd(-1/sqrt(2),-1/sqrt(2)), cd(0,-1), cd( 1/sqrt(2),-1/sqrt(2))}; break;
		    	case 16: v={cd( 1, 0), cd(cos( 1*pi/8),sin( 1*pi/8)), cd(cos(2*pi/8),sin(2*pi/8)), cd(cos(3*pi/8),sin(3*pi/8)),
				    cd( 0, 1), cd(cos( 5*pi/8),sin( 5*pi/8)), cd(cos(6*pi/8),sin(6*pi/8)), cd(cos(7*pi/8),sin(7*pi/8)),
				    cd(-1, 0), cd(cos( 9*pi/8),sin( 9*pi/8)), cd(cos(10*pi/8),sin(10*pi/8)), cd(cos(11*pi/8),sin(11*pi/8)),
				    cd( 0,-1), cd(cos(13*pi/8),sin(13*pi/8)), cd(cos(14*pi/8),sin(14*pi/8)), cd(cos(15*pi/8),sin(15*pi/8))}; break;
			default: v={}; 	break;
		} // end_switch
	}
	else if (mod_scheme==2) {
		//if(!CSV) {cout << "********** " << M << "-ary QAM" << " **********" << endl;}
		switch (M) {
			case 2:  v={cd(-1,0), cd(1,0)}; break;
			case 4:  v={cd(1,0), cd(0,1), cd(-1,0), cd(0,-1)}; break;
			case 8:  v={cd( 3, 1), cd( 1, 1), cd(-1, 1), cd(-3, 1),
				    cd(-3,-1), cd(-1,-1), cd( 1,-1), cd( 3,-1)}; break;
		    	case 16: v={cd(-3, 3), cd(-1, 3), cd( 1, 3), cd( 3, 3),
				    cd( 3, 1), cd( 1, 1), cd(-1, 1), cd(-3, 1),
				    cd(-3,-1), cd(-1,-1), cd( 1,-1), cd( 3,-1),
				    cd( 3,-3), cd( 1,-3), cd(-1,-3), cd(-3,-3)}; break;
			default: v={}; 	break;
		} // end_switch
	}
	else if (mod_scheme==3) {
		//if(!CSV) {cout << M << "-ary Orthogonal" << endl;}
		//v = 1/sqrt(2^log2(M))*hadamard(M);
		v = {};
	}
	else {v = {};}

	double RMS = rms(v);
	for (int n=0; n<M; n++) {v[n]=v[n]/(1.0*RMS);}

	return v;
} // end_sym_table()

/*
vi int2bin(vi in,int k)
{
	// Note: For the purposes of converting in the comms simulation,
	//	 I don't need a binary number exceeds k (i.e. log2(M)).
	//	 Therefore, the needed vector is will be
	//	 K = k*ns = log2(M)*(#_of_symbols).
	// Note: Input = vector<int> of size ns, Output = vector<int> or vector<bool> of size K = k*ns

	int h,rem,dec,ns=in.size();
	vi out(ns*k);

	for (int n=0;n<ns;n++) {
		h=(n+1)*k-1;
		dec=in[n];
		do {rem=dec%2; out[h]=rem; dec=dec/2; h--;}
		while(dec>0); //end_dowhile
	} //end_for
	return out;
}
*/

void display(vi E,vcd st,vi du,vi dc,vi sym,vcd s,vcd n,vcd r,vi sym_est,vi sym_err,vi destu,vi destc,vi berru,vi berrc)
{
	cout << endl;
	cout<<"       Symbol Table: "; for(int n=0;n<st.size();n++){cout<< st[n]<<" ";} cout<<endl;
	//cout<<"                RMS: " << RMS <<endl;
	//cout << endl;
	//cout<<"          Signal TX: "; for (int z=0;z<s.size();z++){cout<<s[z]<<" ";} cout<<endl;
	//cout<<"              Noise: "; for (int z=0;z<n.size();z++){cout<<n[z]<<" ";} cout<<endl;
	//cout<<"          Signal RX: "; for (int z=0;z<r.size();z++){cout<<r[z]<<" ";} cout<<endl;
	cout << endl;
	cout<<"       Sym. Ind. TX: "; for(int z=0;z<sym.size();z++){cout<<sym[z]<<" ";} cout<<endl;
	cout<<"     Sym. Ind. Est.: "; for(int z=0;z<sym_est.size();z++){cout<<sym_est[z]<<" ";} cout<<endl;
	cout<<"         Sym. Error: "; for(int z=0;z<sym_err.size();z++){cout<<sym_err[z]<<" ";} cout<<endl;
	cout << endl;
	cout<<"       Data (Coded): "; for(int z=0;z<dc.size();z++){cout<<dc[z];} cout<<endl;
	cout<<"  Data Est. (Coded): "; for(int z=0;z<destc.size();z++){cout<<destc[z];} cout<<endl;
	cout<<"  Bit Error (Coded): "; for(int z=0;z<berrc.size();z++){cout<<berrc[z];} cout<<endl;
	cout << endl;
	cout<<"     Data (Uncoded): "; for(int z=0;z<du.size();z++){cout<<du[z];} cout<<endl;
	cout<<"Data Est. (Uncoded): "; for(int z=0;z<destu.size();z++){cout<<destu[z];} cout<<endl;
	cout<<"Bit Error (Uncoded): "; for(int z=0;z<berru.size();z++){cout<<berru[z];} cout<<endl;
	cout << endl;
}

double rms(vcd v)
{
	int num = v.size();
	double temp = 0.0;
	vd temp_v(num);
	
	for (int z=0; z<num; z++) {temp_v[z]=pow(abs(v[z]),2);} // t1=abs(v)^2
	temp = accumulate(&temp_v[0],&temp_v[num-1],temp_v[num-1])/(1.0*num); // t2=mean(t1)
	temp = sqrt(temp); // rms=sqrt(t2)

	return temp;
}

vd randn(int K,double mean,double std)
{
	// input: gaussian mean and standard deviation
	// output: a vector of gaussian distributed random numbers

	vd p(K);
	random_device rd;
	mt19937 e2(rd());
	normal_distribution<double> dist(mean,std);

	for (int i=0; i<K; i++) {p[i]=dist(e2);}

	return p;
}

vcd crandn(int K,double mean,double std)
{
	// input: gaussian mean and standard deviation
	// output: a complex vector of gaussian distributed random numbers

	vcd p(K);
	random_device rd;
	mt19937 e2(rd());
	normal_distribution<double> dist(mean,std);

	for (int i=0; i<K; i++) {p[i]=cd(dist(e2),dist(e2));}

	return p;
}

vi bitstream(vd in)
{
	vi out(in.size());
	for (int n=0; n<in.size(); n++) {out[n]=(copysign(1,in[n])+1)/2;}
	return out;
}

vi symbol_map(vi in,int K,int k)
{
	int ns = K/k;
	vi temp(ns);
	vi out(ns);

	switch (k) {
		case 1: out=in; break;
		case 2: for(int m=0;m<ns;m++) {out[m]=2*in[2*m]+1*in[2*m+1];} break;
		case 3: for(int m=0;m<ns;m++) {out[m]=4*in[3*m]+2*in[3*m+1]+1*in[3*m+2];} break;
		case 4: for(int m=0;m<ns;m++) {out[m]=8*in[4*m]+4*in[4*m+1]+2*in[4*m+2]+1*in[4*m+3];} break;
		default: out={}; break;
	} // end_switch
	
	return out;
}


