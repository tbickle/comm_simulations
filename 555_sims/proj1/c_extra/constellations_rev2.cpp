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

// global vars
const bool CSV=1; // simulation outputs in .csv format (0=No,1=Yes)
double RMS=1;

// fcn prototypes
vcd sym_table(int,int);
double rms(vcd);
vd randn(int,double,double);
vi bitstream(vd);
vi symbol_map(vi,int,int);
void display(vi,vcd,vi,vi,vcd,vcd,vcd,vi,vi,vi,vi);
vcd crandn(int K,double mean,double std);
vi int2bin(vi,int);

int main(int argc, char *argv[])
{
	if(!CSV) cout<<endl;
	// INITIALIZATION
	//int MOD; cout << "Enter Modulation Scheme (0=PAM,1=PSK,2=QAM,3=Orth): "; cin >> MOD;
	//int M; cout << "Enter # of Symbols (2,4,8,16): "; cin >> M;
	int MOD=atoi(argv[1]), M=atoi(argv[2]);
	const int ORD=5, MAX=35;
	int k=log2(M), K=k*pow(10,ORD), ns=K/k, err=0, ns2=0, dchk=0, K2=0;
	double N0=0;
	// vectors
	vi EbN0dB(MAX), bs(K), sym(ns), sym_est(ns), sym_err(ns), b_est(K), b_err(K);
	vcd sym_tbl(M), s(ns), n(ns), r(ns);
	// final result vectors
	vd SER(MAX), BER(MAX);


	// SIMULATION
	for(int z=0; z<MAX; z++) {EbN0dB[z] = z;} // orig
	sym_tbl = sym_table(MOD,M);

	for(int z=0; z<MAX; z++) {

		bs = bitstream(randn(K,0,1)); // generate bit stream
		sym = symbol_map(bs,K,k); // assign index
		N0 = 1/(1.0*k)*pow(10,EbN0dB[z]*-1/(1.0*10))*1; // determine N0 (see note)
		n = crandn(K,0,sqrt(N0/2)); // generate noise

		for(int y=0; y<ns; y++) {
			s[y] = sym_tbl[sym[y]]; // assign symbol to tx
			r[y] = s[y] + n[y]; // rx signal (post-channel)
		}

		// ML decision [1]: min._dist=min(abs(sym_tbl-r[x]))
		err=0; ns2=0; fill(sym_err.begin(),sym_err.end(),0); K2=0;
		for (int y=0; y<ns; y++) {
			vd temp(M);
			for (int x=0; x<M; x++) {temp[x] = abs(sym_tbl[x]-r[y]);}
			sym_est[y] = distance(temp.begin(),min_element(temp.begin(),temp.end())); // outputs index of min. value
			sym_err[y] = bool(sym[y]-sym_est[y]);

			//terminate early
			ns2++;
			if(sym_err[y]) err++;
			if(err>=100) break;
		} // end_for

		// symbol error determination [1]
		SER[z] = err/(1.0*ns2);

		// bit error determination [1]
		K2 = ns2*k;
		//b_est.resize(K2), b_err.resize(K2);
		b_est = int2bin(sym_est,k);
		for (int y=0; y<K2; y++) {b_err[y] = bool(bs[y]-b_est[y]);}
		BER[z] = accumulate(&b_err[0],&b_err[K2-1],b_err[K2-1])/(1.0*K2);

		// specific Eb/N0 metrics
		if(!CSV) {cout<<"Eb/N0 (dB): "<<EbN0dB[z]<<"  err: "<<err<<"  ns2: "<<ns2<<"  SER: "<<SER[z]<<"  BER: "<<BER[z]<<endl;}

		// If Es/N0 is high, then terminate early.
		if(err||ns2!=ns) dchk=0; // dchk is a double-check checking for 2 zeros
		if(!err&&ns2==ns&&dchk) break;
		if(!err&&ns2==ns) dchk++;
	} // end_for

	if(!CSV) {if(K<=k*100) display(EbN0dB,sym_tbl,bs,sym,s,n,r,sym_est,sym_err,b_est,b_err);}
	else { // output in .csv format
		  for(int z=0;z<MAX;z++) {cout<<SER[z]<<",";} cout<<endl;
		  for(int z=0;z<MAX;z++) {cout<<BER[z]<<",";} cout<<endl;
	}

	return 0;
} // end_main()


// FUNCTIONS
vcd sym_table(int mod_scheme,int M)
{
	vcd v(M);
	double pi=M_PI;

	if (mod_scheme==0) {
		if(!CSV) {cout << "********** " << M << "-ary PAM" << " **********" << endl;}
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
		if(!CSV) {cout << "********** " << M << "-ary PSK" << " **********" << endl;}
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
		if(!CSV) {cout << "********** " << M << "-ary QAM" << " **********" << endl;}
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
		if(!CSV) {cout << M << "-ary Orthogonal" << endl;}
		//v = 1/sqrt(2^log2(M))*hadamard(M);
		v = {};
	}
	else {v = {};}

	RMS = rms(v);
	for (int n=0; n<M; n++) {v[n]=v[n]/(1.0*RMS);}

	return v;
} // end_sym_table()

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

void display(vi E,vcd st,vi bs,vi sym,vcd s,vcd n,vcd r,vi sym_est,vi sym_err,vi b_est,vi b_err)
{
	cout << endl;
	cout<<"  Symbol Table: "; for(int n=0;n<st.size();n++){cout<< st[n]<<" ";} cout<<endl;
	cout<<"           RMS: " << RMS <<endl;
	cout << endl;
	cout<<"     Signal TX: "; for (int z=0;z<s.size();z++){cout<<s[z]<<" ";} cout<<endl;
	cout<<"         Noise: "; for (int z=0;z<s.size();z++){cout<<n[z]<<" ";} cout<<endl;
	cout<<"     Signal RX: "; for (int z=0;z<s.size();z++){cout<<r[z]<<" ";} cout<<endl;
	cout << endl;
	cout<<"  Sym. Ind. TX: "; for(int z=0;z<sym.size();z++){cout<<sym[z]<<" ";} cout<<endl;
	cout<<"Sym. Ind. Est.: "; for(int z=0;z<sym_est.size();z++){cout<<sym_est[z]<<" ";} cout<<endl;
	cout<<"    Sym. Error: "; for(int z=0;z<sym_err.size();z++){cout<<sym_err[z]<<" ";} cout<<endl;
	cout << endl;
	cout<<"    Bit Stream: "; for(int z=0;z<bs.size();z++){cout<<bs[z];} cout<<endl;
	cout<<"      Bit Est.: "; for(int z=0;z<b_est.size();z++){cout<<b_est[z];} cout<<endl;
	cout<<"     Bit Error: "; for(int z=0;z<b_err.size();z++){cout<<b_err[z];} cout<<endl;
	cout << endl;
}


// REFERENCES
// [1] MPSK, by K. Bell (11/22/99)
