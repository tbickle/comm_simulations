// compile: g++ -o t constellations.cpp -std=c++11
// note: I have a scaling issue. My plots are trending correctly, but they are off from what they should be.
// update: scaling issue seems to have been resolved by scaling N0 by 1/sqrt(2^k)

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
double pi = M_PI;
double RMS = 1;

// fcn prototypes
vcd sym_table(int,int);
double rms(vcd);
vd randn(int,double,double);
vi bitstream(vd);
vi symbol_map(vi,int,int);
void display(vi,vcd,vi,vi,vcd,vcd,vcd,vi,vi);
vcd crandn(int K,double mean,double std);

int main(int argc, char *argv[])
{
	// INITIALIZATION
	cout<<endl;
	//int MOD; cout << "Enter Modulation Scheme (0=PAM,1=PSK,2=QAM,3=Orth): "; cin >> MOD;
	//int M; cout << "Enter # of Symbols (2,4,8,16): "; cin >> M;
	int MOD = atoi(argv[1]); int M = atoi(argv[2]);

	const int ORD = 5;
	const int MAX = 35;
	int k = log2(M);
	int K = k*pow(10,ORD);
	int ns = K/k;
	bool print = 1;
	int err = 0;
	int ns2 = 0;
	int dchk = 0;
	// vectors
	vi EbN0dB(MAX);
	vcd sym_tbl(M);
	vi bs(K);
	vi sym(ns);
	vcd s(ns); // signal
	double N0 = 0;
	vcd n(ns); // noise
	vcd r(ns); // received signal
	vi sym_est(ns); // signal estimation
	vi sym_err(ns);
	// final result vectors
	vd SER(MAX);
	vd BER(MAX);


	// SIMULATION
	for(int z=0; z<MAX; z++) {EbN0dB[z] = z;} // orig
	sym_tbl = sym_table(MOD,M);

	for(int z=0; z<MAX; z++) {

		// generate bit stream
		bs = bitstream(randn(K,0,1));

		// transmit
		sym = symbol_map(bs,K,k); // assign index
		for(int y=0; y<ns; y++) {s[y] = sym_tbl[sym[y]];} // assign symbol to tx

		// channel
		//N0 = pow(10,EbN0dB[z]*-1/(1.0*10))*1; // orig (EbN0dB here is actually EsN0dB)
		N0 = 1/sqrt(pow(2,k))*pow(10,EbN0dB[z]*-1/(1.0*10))*1; // orig2 EbN0dB here is actually EsN0dB <- This seems to work!

		n = crandn(K,0,sqrt(N0/2));

		// receive: r = s + n
		for(int y=0; y<ns; y++) {r[y] = s[y] + n[y];}

		// MAP decision [1]: min. dist = min(abs(sym_tbl-r[x]))
		err = 0; ns2 = 0;
		fill(sym_err.begin(),sym_err.end(),0); // new
		for (int y=0; y<ns; y++) {
			vd temp(M); // original
			//int test=0; // test
			for (int x=0; x<M; x++) {temp[x] = abs(sym_tbl[x]-r[y]);} // original
			//for(int x=0;x<M;x++){temp[x]=abs(sym_tbl[x]-r[y]); if(bool(x)){if(temp[x]<=temp[x-1]) test=x;}} // test
			sym_est[y] = distance(temp.begin(),min_element(temp.begin(),temp.end())); // original
			//sym_est[y] = test; // test
			sym_err[y] = bool(sym[y]-sym_est[y]); // original

			//terminate early
			ns2++;
			if(sym_err[y]) err++;
			if(err>=100) break;
		} // end_for

		// symbol error determination [1]
		//SER[z] = accumulate(&sym_err[0],&sym_err[ns-1],sym_err[ns-1])/(1.0*ns2);
		SER[z] = err/(1.0*ns2);
		cout << "Eb/N0 (dB): " << EbN0dB[z] << "  err: " << err << "  ns2: " << ns2 << "  SER: " << SER[z] << endl;

		// bit error determination [1]
		//BER(i) = sum(b_err_v)/(k*col);

		// If Es/N0 is high, then terminate early.
		if(err||ns2!=ns) dchk=0; // dchk is a double-check checking for 2 zeros
		if(!err&&ns2==ns&&dchk) break;
		if(!err&&ns2==ns) dchk++;
	} // end for

	if(print&&K<=k*100) display(EbN0dB,sym_tbl,bs,sym,s,n,r,sym_est,sym_err);
	//cout<<"SER: "; for(int z=0;z<MAX;z++) {cout<<SER[z]<<" ";} cout<<endl<<endl;

	return 0;
} // } main()


// FUNCTIONS
vcd sym_table(int mod_scheme,int M)
{
	vcd v(M);

	if (mod_scheme==0) {
		cout << "********** " << M << "-ary PAM" << " **********" << endl;
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
		} // end switch
	}
	else if (mod_scheme==1) {
		cout << "********** " << M << "-ary PSK" << " **********" << endl;
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
		} // end switch
	}
	else if (mod_scheme==2) {
		cout << "********** " << M << "-ary QAM" << " **********" << endl;
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
		} // end switch
	}
	else if (mod_scheme==3) {
		cout << M << "-ary Orthogonal" << endl;
		//v = 1/sqrt(2^log2(M))*hadamard(M);
		v = {};
	}
	else {v = {};}

	RMS = rms(v); // orig
	for (int n=0; n<M; n++) {v[n]=v[n]/(1.0*RMS);}

	return v;
} // end sym_table()

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
	} // end switch
	
	return out;
}

void display(vi E,vcd st,vi bs,vi sym,vcd s,vcd n,vcd r,vi sym_est,vi sym_err)
{
	cout << endl;
	//cout<<"    Eb/N0 (dB): "; for(int n=0;n<E.size();n++){cout<<E[n]<<" ";} cout<<endl;
	cout<<"  Symbol Table: "; for(int n=0;n<st.size();n++){cout<< st[n]<<" ";} cout<<endl;
	cout<<"           RMS: " << RMS <<endl;
	cout<<"    Bit Stream: "; for(int z=0;z<bs.size();z++){cout<<bs[z];} cout<<endl;
	//cout<<"Symbol Indexes: "; for(int z=0;z<sym.size();z++){cout<<sym[z]<<" ";} cout<<endl;
	cout << endl;
	//cout<<"     Signal TX: "; for (int z=0;z<s.size();z++){cout<<s[z]<<" ";} cout<<endl;
	//cout<<"         Noise: "; for (int z=0;z<s.size();z++){cout<<n[z]<<" ";} cout<<endl;
	//cout<<"     Signal RX: "; for (int z=0;z<s.size();z++){cout<<r[z]<<" ";} cout<<endl;
	//cout << endl;
	cout<<"  Sym. Ind. TX: "; for(int z=0;z<sym.size();z++){cout<<sym[z]<<" ";} cout<<endl;
	cout<<"Sym. Ind. Est.: "; for(int z=0;z<sym_est.size();z++){cout<<sym_est[z]<<" ";} cout<<endl;
	cout<<"    Sym. Error: "; for(int z=0;z<sym_err.size();z++){cout<<sym_err[z]<<" ";} cout<<endl;
	cout << endl;
}

