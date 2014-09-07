// compile: g++ -static -o t constellation_gsl.cpp -lgsl -lgslcblas -lm -std=c++11

#include <iostream>
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <random>
using namespace std;

// global vars
double pi = M_PI;

typedef gsl_vector gv;
typedef gsl_complex gc;
typedef gsl_vector_complex gvc;

// fcn prototypes
gvc* sym_table(int,int);
gc cr(double,double); // gsl_complex_rect()

void set(gv*,int,double);
void cset(gvc*,int,gc);
double get(gv*,int); 
gc cget(gvc*,int);

void vcprint(gvc*,int);
gvc* randn(int,double,double);
gv* bitstream(int);
gv* symbol_map(gv*,int,int);

int main(void)
{
	// SIMULATION PARAMETERS
	//int MOD=1; int M=4;
	int MOD; cout << "Enter Modulation Scheme (0=PAM,1=PSK,2=QAM,3=Orth): "; cin >> MOD;
	int M; cout << "Enter # of Symbols (2,4,8,16): "; cin >> M;
	const int ORD=3;
	const int MAX=20;
	int k=log2(M);
	int K=k*pow(10,ORD);
	int num_syms=K/k;
	gv* EbN0dB=gsl_vector_calloc(MAX);
	gv* SER=gsl_vector_calloc(MAX);
	gv* BER=gsl_vector_calloc(MAX);
	for(int n=0;n<MAX;n++) {set(EbN0dB,n,n);}
	//for(int n=0;n<MAX;n++) {cout<<get(EbN0dB,n);} cout<<endl; // test

	gvc* sym_tbl = sym_table(MOD,M); // M
	//vcprint(sym_tbl,M); // test

	// SIMULATION
	for(int i=0; i<MAX; i++) {
		cout<<".";
    
		// generate bit stream
		gv* b = gsl_vector_calloc(K);
		b = bitstream(K);
		//for(int n=0;n<K;n++) {cout<<get(b,n);} cout<<endl; // test

		// transmit
		gv* sym = symbol_map(b,K,k); // num_syms
		gvc* s = gsl_vector_complex_calloc(num_syms);
		for (int n=0;n<num_syms;n++) {cset(s,n,cget(sym_tbl,get(sym,n)));}
		//for (int n=0;n<num_syms;n++) {cout<<get(sym,n)<<" ";} cout<<endl; // test
		//vcprint(s,num_syms); // test

		// channel
		double N0 = pow(10,get(EbN0dB,i)*-k/10);
		//gvc* noise = gsl_vector_complex_calloc(num_syms);
		gvc* noise = randn(num_syms,0,sqrt(N0/2));
		//vcprint(noise,num_syms); // test

		// receive = signal + noise
		gvc* r = gsl_vector_complex_calloc(num_syms);
		gsl_blas_zaxpy(cr(1,0),s,noise); // noise = s+noise
		gsl_blas_zcopy(noise,r); // r = noise = noise + signal
		//vcprint(r,num_syms); // test

		// MAP decision [1]
		//for col = 1:num_syms [ee ind]=min(abs(sym_tbl-r(col)));        // find minimum distance [ee=min. dist. value, ind=index of min. dist.]
		gv* s_est = gsl_vector_calloc(num_syms);
		int nsyms = 1;
		for (int n=0;n<num_syms;n++) {
			gv* temp = gsl_vector_calloc(M);
			for (int m=0;m<M;m++) {set(temp,m,gsl_complex_abs(gsl_complex_sub(cget(sym_tbl,m),cget(r,n))));}
			int ind = gsl_vector_min_index(temp);
			set(s_est,n,ind);

			//s_est(:,col) = MAP(ind,M);
			// terminate early
			//temp = abs(symbols(:,col)-s_est(:,col));
			//if(M~=2), temp = sum(temp); }
			//if(temp>0), sym_err=sym_err+1; }
			//if (sym_err>=100), s_est=s_est(:,1:col); break; }
			gsl_vector_free(temp);
			nsyms++;
		}
		//for(int n=0;n<num_syms;n++) {cout<<get(s_est,n)<<" ";} cout<<endl; // test

		// symbol error determination [1]
		gv* sym_err = gsl_vector_calloc(num_syms); // nsyms
		gsl_vector_sub(s_est,sym);
		gsl_vector_memcpy(sym_err,s_est);
		for (int n=0;n<num_syms;n++) {if(get(sym_err,n)!=0) {set(sym_err,n,1);}}
		//for (int n=0;n<num_syms;n++) {cout<<get(sym_err,n)<<" ";} cout<<endl; // test
		set(SER,i,gsl_blas_dasum(sym_err)/num_syms);
/*
		// bit error determination [1]
		b_est = reshape(s_est,1,[]);                        // reshape symbol matrix into bit stream
		b_err_v = abs(b(1,1:k*col)-b_est);                    // a vector of bit errors (0 = no error, 1 = error)
		BER(i) = sum(b_err_v)/(k*col);
*/
		//if (BER(i) < 100/K), break; }

		gsl_vector_free(b);
		gsl_vector_free(sym);
		gsl_vector_complex_free(s);
		gsl_vector_complex_free(noise);
		gsl_vector_free(s_est);

	} cout<<endl; // end for
	for(int n=0;n<MAX;n++) {cout<<get(SER,n)<<" ";} cout<<endl;

	gsl_vector_free(EbN0dB);
	gsl_vector_free(SER);
	gsl_vector_free(BER);
	gsl_vector_complex_free(sym_tbl);
	return 0;

} // } main()




// FUNCTIONS
gvc* sym_table(int mod_scheme,int M)
{
	gvc* v = gsl_vector_complex_calloc(M);

	if (mod_scheme==0) {
		cout << M << "-ary PAM" << endl;
		switch (M) {
		case 2:  cset(v,0,cr(-1,0)); cset(v,1,cr(1,0)); break;
		case 4:  cset(v,0,cr(-3,0)); cset(v,1,cr(-1,0)); cset(v,2,cr(1,0)); cset(v,3,cr(3,0)); break;
		case 8:  cset(v,0,cr(-7,0)); cset(v,1,cr(-5,0)); cset(v,2,cr(-3,0)); cset(v,3,cr(-1,0));
			 cset(v,4,cr(1,0)); cset(v,5,cr(3,0)); cset(v,6,cr(5,0)); cset(v,7,cr(7,0)); break;
		case 16: cset(v,0,cr(-15,0)); cset(v,1,cr(-13,0)); cset(v,2,cr(-11,0)); cset(v,3,cr(-9,0));
			 cset(v,4,cr( -7,0)); cset(v,5,cr( -5,0)); cset(v,6,cr( -3,0)); cset(v,7,cr(-1,0));
			 cset(v,8,cr(  1,0)); cset(v,9,cr(  3,0)); cset(v,10,cr(  5,0)); cset(v,11,cr( 7,0));
			 cset(v,12,cr(  9,0)); cset(v,13,cr( 11,0)); cset(v,14,cr( 13,0)); cset(v,15,cr(15,0)); break;
		} // end switch
	}
	else if (mod_scheme==1) {
		cout << M << "-ary PSK" << endl;
		switch (M) {
		case 2:  cset(v,0,cr(-1,0)); cset(v,1,cr(1,0)); break;
		case 4:  cset(v,0,cr(1,0)); cset(v,1,cr(0,1)); cset(v,2,cr(-1,0)); cset(v,3,cr(0,-1)); break;
		case 8:  cset(v,0,cr(1,0)); cset(v,1,cr(1/sqrt(2),1/sqrt(2))); cset(v,2,cr(0,1)); cset(v,3,cr(-1/sqrt(2),1/sqrt(2)));
			 cset(v,4,cr(-1,0)); cset(v,5,cr(-1/sqrt(2),-1/sqrt(2))); cset(v,6,cr(0,-1)); cset(v,7,cr(1/sqrt(2),-1/sqrt(2))); break;
		case 16: cset(v,0,cr(1,0)); cset(v,1,cr(cos(1*pi/8),sin(1*pi/8))); cset(v,2,cr(cos(2*pi/8),sin(2*pi/8))); cset(v,3,cr(cos(3*pi/8),sin(3*pi/8)));
			 cset(v,4,cr(0,1)); cset(v,5,cr(cos(5*pi/8),sin(5*pi/8))); cset(v,6,cr(cos(6*pi/8),sin(6*pi/8))); cset(v,7,cr(cos(7*pi/8),sin(7*pi/8)));
			 cset(v,8,cr(-1, 0)); cset(v,9,cr(cos( 9*pi/8),sin( 9*pi/8))); cset(v,10,cr(cos(10*pi/8),sin(10*pi/8))); cset(v,11,cr(cos(11*pi/8),sin(11*pi/8)));
			 cset(v,12,cr(0,-1)); cset(v,13,cr(cos(13*pi/8),sin(13*pi/8))); cset(v,14,cr(cos(14*pi/8),sin(14*pi/8))); cset(v,15,cr(cos(15*pi/8),sin(15*pi/8)));
			 break;
		} // end switch
	}
	else if (mod_scheme==2) {
		cout << M << "-ary QAM" << endl;
		switch (M) {
			case 2:  cset(v,0,cr(-1,0)); cset(v,1,cr(1,0)); break;
			case 4:  cset(v,0,cr(1,0)); cset(v,1,cr(0,1)); cset(v,2,cr(-1,0)); cset(v,3,cr(0,-1)); break;
			case 8:  cset(v,0,cr(3,1)); cset(v,1,cr(1,1)); cset(v,2,cr(-1, 1)); cset(v,3,cr(-3, 1));
				 cset(v,4,cr(-3,-1)); cset(v,5,cr(-1,-1)); cset(v,6,cr(1,-1)); cset(v,7,cr( 3,-1)); break;
		    	case 16: cset(v,0,cr(-3,3)); cset(v,1,cr(-1,3)); cset(v,2,cr(1,3)); cset(v,3,cr(3,3));
				 cset(v,4,cr(3,1)); cset(v,5,cr(1,1)); cset(v,6,cr(-1,1)); cset(v,7,cr(-3, 1));
				 cset(v,8,cr(-3,-1)); cset(v,9,cr(-1,-1)); cset(v,10,cr(1,-1)); cset(v,11,cr( 3,-1));
				 cset(v,12,cr(3,-3)); cset(v,13,cr(1,-3)); cset(v,14,cr(-1,-3)); cset(v,15,cr(-3,-3)); break;
		} // end switch
	}
	else if (mod_scheme==3) {
		cout << M << "-ary Orthogonal" << endl;
		//v = 1/sqrt(2^log2(M))*hadamard(M);
		//v = {};
	}

	gv* temp = gsl_vector_calloc(M);
	for (int n=0;n<M;n++) {set(temp,n,gsl_complex_abs2(cget(v,n)));} // temp[n] = abs(v[n])^2
	double RMS = sqrt(gsl_blas_dasum(temp)/M); // sqrt(sum(v)/M)
	//cout << "RMS = " << RMS << endl;
	gsl_blas_zdscal(1/RMS,v);

	//double RMS = gsl_blas_dznrm2(v);
	//gsl_blas_zdscal(sqrt(pow(2,log2(M)))/RMS,v);

	return v;
	gsl_vector_complex_free(v);
	gsl_vector_free(temp);
} // end sym_table()

gv* symbol_map(gv* in,int K,int k)
{
	int num_syms = K/k;
	gv* temp = gsl_vector_calloc(num_syms);
	gv* out = gsl_vector_calloc(num_syms);

	switch (k) {
	case 1: gsl_vector_memcpy(out,in); break;
	case 2: for(int n=0;n<num_syms;n++) {
			set(out,n,2*get(in,2*n)+1*get(in,2*n+1));
		} break;
	case 3: for(int n=0;n<num_syms;n++) {
			set(out,n,4*get(in,3*n)+2*get(in,3*n+1)+1*get(in,3*n+2));
		} break;

	case 4: for(int n=0;n<num_syms;n++) {
			set(out,n,8*get(in,4*n)+4*get(in,4*n+1)+2*get(in,4*n+2)+1*get(in,4*n+3));
		} break;
	} // end switch
	
	return out;
	gsl_vector_free(temp);
	gsl_vector_free(out);
}

gc cr(double real,double imag) // gsl_complex_rect()
{
	gc out = gsl_complex_rect(real,imag);
	return out;
}

void set(gv* v,int i,double x) {gsl_vector_set(v,i,x);}
void cset(gvc* v,int i,gc x) {gsl_vector_complex_set(v,i,x);}
double get(gv* v, int i) {return gsl_vector_get(v,i);}
gc cget(gvc* v,int i) {return gsl_vector_complex_get(v,i);}

void vcprint(gvc* v,int n)
{
	for (int i=0; i<n; i++) {
		if(GSL_IMAG(cget(v,i))<0) {cout<< GSL_REAL(cget(v,i)) << GSL_IMAG(cget(v,i)) << "j ";}
		else {cout<< GSL_REAL(cget(v,i)) << "+" << GSL_IMAG(cget(v,i)) << "j ";}
	}
	cout << endl;
}	

gvc* randn(int K,double mean,double std)
{
	// input: gaussian mean and standard deviation
	// output: a vector of gaussian distributed random numbers

	gvc* out = gsl_vector_complex_calloc(K);
	random_device rd;
	mt19937 e2(rd());
	normal_distribution<double> dist(mean,std); 
	for (int n=0;n<K;n++) {cset(out,n,cr(dist(e2),dist(e2)));}
	//vcprint(out,K); // test

	return out;
	gsl_vector_complex_free(out);
}

gv* bitstream(int K)
{
	double mean=0.0; double std=1.0;
	gv* out = gsl_vector_calloc(K);
	random_device rd;
	mt19937 e2(rd());
	normal_distribution<double> dist(mean,std);
	for (int i=0; i<K; i++) {
		if(dist(e2)>0) {set(out,i,1);}
		else {set(out,i,0);}
	}
	return out;
	gsl_vector_free(out);
}

