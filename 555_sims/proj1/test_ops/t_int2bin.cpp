#include <iostream>
#include <vector>
using namespace std;

typedef vector<int> vi;
vi int2bin(vi,int);

int main(void)
{
	// init
	int ns = 16; // # of symbols at input (i.e. channel use)
	int k = 4; // bits/symbol
	vi in(ns);
	vi out(ns*k);
	for(int n=0;n<ns;n++) {in[n]=n;} // vector of symbol decimal values, test range 0-15

	// calc
	out = int2bin(in,k);

	// disp
	cout<<"integer: "; for(int n=0;n<ns;n++) {cout<<in[n]<<" ";} cout<<endl;
	cout<<" binary: "; for(int n=0;n<ns;n++) {for(int m=0;m<k;m++) {cout<<out[n*k+m];} cout<<" ";} cout<<endl;
	return 0;
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

/*
// ALTERNATIVE
#include <iostream>
#include <algorithm>
using namespace std;

void binary(int);

int main(int argc, char *argv[])
{
	int number = atoi(argv[1]);
	//int number;
	//cout << "Please enter a positive integer: "; cin >> number;
	//if (number < 0) {cout<<"That is not a positive integer.\n";}
	//else {cout<<number<<" converted to binary is: "; binary(number); cout<<endl;}
	binary(number);
	return 0;
}

void binary(int number) {
	int remainder;
	int inif=0;
	int outif=0;
	cout<<"int: "<<number<<endl;

	if(number <= 1) {
		inif++;
		cout << "in if: " << inif << endl;
		cout << number;
		return;
	}

	outif++;
	cout << "out if: " << outif << endl;
	remainder = number%2;
	binary(number >> 1);    
	cout<<remainder<<endl; // not the full output
}
*/


