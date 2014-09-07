#include <vector>
#include "../include/fsm.h"
#include "../include/viterbi.h"
using namespace gr;
using namespace trellis;

vector<int> VA_encode(gr::trellis::fsm cc,vector<int> data,int n)
{
	int state, in, next;
	vector<int> out(data.size());
	vector<int> outb(data.size()*n);
	vector<int> ns = cc.NS();
	vector<int> os = cc.OS();
	state=0;				// inital state is 0
	for(int i=0;i<data.size();i++) {
		in=data[i]; out[i]=os[state*cc.I()+in]; next=ns[state*cc.I()+in];
		//printf("state:%d, in:%d -> out:%d, next_state:%d\n",state,in,out[i],next);
		state=next;
	}
	outb = int2bin(out,n);
	return outb;
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
