#include <iostream>
#include <vector>
//#include "base.h"
#include "fsm.h"
using namespace std;
using namespace gr;
using namespace trellis;

typedef vector<int> vi;

int main()
{
	//vector<vi> G;
	/*
	int row = 3; // # of rows
	int col = 3; // # of columns
	G.resize(row,vi(col));

	for (int a=0; a<row; a++) {
		for (int b=0; b<col; b++) {cout<< G[a][b] << " ";}
		cout<<endl;
	}
	*/
	vi G(3);
	G = {7,7,5};
	fsm(1,3,G);
	return 0;
}
