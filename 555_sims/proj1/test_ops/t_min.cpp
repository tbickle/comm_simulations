// compile: g++ -o t t_min.cpp -std=c++11
// min_element/max_element example
#include <iostream>     // std::cout
#include <algorithm>    // std::min_element, std::max_element
#include <iterator>
#include <vector>
using namespace std;

//bool myfn(int i,int j) {return i<j;}
//struct myclass {bool operator() (int i,int j) {return i<j;}} myobj;

int main () {

	vector<int> v(7);
	v = {3,7,5,2,6,4,9};

	for(int n=0;n<v.size();n++) {cout << v[n];} cout<<endl;
	cout<<"smallest element: "<< *min_element(&v[0],&v[v.size()]) <<'\n';
	double acc = accumulate(v.begin(),v.end(),0);	//double test = accumulate(&v[0],&v[v.size()],0);
	cout << "accumulated: " << acc << endl; // total = 36

	int min_pos = distance(v.begin(),min_element(v.begin(),v.end()));
	cout << "index of minimum value: " << min_pos << endl;

	cout << bool(4) << endl;
	cout << bool(4.1) << endl;
	cout << bool(-1.1) << endl;
	cout << bool(0) << endl;

	// PROVIDED TEST
	//int myints[] = {3,7,2,5,6,4,9};

	// using default comparison:
	//cout<<"The smallest element is "<< *min_element(myints,myints+7) <<'\n';
	//cout<<"The largest element is "<< *max_element(myints,myints+7) <<'\n';

	// using function myfn as comp:
	//std::cout << "The smallest element is " << *std::min_element(myints,myints+7,myfn) << '\n';
	//std::cout << "The largest element is "  << *std::max_element(myints,myints+7,myfn) << '\n';

	// using object myobj as comp:
	//std::cout << "The smallest element is " << *std::min_element(myints,myints+7,myobj) << '\n';
	//std::cout << "The largest element is "  << *std::max_element(myints,myints+7,myobj) << '\n';

	return 0;
}
