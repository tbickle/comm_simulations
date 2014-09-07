// passing vectors in/out of fuctions test

#include <iostream>
#include <vector>

using namespace std;

vector<int> function1(void);

int main()
{
	vector<int> v;
	v = function1();

	cout<<"size:"<<v.size()<<endl;
	cout<<"contents: ";
	for(int i=0; i<v.size(); i++) {cout<<v[i]<<" ";}
	cout<<endl;


	return 0;
}

vector<int> function1(void)
{
	vector<int> v(20);
	for(int i=0; i<v.size(); i++) {v[i]=i*3;}
	return v;
}
