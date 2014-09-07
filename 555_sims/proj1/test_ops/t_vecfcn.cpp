// passing vectors in/out of fuctions test

#include <iostream>
#include <vector>

using namespace std;

vector<int> function1(vector<int> v);

int main()
{
	vector<int> test;
	vector<int> test2;
	test2 = function1(test);

	cout<<"size:"<<test2.size()<<endl;
	cout<<"contents: ";
	for(int i = 0; i < test2.size(); i ++)
	{
		cout<<test2[i]<<" ";
	}
	cout<<endl;
	return 0;
}

vector<int> function1(vector<int> v)
{
	v.resize(20);
	for(int i = 0; i < 20; i++)
	{
		v[i] = i * 3;
	}

	return v;
}
