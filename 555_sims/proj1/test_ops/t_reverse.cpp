#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

int main(void)
{
	vector<int> v(3);
	v = {1,2,3};
	reverse(v.begin(),v.end());
	for(int i=0;i<v.size();i++) cout<<v[i]; cout<<endl;
	return 0;
}
