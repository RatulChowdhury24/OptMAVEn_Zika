#include "UTILITIES2.h"
#include <iostream>
using namespace std;

int main()
{
vector<float> x;
int i;
float nums[6] = {3.1, 2.9, 0.1, -9.8, 10.3, 0.0};
for (i = 0; i < 6; i++)
	x.push_back(nums[i]);
vector<int> y = argsort(x);
// Should print 3 5 2 1 0 4
for (i = 0; i < y.size(); i++)
	cout << y[i] << ", ";
cout << endl;
}
