#include <iostream>
#include <vector>
using namespace std;

// Figure out how to create a vector.
int main()
{
	vector<int[3]*> ints;
	for (int i = 0; i < 10; i++)
	{
		int* item[3] = {i, i*i, i*i*i};
		//item[0] = i; item[1] = i*i; item[2] = i*i*i;
		ints.push_back(item);
	}
	/*
	for (vector<int[3]>::iterator it = ints.begin(); it != ints.end(); it++)
	{
		for (int j = 0; j < 3; j++)
		{
			//cout << (*it) << ", ";
		}
		cout << endl;
	}
	*/

	return 0;
}
