#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;

#ifndef MOLVEC
#define MOLVEC

class Point3D
{
public:
	float coords[3];
	Point3D(float x, float y, float z)
	{
		coords[0] = x; coords[1] = y; coords[2] = z;
	};
};


class MolVec
{
private:
	int npoints;
public:
	MolVec() {};
	vector<Point3D> coords;
	void load(string file_name)
	{
		coords.clear();
		npoints = 0;
		ifstream in(file_name.c_str());
		string line;
		float x, y, z;
		getline(in, line);
		while (!in.eof())
		{
			stringstream data(line);
			data >> x >> y >> z;
			coords.push_back(Point3D(x, y, z));
			npoints++;
			getline(in, line);
		}
	}
	int get_size()
	{
		return npoints;
	}
	void print()
	{
		for (vector<Point3D>::iterator it = coords.begin(); it != coords.end(); it++)	
		{
			cout << "(";
			for (int i = 0; i < 3; i++)
			{
				cout << (*it).coords[i];
				if (i < 2)
				{
					cout << ", ";
				}
			}
			cout << ")" << endl;
		}
	}
};

#endif
