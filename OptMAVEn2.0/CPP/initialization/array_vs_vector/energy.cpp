#include <fstream>
#include "molvec.cpp"
#include <cmath>
using namespace std;


float interactionEnergy(MolVec mv1, MolVec mv2)
{
	float energy = 0.0;
	float dist2;
	int i, j, k;
	for (i = 0; i < mv1.get_size(); i++)
	{
		for (j = 0; j < mv2.get_size(); j++)
		{
			dist2 = 0.0;
			for (k = 0; k < 3; k++)
			{
				dist2 += pow(mv1.coords.at(i).coords[k] - mv2.coords.at(j).coords[k], 2);
			}
			energy += pow(dist2, -0.5);
		}
		if (i % 1000 == 0)
		{
			cout << i << endl;
		}
	}
	return energy;
}
