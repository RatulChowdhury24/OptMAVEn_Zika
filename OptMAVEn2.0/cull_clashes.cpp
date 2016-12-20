/* This program performs a grid search to find locations of the antigen in which it does not clash with the antibody.
*/

#include <fstream>
#include <utility>
#include <vector>
#include "ENERGY2.h"
#include "MOLECULES2.h"
#include "STANDARDS2.h"


vector<array<float, 3>, 2> test_all(Molecule mol1, Molecule mol2, vector<pair<float, array<float, 3>>> movements)
{
	// Move the antigen to the first position to be tested.
	array<float, 3> translation = {x_trans[0], y_trans[0], z_trans[0]};
	array<float, 3> negz = {0.0, 0.0, -1.0};
	antigen.rotate(negz, z_rot[0]);
	antigen.translate(translation);
	// Assume that atoms in mol1 that are closest to the center of mol2 are most likely to clash, and vice versa, so search them for clashes first. Hope that this heuristic reduces pairwise clash-finding time from O(N^2) to O(NlogN), i.e. the sorting time. Begin by sorting all atoms based on their distance to the center of the other Molecule.
	//cout << endl << "sorting ";
	array<float, 3> center1 = mol1.get_center(); array<float, 3> center2 = mol2.get_center();  // Locate the center of each molecule.
	// Compute the distance from each atom in each Molecule to the center of the other Molecule.
	int i; vector<float> dist1(mol1.size()); vector<float> dist2(mol2.size());
	i = 0; generate(dist1.begin(), dist1.end(), [&]{return dist3(mol1.c[i++], center2);});
	i = 0; generate(dist2.begin(), dist2.end(), [&]{return dist3(mol2.c[i++], center1);});
	vector<int> order1 = argsort(dist1);  // Sort atoms in mol1 by distance.
	vector<int> order2 = argsort(dist2);  // Sort atoms in mol2 by distance.
}


int main(int argc, char* argv[])
{
	//string id = string(argv[1]);  // Get the ID of this process.
	//string file_path = 
	// FIXME: make this part customizable
	vector<float> z_rot, x_trans, y_trans, z_trans;
	float v;
	for (v = 0.0; v < 360.0; v += 20.0)
		z_rot.push_back(v);
	for (v = -5.0; v <= 10.0; v += 1.25)
		x_trans.push_back(v);
		y_trans.push_back(v - 5.0);
	for (v = -16.25; v <= -3.75; v += 1.25)
		z_trans.push_back(v);
	// Determine all possible movements to test.
	vector<pair<float, array<float, 3>>> movements;
	for (vector<float>::iterator itzr = z_rot.begin(); itzr != z_rot.end(); itzr++)
	{
		for (vector<float>::iterator itxt = x_trans.begin(); itxt != x_trans.end(); itxt++)
		{
			for (vector<float>::iterator ityt = y_trans.begin(); ityt != y_trans.end(); ityt++)
			{
				for (vector<float>::iterator itzt = z_trans.begin(); itzt != z_trans.end(); itzt++)
				{
					pair<float, array<float, 3>> movement; movement.first = *itzr;
					array<float, 3> trans = {*itxt, *ityt, *itzt}; movement.second = trans;
					movements.push_back(movement);
				}
			}
		}
	}
	// Load the Molecules.
	Antigen antigen; antigen.load(AgInit());
	Molecule IgH; IgH.load("MoleculeH.pdb"); // FIXME: make names self-consistent!
	Molecule IgK; IgK.load("MoleculeK.pdb");
	//cout << mol1.c.size() << " " << mol2.c.size() << endl;
	
	
	
}
