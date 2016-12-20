/* This file contains functions for computing energies and finding clashes in OptMAVEn 2.0.
*/

#ifndef ENERGY_H
#define ENERGY_H

#include <cmath>
#include "MOLECULES2.h"

// Check for a clash between two Molecules.
bool clashing(Molecule mol1, vector<int> order1, Molecule mol2, vector<int> order2)
{
	int i;	
	// Begin by sorting all atoms based on their distance to the center of the other Molecule.
	//cout << endl << "sorting ";
	array<float, 3> center1 = mol1.get_center(); array<float, 3> center2 = mol2.get_center();  // Locate the center of each molecule.
	vector<float> dist1(mol1.size()); vector<float> dist2(mol2.size());
	i = 0; generate(dist1.begin(), dist1.end(), [&]{return dist3(mol1.c[i++], center2);});  // Compute distances between atoms in mol1 and center2.
	i = 0; generate(dist2.begin(), dist2.end(), [&]{return dist3(mol2.c[i++], center1);});  // Compute distances between atoms in mol2 and center1.
	vector<int> closest1 = argsort(dist1);  // Sort atoms in mol1 by distance.
	vector<int> closest2 = argsort(dist2);  // Sort atoms in mol1 by distance.
	//cout << mol1.c.size() << " " << mol2.c.size() << endl;
	// Assume that atoms in mol1 that are closest to the center of mol2 are most likely to clash, and vice versa, so search them for clashes first. Hope that this heuristic reduces pairwise clash-finding time from O(N^2) to O(NlogN), i.e. the sorting time.
	//cout << "searching ";
	//i = 0;
	for (vector<int>::iterator i1 = closest1.begin(); i1 != closest1.end(); i1++)
	{
		for (vector<int>::iterator i2 = closest2.begin(); i2 != closest2.end(); i2++)
		{
			// Compare squared distances to avoid the expense of calculating the square root.
			//cout << *i1 << " " << *i2 << " (" << mol1.c[*i1][0] << "," << mol1.c[*i1][1] << "," << mol1.c[*i1][2] << "," << mol2.c[*i2][0] << "," << mol2.c[*i2][1] << "," << mol2.c[*i2][2] << ") " << dist3(mol1.c[*i1], mol2.c[*i2], false) << endl;
			//char dummy;
			//cin >> dummy;
			//i++;
			if (dist3(mol1.c[*i1], mol2.c[*i2], true) < 1.0) // FIXME: base clash distance on square of atomic radii.
			{
				// Because a clash was found between these two atoms, it is likely that clashes will be found between either of these atoms in other translations, so move them to the beginning of the search list.
				iter_swap(closest1.begin(), i1);
				iter_swap(closest2.begin(), i2);
				return true;  // A clash was found.
			}
		}
	}
	//cout << i << endl;
	return false;  // No clashes were found.
}

#endif
