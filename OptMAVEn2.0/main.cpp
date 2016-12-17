#include <iostream>
#include <string>
#include "MOLECULES2.cpp"
#include "UTILITIES.h"
using namespace std;

int main()
{
	Molecule mol;
	string fn = "MoleculeK.pdb";
	mol.load(fn, "PDB");
	float rmatrix[3][3];
	float axis[3] = {0,0,1};
	float angle = 20;
	rotation_matrix3(rmatrix, axis, angle);
	char file_name[8];
	for (int i = 0; i < 18; i++)
	{
		cout << i << endl;
		mol.rotate(rmatrix);
		sprintf(file_name, "rot%u.pdb", i);
		string fn(file_name);
		mol.write(file_name, "PDB");
	}
}
