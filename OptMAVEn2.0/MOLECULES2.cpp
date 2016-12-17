/* This file is the source for the Molecule object of OptMAVEn 2.0
*/

// Include MOLECULES.cpp only once.
#ifndef MOLECULES_CPP
#define MOLECULES_CPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include "MOLECULES2.h"
#include "UTILITIES.h"
using namespace std;


// Ensure that the x, y, and z vectors have the same length.
bool Molecule::checkXYZ()
{
	return ((c[0].size() == c[1].size()) && (c[1].size() == c[2].size()));
}

// Get the number of atoms in the Molecule.
int Molecule::size()
{
	if (checkXYZ())
	{
		return c[0].size();
	}
	else
	{
		return -1;
	}
}

// Load molecular data from a file.
bool Molecule::load(string file_name, string format)
{
	if (format == "PDB")
	{
		// Load from a PDB file.
		PDB parser;
		//	return parser.parse(file_name, chain, ai, an, ri, rn, c);
		return parser.parse(file_name, *this);
	}
	else
	{
		// Cannot load from an unrecognized file format.
		return false;
	}
}

// Print the molecular data.
bool Molecule::write(string file_name, string format)
{
	// Ensure the coordinates are consistent.
	if (!checkXYZ())
		return false;
	if (format == "PDB")
	{
		// Write a PDB file.
		PDB writer;
		writer.write(file_name, *this);
	}
	else
	{
		// Cannot write to an unrecognized file format.
		return false;
	}
	return true;
}

// Translate the molecule using a translation vector.
void Molecule::translate(float trans[3])
{
	int i, j;  // Declate index variables.
	// Translate all points.
	for (i = 0; i < size(); i++)
	{
		for (j = 0; j < 3; j++)
			c[j][i] += trans[j];
	}
}

// Rotate the molecule using a rotation matrix.
void Molecule::rotate(float mat[3][3])
{
	int i, j, k;  // Declare index variables for the matrix multiplication.
	// Declare an array to hold the old coordinates of the point.
	float old[3];
	// Rotate all points.
	for (i = 0; i < size(); i++)
	{
		// Store the old coordinates of the point.
		for (j = 0; j < 3; j++)
			old[j] = c[j][i];
		// Multiply the old coordinates and the rotation matrix.
		for (j = 0; j < 3; j++)
		{
			c[j][i] = 0.0;  // Clear the old coordinate.
			for (k = 0; k < 3; k++)
			{
				c[j][i] += old[k] * mat[k][j];
			}
		}
	}
}

// (Overloaded) Calculate the center of the Molecule or a subset of its residues.
void Molecule::get_center(float coords[3])
{
	int i, j;  // Declare index variables.
	// Perform for all three dimensions.
	for (j = 0; j < 3; j++)
		coords[j] = 0.0;  // Clear the coordinates.
		// Sum the coordinates in that dimension.
		for (i = 0; i < size(); i++)
			coords[j] += c[j][i];
		// Average the coordinates.
		coords[j] /= size();
}

void Molecule::get_center(float coords[3], vector<int> residues)
{
	cout << "Residues: ";
	for (int i = 0; i < residues.size(); i++)
		cout << residues[i] << ", ";
	cout << endl;
	int i, j;  // Declare index variables.
	// Make a set from the residues for faster searching.
	unordered_set<int> resiset;
	resiset.insert(residues.begin(), residues.end());
	// Clear the data.
	int atoms_selected = 0;
	for (j = 0; j < 3; j++)
		coords[j] = 0.0;
	// Examine each atom in the Molecule.
	for (i = 0; i < size(); i++)
	{
		// Test for membership in the selection.
		if (resiset.find(ri[i]) != resiset.end())
		{
			atoms_selected++;  // Count the atom.
			// Add the coordinates of the atom.
			for (j = 0; j < 3; j++)
			{
				coords[j] += c[j][i];
			}
		}
	}
	// Average the coordinates.
	for (j = 0; j < 3; j++)
	{
		coords[j] / atoms_selected;
	}
}

// Center the Molecule on the origin.
void Molecule::translate_center()
{
	float trans[3];
	get_center(trans);
	for (int i = 0; i < 3; i++)
		trans[i] *= -1;
	translate(trans);
}

// Translate the Molecule so that its geometric center lies at the given point.
void Molecule::translate_center(float point[3])
{
	float trans[3];
	get_center(trans);
	// The translation vector is the desired center point minus the current center point.
	for (int i = 0; i < 3; i++)
		trans[i] = point[i] - trans[i];
	translate(trans);
}

// Select a subset of residues to be the epitope.
void Antigen::set_epitope(vector<int> res)
{
	epitope.clear();  // Erase any previous epitope residues.
	// Copy all residue numbers into the epitope.
	for (vector<int>::iterator it = res.begin(); it != res.end(); it++)
		epitope.push_back(*it);
}

// Calculate the center of geometry of the epitope.
void Antigen::epicenter(float coords[3])
{
	get_center(coords, epitope);
}

// Rotate the Antigen so that the z coordinates of its epitope are minimized.
void Antigen::epitope_zmin()
{
	cout << "Centering" <<endl;
	translate_center();  // Move the Antigen to the origin so that the center of geometry does not compromise the rotation.
	// Get the epitope's center of geometry.
	float epivec[3];
	epicenter(epivec);
	// Compute the angle between epivec and the negative z axis (equal to -k, where k is the unit vector <0, 0, 1>) using the formula angle = acos(a * b / (|a||b|)), where a is epivec and b is -k.
	cout << "Axis and angle" << endl;
	float angle = acos(-epivec[2] / norm3(epivec));  // a * b = epivec[0] * 0 + epivec[1] * 0 + epivec[2] * -1; |a| = epinorm and |b| = 1.
	// Compute the axis of rotation, which is the cross product of epivec and -k, or a x b / |a x b|.
	float axis[3] = {-epivec[1], epivec[0], 0.0};  // a x b = <epivec[1] * -1 - epivec[2] * 0, epivec[2] * 0 - epivec[0] * -1, epivec[0] * 0 - epivec[1] * 0>
	// Compute the rotation matrix using the axis and angle.
	float matrix[3][3];
	cout << "Generate matrix" << endl;
	rotation_matrix3(matrix, axis, angle);
	cout << "Rotating" << endl;
	rotate(matrix);  // Rotate the Antigen.	
}

// Parse a PDB file.
bool PDB::parse(string file_name, Molecule& mol)
{
	ifstream in(file_name.c_str());  // Open the file as a filestream.
	// Declare containers for data that will be loaded.
	string line;
	float c[3];
	int ai, ri;
	string chain, type, an, rn;
	int i;
	getline(in, line);  // Read the first line.
	// Read until the end of the file.
	while (!in.eof())
	{
		stringstream lstream(line);  // Make a stringstream from the line.
		lstream >> type;  // Read the first datum in the line: the type.
		// Only read data from ATOMs and HETATMs.
		if (type == "ATOM" || type == "HETATM")
		{
			lstream >> ai >> an >> rn >> chain >> ri >> c[0] >> c[1] >> c[2];  // Read the whole line.
			// Add the data from the line to the Molecule's data vectors.
			mol.chain.push_back(chain); mol.ai.push_back(ai); mol.an.push_back(an); mol.ri.push_back(ri); mol.rn.push_back(rn);
			for (i = 0; i < 3; i++)
				mol.c[i].push_back(c[i]);
		}
		getline(in, line);  // Read the next line.
	}
}

// Write a PDB file.
void PDB::write(string file_name, Molecule& mol)
{
	ofstream of(file_name.c_str());  // Open the file as a filestream.
	char line[81];  // A "string" (really an array of chars) to store the text to be written.
	// Write every atom in the molecule.
	for (int i = 0; i < mol.size(); i++)
	{
		sprintf(line, "ATOM  %5d %4s%1c%3s %1s%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f\n", mol.ai[i], mol.an[i].c_str(), ' ', mol.rn[i].c_str(), mol.chain[i].c_str(), mol.ri[i], ' ', mol.c[0][i], mol.c[1][i], mol.c[2][i], 1.00, 0.00);  // Write the atomic data to the line.
		of << line;  // Write the line to the file.
	}
}

#endif
