/* This file is the source for the Molecule object of OptMAVEn 2.0
*/

// Include MOLECULES.cpp only once.
#ifndef MOLECULES_CPP
#define MOLECULES_CPP

#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include "MOLECULES2.h"
#include "UTILITIES2.h"
using namespace std;


// Get the number of atoms in the Molecule.
int Molecule::size()
{
	return c.size();
}

// Clear all atoms from the Molecule.
void Molecule::clear()
{
	c.clear(); ai.clear(); an.clear(); ri.clear(); rn.clear(); chain.clear();
}

// Load molecular data from a file.
void Molecule::load(string file_name, string format)
{
	clear();  // Clear any previous atoms.
	center_is_known = false;
	if (format == "PDB")
	{
		// Load from a PDB file.
		PDB parser;
		parser.parse(file_name, *this);
	}
	else
	{
		// Cannot load from an unrecognized file format.
		char message[81];
		snprintf(message, 81, "Load failed (invalid file format: %s)", format.c_str());
		throw message;
	}
}

// Write the molecular data to a file.
void Molecule::write(string file_name, string format)
{
	if (format == "PDB")
	{
		// Write a PDB file.
		PDB writer;
		writer.write(file_name, *this);
	}
	else
	{
		// Cannot write to an unrecognized file format.
		char message[81];
		snprintf(message, 81, "Write failed (invalid file format: %s)", format.c_str());
		throw message;
	}
}

// Translate the molecule using a translation vector.
void Molecule::translate(array<float, 3> trans)
{
	// Trying to translate a Molecule with no atoms suggests an error in the workflow.
	if (size() == 0)
		throw "Cannot translate a Molecule with 0 atoms.";
	int i, j;  // Declate index variables.
	// Translate all points.
	for (i = 0; i < size(); i++)
	{
		for (j = 0; j < 3; j++)
			c[i][j] += trans[j];
	}
}

// Rotate the molecule by an angle around an axis.
void Molecule::rotate(array<float, 3> axis, float angle, bool degrees=true)
{
	// Trying to rotate a Molecule with no atoms suggests an error in the workflow.
	if (size() == 0)
		throw "Cannot rotate a Molecule with 0 atoms.";
	array<array<float, 3>, 3> mat = rotation_matrix3(axis, angle, degrees);  // Make the rotation matrix.
	int i, j, k;  // Declare index variables for the matrix multiplication.
	array<float, 3> initial;  // Declare an array to hold the initial coordinates.
	for (i = 0; i < size(); i++)
	{
		// Copy the initial coordinates.
		for (j = 0; j < 3; j++)
			initial[j] = c[i][j];
		// Multiply the old coordinates and the rotation matrix.
		for (j = 0; j < 3; j++)
		{
			// Set the new coordinate to the dot product of the row and column.
			c[i][j] = 0.0;  // Clear the old coordinate.
			for (k = 0; k < 3; k++)
				c[i][j] += initial[k] * mat[k][j];
		}
	}
}

// (Overloaded) Calculate the center of the Molecule or a subset of its residues.
array<float, 3> Molecule::get_center()
{
	// To save time, return the central coordinates if they are known.
	if (center_is_known)
		return known_center;
	// Cannot get center if there are no atoms.
	if (size() == 0)
		throw "Cannot get center of Molecule with 0 atoms.";
	int i, j;  // Declare index variables.
	array<float, 3> center = {0.0, 0.0, 0.0};  // Declare an array to hold the coordinates.
	// Loop through all atoms.
	for (i = 0; i < size(); i++)
	{
		// Add the coordinates of that atom.
		for (j = 0; j < 3; j++)
			center[j] += c[i][j];
	}
	// Average the coordinates.
	for (j = 0; j < 3; j++)
	{
		center[j] /= size();
		known_center[j] = center[j];  // Memorize the coordinates.
	}
	center_is_known = true;
	return center;
}

// Translate the Molecule so that its geometric center lies at the given point.
void Molecule::translate_center(array<float, 3> point)
{
	array<float, 3> center = get_center();  // Compute the center of the Molecule.
	array<float, 3> trans;
	// The translation vector is the desired center point minus the current center point.
	for (int i = 0; i < 3; i++)
		trans[i] = point[i] - center[i];
	translate(trans);
}

// Center the Molecule on the origin.
void Molecule::translate_center()
{
	array<float, 3> origin = {0.0, 0.0, 0.0};
	translate_center(origin);
}

// Select a group of residues.
void Molecule::select_residues(vector<int> res)
{
	// Ensure that all specified residues are in the Molecule.
	unordered_set<int> mol_res;
	mol_res.insert(ri.begin(), ri.end());
	for (vector<int>::iterator it = res.begin(); it != res.end(); it++)
	{
		if (mol_res.find(*it) == mol_res.end())
		{
			char message[81];
			snprintf(message, 81, "Cannot select non-existant residue %u.", *it);
			throw message;
		}
	}
	// If all residues are valid, add them to the selections.
	residue_selections.push_back(res);
}

// Get a group of selected residues.
vector<int> Molecule::get_residue_selection(int selection)
{
	// Ensure the selection is valid.
	if (selection >= residue_selections.size())
	{
		char message[81];
		snprintf(message, 81, "Cannot select non-existant residue group %u.", selection);
		throw message;
	}
	return residue_selections[selection];
}

// Calculate the geometric center of a group of selected residues.
array<float, 3> Molecule::get_residue_center(int selection)
{
	vector<int> residues = get_residue_selection(selection);  // Get the selected residues.
	int i, j;  // Declare index variables.
	array<float, 3> center = {0.0, 0.0, 0.0};  // Declare an array to hold the coordinates.
	// Make a set from the residues for faster searching.
	unordered_set<int> resiset;
	resiset.insert(residues.begin(), residues.end());
	int atoms_selected = 0;
	// Examine each atom in the Molecule.
	for (i = 0; i < size(); i++)
	{
		// Test for membership in the selection.
		if (resiset.find(ri[i]) != resiset.end())
		{
			atoms_selected++;  // Count the atom.
			// Add the coordinates of the atom.
			for (j = 0; j < 3; j++)
				center[j] += c[i][j];
		}
	}
	// Cannot get center if no atoms are in the selection.
	if (atoms_selected == 0)
		throw "Cannot get center of a selection that contains 0 atoms.";
	// Average the coordinates.
	for (j = 0; j < 3; j++)
		center[j] / atoms_selected;
	return center;
}

// Translate a residue selection to the origin.
void Molecule::translate_residues_center(int selection)
{
	array<float, 3> residue_center = get_residue_center(selection);  // Find the center of the residue group.
	array<float, 3> difference;
	// Compute the vector from the residue center to the Molecule center.
	for (int i = 0; i < 3; i++)
		difference[i] = get_center()[i] - residue_center[i];
	// Translating the Molecule center to the point represented by "difference" translates the residue center to the origin.
	translate_center(difference);	
}

// Select a subset of residues to be the epitope.
void Antigen::define_epitope(vector<int> res)
{
	residue_selections.clear();  // Erase the previous epitope, if there was one. So far, only let an Antigen have one epitope at a time.
	select_residues(res);
}

// Calculate the center of geometry of the epitope.
array<float, 3> Antigen::get_epicenter()
{
	// Make sure there is an epitope.
	if (residue_selections.size() == 0)
		throw "Cannot find center of epitope that was not defined.";
	return get_residue_center(0);  // The epitope is at position 0.
}

// Translate the Antigen's epitope to the origin.
void Antigen::translate_epicenter()
{
	translate_residues_center(0);
}

// Rotate the Antigen so that the z coordinates of its epitope are minimized.
void Antigen::epitope_zmin()
{
	translate_center();  // Move the Antigen to the origin so that the center of geometry does not compromise the rotation.
	// Get the epitope's center of geometry.
	array<float, 3> epivec = get_epicenter();
	// Compute the angle between epivec and the negative z axis (equal to -k, where k is the unit vector <0, 0, 1>) using the formula angle = acos(a * b / (|a||b|)), where a is epivec and b is -k.
	float angle = acos(-epivec[2] / norm3(epivec));  // a * b = epivec[0] * 0 + epivec[1] * 0 + epivec[2] * -1; |a| = epinorm and |b| = 1.
	// Compute the axis of rotation, which is the cross product of epivec and -k, or epivec x -k.
	array<float, 3> axis = {-epivec[1], epivec[0], 0.0};  // epivec x -k = <epivec[1] * -1 - epivec[2] * 0, epivec[2] * 0 - epivec[0] * -1, epivec[0] * 0 - epivec[1] * 0>
	rotate(axis, angle, false);  // Rotate.
}

// Parse a PDB file.
void PDB::parse(string file_name, Molecule& mol)
{
	ifstream in(file_name.c_str());  // Open the file as a filestream.
	// Declare containers for data that will be loaded.
	float x, y, z;
	string line;
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
			lstream >> ai >> an >> rn >> chain >> ri >> x >> y >> z;  // Read the whole line.
			// Add the data from the line to the Molecule's data vectors.
			mol.chain.push_back(chain); mol.ai.push_back(ai); mol.an.push_back(an); mol.ri.push_back(ri); mol.rn.push_back(rn);
			array<float, 3> xyz = {x, y, z}; mol.c.push_back(xyz);
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
		sprintf(line, "ATOM  %5d %4s%1c%3s %1s%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f\n", mol.ai[i], mol.an[i].c_str(), ' ', mol.rn[i].c_str(), mol.chain[i].c_str(), mol.ri[i], ' ', mol.c[i][0], mol.c[i][1], mol.c[i][2], 1.00, 0.00);  // Write the atomic data to the line.
		of << line;  // Write the line to the file.
	}
}

#endif
