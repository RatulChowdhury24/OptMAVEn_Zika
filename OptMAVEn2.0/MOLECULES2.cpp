/* This file is the source for the Molecule object of OptMAVEn 2.0
*/

// Include MOLECULES.cpp only once.
#ifndef MOLECULES_CPP
#define MOLECULES_CPP

#include <algorithm>
#include <array>
#include <boost/lambda/lambda.hpp>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include "MOLECULES2.h"
#include "UTILITIES2.h"


// Get the number of atoms in the Molecule.
int Molecule::size()
{
	return c.size();
}

// Clear all atoms from the Molecule.
void Molecule::clear()
{
	c.clear(); ai.clear(); an.clear(); ri.clear(); rn.clear(); chain.clear(); residue_selection_center_of_geometry.clear();
}

// Load molecular data from a file.
void Molecule::load(std::string file_name, std::string format)
{
	clear();  // Clear any previous atoms.
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
		snprintf(message, 81, "Cannot load file format: %s", format.c_str());
		throw message;
	}
	// Error if the Molecule loaded no atoms.
	if (size() == 0)
	{
		throw "Loading error: Molecule contains zero atoms.";
	}
	center_of_geometry = get_center();  // Calculate the center of geometry.
}

// Write the molecular data to a file.
void Molecule::write(std::string file_name, std::string format)
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

// Translate the molecule using a translation std::vector.
void Molecule::translate(std::array<float, 3> trans)
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
	// Also translate the centers of geometry.
	for (j = 0; j < 3; j++)
		center_of_geometry[j] += trans[j];
	for (i = 0; i < residue_selection_center_of_geometry.size(); i++)
	{
		for (j = 0; j < 3; j++)
			residue_selection_center_of_geometry[i][j] += trans[j];
	}
}

// Rotate the molecule by an angle around an axis.
void Molecule::rotate(std::array<float, 3> axis, float angle, bool degrees=true)
{
	// Trying to rotate a Molecule with no atoms suggests an error in the workflow.
	if (size() == 0)
		throw "Cannot rotate a Molecule with 0 atoms.";
	std::array<std::array<float, 3>, 3> mat = rotation_matrix3(axis, angle, degrees);  // Make the rotation matrix.
	int i, j, k;  // Declare index variables for the matrix multiplication.
	std::array<float, 3> initial;  // Declare an std::array to hold the initial coordinates.
	for (i = 0; i < size(); i++)
	{
		// Copy the initial coordinates.
		//for (j = 0; j < 3; j++)
		initial = c[i];  // initial[j] = c[i][j];
		// Multiply the old coordinates and the rotation matrix.
		for (j = 0; j < 3; j++)
		{
			// Set the new coordinate to the dot product of the row and column.
			c[i][j] = 0.0;  // Clear the old coordinate.
			for (k = 0; k < 3; k++)
				c[i][j] += initial[k] * mat[k][j];
		}
	}
	// Also rotate the centers of geometry in the same way.
	initial = center_of_geometry;
	center_of_geometry = {0.0, 0.0, 0.0};
	for (j = 0; j < 3; j++)
	{
		for (k = 0; k < 3; k++)
			center_of_geometry[k] += initial[k] * mat[k][j];
	}
	for (i = 0; i < residue_selection_center_of_geometry.size(); i++)
	{
		initial = residue_selection_center_of_geometry[i];
		residue_selection_center_of_geometry[i] = {0.0, 0.0, 0.0};
		for (j = 0; j < 3; j++)
		{
			for (k = 0; k < 3; k++)
				residue_selection_center_of_geometry[i][k] += initial[k] * mat[k][j];
		}
	}
}

// Calculate the center of the Molecule or a subset of its residues.
std::array<float, 3> Molecule::get_center()
{
	// Cannot get center if there are no atoms.
	if (size() == 0)
		throw "Cannot get center of Molecule with 0 atoms.";
	int i, j;  // Declare index variables.
	std::array<float, 3> center = {0.0, 0.0, 0.0};  // Declare an std::array to hold the coordinates.
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
	}
	return center;
}

// (Overloaded) Translate the Molecule so that its geometric center lies at the given point.
void Molecule::translate_center(std::array<float, 3> point)
{
	std::array<float, 3> center = get_center();  // Compute the center of the Molecule.
	std::array<float, 3> trans;
	// The translation std::vector is the desired center point minus the current center point.
	for (int i = 0; i < 3; i++)
		trans[i] = point[i] - center[i];
	translate(trans);
}

// (Overloaded) Center the Molecule on the origin.
void Molecule::translate_center()
{
	std::array<float, 3> origin = {0.0, 0.0, 0.0};
	translate_center(origin);
}

// Select a group of residues.
void Molecule::select_residues(std::vector<int> res)
{
	std::cout << "Select residues" << std::endl;
	// Ensure that all specified residues are in the Molecule.
	std::unordered_set<int> mol_res;  // Make an unordered set of residue numbers for faster searching.
	mol_res.insert(ri.begin(), ri.end());
	std::cout << "Verifying residues:";
	for (std::vector<int>::iterator it = res.begin(); it != res.end(); it++)
	{
		std::cout << " " << *it;
		if (mol_res.find(*it) == mol_res.end())
		{
			char message[81];
			snprintf(message, 81, "Cannot select non-existant residue %u.", *it);
			throw message;
		}
	}
	// If all residues are valid, add them to the selections.
	residue_selections.push_back(res);
	// Calculate the center of geometry of those residues.
	residue_selection_center_of_geometry.push_back(get_residue_center(res));
	std::cout << std::endl;
}

// Get a group of selected residues.
std::vector<int> Molecule::get_residue_selection(int selection)
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

// (Overloaded) Calculate the geometric center of a group of residues.
std::array<float, 3> Molecule::get_residue_center(std::vector<int> residues)
{
	std::cout << "Calculating center of residues." << std::endl;
	int i, j;  // Declare index variables.
	std::array<float, 3> center = {0.0, 0.0, 0.0};  // Declare an std::array to hold the coordinates.
	// Make an unordered set from the residues for faster searching.
	std::unordered_set<int> resiset;
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
		center[j] /= atoms_selected;
	std::cout << "Center: " << center[0] <<", "<<center[1]<<", "<<center[2]<<std::endl;
	return center;
}

// (Overloaded) Get the center of geometry of a selection of residues.
std::array<float, 3> Molecule::get_residue_center(int selection)
{
	std::cout << "Getting center of selection " << selection << std::endl;
	if (selection >= residue_selection_center_of_geometry.size())
	{
		throw "Cannot get center of invalid residue group number.";
	}
	return residue_selection_center_of_geometry[selection];
}

// Translate a residue selection to the origin.
void Molecule::translate_residues_center(int selection)
{
	std::array<float, 3> residue_center = get_residue_center(selection);  // Find the center of the residue group.
	std::array<float, 3> difference;
	// Compute the std::vector from the residue center to the Molecule center.
	for (int i = 0; i < 3; i++)
		difference[i] = get_center()[i] - residue_center[i];
	// Translating the Molecule center to the point represented by "difference" translates the residue center to the origin.
	translate_center(difference);
}

// Select a subset of residues to be the epitope.
void Antigen::define_epitope(std::vector<int> res)
{
	std::cout << "Define epitope" << std::endl;
	residue_selections.clear();  // Erase the previous epitope, if there was one. So far, only let an Antigen have one epitope at a time.
	select_residues(res);
}

// Calculate the center of geometry of the epitope.
std::array<float, 3> Antigen::get_epicenter()
{
	// Make sure there is an epitope.
	if (residue_selections.size() == 0)
		throw "Cannot find center of epitope that was not defined.";
	return get_residue_center(0);  // The epitope is at position 0.
}

// Translate the Antigen's epitope to the origin.
void Antigen::translate_epicenter()
{
	translate_residues_center(0);  // The epitope is at position 0.
}

// Rotate the Antigen so that the z coordinates of its epitope are minimized.
void Antigen::epitope_zmin()
{
	write("mount0.pdb", "PDB");
	translate_center();  // Move the Antigen to the origin so that the center of geometry does not compromise the rotation.
	write("mount1.pdb", "PDB");
	// Get the epitope's center of geometry.
	std::array<float, 3> epivec = get_epicenter();
	print_array(epivec, "Center of epitope in mount1");
	// Compute the angle between epivec and the negative z axis (equal to -k, where k is the unit std::vector <0, 0, 1>) using the formula angle = acos(a * b / (|a||b|)), where a is epivec and b is -k.
	float angle = acos(-epivec[2] / norm3(epivec));  // a * b = epivec[0] * 0 + epivec[1] * 0 + epivec[2] * -1; |a| = epinorm and |b| = 1.
	// Compute the axis of rotation, which is the cross product of -k and epivec, or -k x epivec.
	std::array<float, 3> axis = {epivec[1], -epivec[0], 0.0};  // -k x epivec = <0 * epivec[2] - -1 * epivec[1], 0 * epivec[1] - 1 * epivec[0], 0 * epivec[1] - 0 * epivec[0]>
	print_array(axis, "Axis of rotation");
	rotate(axis, angle, false);  // Rotate.
	write("mount2.pdb", "PDB");
	print_array(residue_selection_center_of_geometry[0], "Center of epitope in mount2");
	print_array(get_residue_center(residue_selections[0]), "Center of epitope in mount2");
}

// Rotate the antigen so that the epitope points towards the negative z axis and center the epitope at the origin.
void Antigen::mount_epitope()
{
	epitope_zmin();  // Minimize the z coordinates of the epitope.
	translate_epicenter();  // Move the center of the epitope to the origin.	
}

// Parse a PDB file.
void PDB::parse(std::string file_name, Molecule& mol)
{
	// Open the file in read mode.
	FILE* file;
	file = fopen(file_name.c_str(), "r");
	// Make sure the file was properly opened.
	if (file == NULL)
	{
		char message[128];
		snprintf(message, 128, "Cannot open file: %s", file_name.c_str());
		throw std::string(message);
	}
	// Declare containers for data that will be loaded.
	const int num_fields = 11;  // Get 11 data from each line of the PDB.
	int ai, ri;	
	float x, y, z, o, b;	
	char line[81];
	char type[6], chain[2], an[6], rn[4], alt1[1], alt2[1];
	std::string stype;
	// Read all of the lines.
	while (fgets(line, 81, file))
	{
		// Skip blank lines.
		if (line[0] == '\n')
			continue;
		int status = sscanf(line, "%6s %5d %6s %3s %2s %4d %8f %8f %8f %6f %6f", type, &ai, an, rn, chain, &ri, &x, &y, &z, &o, &b);  // Parse the line and extract the data.
		// Skip lines that do not start with ATOM or HETATM.
		stype = std::string(type);
		if (stype != "ATOM" && stype != "HETATM")
			continue;
		// If the line was not skipped, ensure the read was good.
		if (status != num_fields)
		{
			char message[128];
			snprintf(message, 128, "Bad line for PDB: %s", line);
			fclose(file);  // Close the file.
			throw std::string(message);			
		}
		// Add the data from the line to the Molecule's data vectors.
		mol.chain.push_back(chain); mol.ai.push_back(ai); mol.an.push_back(an); mol.ri.push_back(ri); mol.rn.push_back(rn);
		std::array<float, 3> xyz = {x, y, z}; mol.c.push_back(xyz);
	}
	fclose(file);  // Close the file.
}

// Write a PDB file.
void PDB::write(std::string file_name, Molecule& mol)
{
	std::ofstream of(file_name.c_str());  // Open the file as a filestream.
	char line[81];  // A "std::string" (really an std::array of chars) to store the text to be written.
	// Write every atom in the molecule.
	for (int i = 0; i < mol.size(); i++)
	{
		sprintf(line, "ATOM  %5d %4s%1c%3s %1s%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f\n", mol.ai[i], mol.an[i].c_str(), ' ', mol.rn[i].c_str(), mol.chain[i].c_str(), mol.ri[i], ' ', mol.c[i][0], mol.c[i][1], mol.c[i][2], 1.00, 0.00);  // Write the atomic data to the line.
		of << line;  // Write the line to the file.
	}
}

#endif
