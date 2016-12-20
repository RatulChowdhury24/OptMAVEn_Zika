/* This file is the header for the Molecule object of OptMAVEn 2.0
*/

#ifndef MOLECULES_H
#define MOLECULES_H

#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;


// The Molecule class stores atom information as vectors.
class Molecule
{
public:
	// Declare containers for the following information for the Molecule:
	vector<array<float, 3>> c;  // atomic coordinates
	vector<float> r, q;  // numeric properties (r: radius, q: partial charge)
	vector<int> ai, ri;  // atomic (ai) and residue (ri) indices
	vector<string> chain, an, rn;  // chain, atomic (an) and residue (rn) names
	// Get the number of atoms in the Molecule.
	int size();
	// Clear all atoms.
	void clear();
	// Load molecular data from a file.
	void load(string file_name, string format);
	// Print the molecular data.
	void write(string file_name, string format);
	// Translate the molecule using a translation vector.
	void translate(array<float, 3>);
	// Rotate the molecule using an angle and axis of rotation.
	void rotate(array<float, 3>, float, bool);
	// Calculate the center of geometry of the whole molecule.
	array<float, 3> get_center();
	// Calculate the center of geometry of a group of residues.
	array<float, 3> get_residue_center(int);
	// Translate the Molecule to the origin.
	void translate_center();
	// Translate the Molecule to a specific point.
	void translate_center(array<float, 3>);
	// Is the center known?
	bool center_is_known;
	array<float, 3> known_center;
	// Select residues in the Molecule.
	vector<vector<int>> residue_selections;
	void select_residues(vector<int>);
	// Get a group of selected residues.
	vector<int> get_residue_selection(int);
	// Center a residue selection at the origin.
	void translate_residues_center(int);
	// Default constructor.
	Molecule() : center_is_known(false) {residue_selections.clear();};
};


// A type of Molecule against which to design an antibody. It has special members relating to its epitope.
class Antigen : public Molecule
{
public:
	// Select a subset of residues to be the epitope.
	void define_epitope(vector<int>);
	// Get the number of atoms in the epitope.
	int epitope_size();
	// Get the geometric center of the epitope.
	array<float, 3> get_epicenter();
	// Translate the epitope center to the origin.
	void translate_epicenter();
	// Rotate the epitope so that the z coordinates of its residues are minimized.
	void epitope_zmin();
	// Is the epicenter known?
	//bool epicenter_is_known;
	//array<float, 3> known_epicenter;
	// Default constructor.
	Antigen() {};
};


class FileFormat
{
public:
	// Load a Molecule from a file.
	virtual void parse(string, Molecule&) = 0;
	// Write a Molecule to a file.
	virtual void write(string, Molecule&) = 0;
};


class PDB : public FileFormat
{
public:
	// Parse a PDB file.
	void parse(string, Molecule&);
	//bool parse(string, vector<string>&, vector<int>&, vector<string>&, vector<int>&, vector<string>&, vector<float>);
	// Write a PDB file.
	void write(string, Molecule&);
};

#endif
