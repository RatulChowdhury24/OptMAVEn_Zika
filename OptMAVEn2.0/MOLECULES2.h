/* This file is the header for the Molecule object of OptMAVEn 2.0
*/

#ifndef MOLECULES_H
#define MOLECULES_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


// The Molecule class stores atom information as vectors.
class Molecule
{
private:
	// Ensure that the data vectors have the same length.
	bool checkXYZ();
public:
	// Declare containers for the following information for the Molecule:
	vector<float> c[3];  // atomic coordinates
	vector<float> r, q;  // numeric properties (r: radius, q: partial charge)
	vector<int> ai, ri;  // atomic (ai) and residue (ri) indices
	vector<string> chain, an, rn;  // chain, atomic (an) and residue (rn) names
	// Get the number of atoms in the Molecule.
	int size();
	// Load molecular data from a file.
	bool load(string file_name, string format);
	// Print the molecular data.
	bool write(string file_name, string format);
	// Translate the molecule using a translation vector.
	void translate(float[3]);
	// Rotate the molecule using a rotation matrix.
	void rotate(float[3][3]);
	// Calculate the center of geometry of the whole molecule.
	void get_center(float[3]);
	// Calculate the center of geometry of a subset of residues.
	void get_center(float[3], vector<int>);
	// Translate the Molecule to the origin.
	void translate_center();
	// Translate the Molecule to a specific point.
	void translate_center(float[3]);
	// Default constructor.
	Molecule() {};
};


// A type of Molecule against which to design an antibody. It has special members relating to its epitope.
class Antigen : public Molecule
{
public:
	vector<int> epitope;  // Identify a subset of residues as the epitope.
	// Select a subset of residues to be the epitope.
	void set_epitope(vector<int>);
	// Get the number of atoms in the epitope.
	int epitope_size();
	// Get the geometric center of the epitope.
	void epicenter(float[3]);
	// Rotate the epitope so that the z coordinates of its residues are minimized.
	void epitope_zmin();
};


class FileFormat
{
public:
	// Load a Molecule from a file.
	virtual bool parse(string, Molecule&) = 0;
	// Write a Molecule to a file.
	virtual void write(string, Molecule&) = 0;
};


class PDB : public FileFormat
{
public:
	// Parse a PDB file.
	bool parse(string, Molecule&);
	//bool parse(string, vector<string>&, vector<int>&, vector<string>&, vector<int>&, vector<string>&, vector<float>);
	// Write a PDB file.
	void write(string, Molecule&);
};

#endif
