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


// The Molecule class stores atom information as vectors.
class Molecule
{
public:
	// Declare containers for the following information for the Molecule:
	std::vector<std::array<float, 3>> c;  // atomic coordinates
	std::vector<float> r, q;  // numeric properties (r: radius, q: partial charge)
	std::vector<int> ai, ri;  // atomic (ai) and residue (ri) indices
	std::vector<std::string> chain, an, rn;  // chain, atomic (an) and residue (rn) names
	// Get the number of atoms in the Molecule.
	int size();
	// Clear all atoms.
	void clear();
	// Load molecular data from a file.
	void load(std::string file_name, std::string format);
	// Print the molecular data.
	void write(std::string file_name, std::string format);
	// Translate the molecule using a translation std::vector.
	void translate(std::array<float, 3>);
	// Rotate the molecule using an angle and axis of rotation.
	void rotate(std::array<float, 3>, float, bool);
	// Calculate the center of geometry of the whole molecule.
	std::array<float, 3> get_center();
	// Calculate the center of geometry of a group of residues.
	std::array<float, 3> get_residue_center(std::vector<int> residues);
	std::array<float, 3> get_residue_center(int);  // And for a selection of residues.
	// Translate the Molecule to the origin.
	void translate_center();
	// Translate the Molecule to a specific point.
	void translate_center(std::array<float, 3>);
	// The center of geometry.
	std::array<float, 3> center_of_geometry;
	// Select residues in the Molecule.
	std::vector<std::vector<int>> residue_selections;
	void select_residues(std::vector<int>);
	// The center of geometry of each selection.
	std::vector<std::array<float, 3>> residue_selection_center_of_geometry;
	// Get a group of selected residues.
	std::vector<int> get_residue_selection(int);
	// Center a residue selection at the origin.
	void translate_residues_center(int);
	// Default constructor.
	Molecule() {};
};


// A type of Molecule against which to design an antibody. It has special members relating to its epitope.
class Antigen : public Molecule
{
public:
	// Select a subset of residues to be the epitope.
	void define_epitope(std::vector<int>);
	// Get the number of atoms in the epitope.
	int epitope_size();
	// Get the geometric center of the epitope.
	std::array<float, 3> get_epicenter();
	// Translate the epitope center to the origin.
	void translate_epicenter();
	// Rotate the epitope so that the z coordinates of its residues are minimized.
	void epitope_zmin();
	// Rotate the antigen so that the epitope points towards the negative z axis and center the epitope at the origin.
	void mount_epitope();
	// Default constructor.
	Antigen() {};
};


class FileFormat
{
public:
	// Load a Molecule from a file.
	virtual void parse(std::string, Molecule&) = 0;
	// Write a Molecule to a file.
	virtual void write(std::string, Molecule&) = 0;
};


class PDB : public FileFormat
{
public:
	// Parse a PDB file.
	void parse(std::string, Molecule&);
	//bool parse(std::string, std::vector<std::string>&, std::vector<int>&, std::vector<std::string>&, std::vector<int>&, std::vector<std::string>&, std::vector<float>);
	// Write a PDB file.
	void write(std::string, Molecule&);
};

#endif
