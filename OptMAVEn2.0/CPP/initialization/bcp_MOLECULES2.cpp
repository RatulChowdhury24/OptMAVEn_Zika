/* This file is the header for the Molecule object of OptMAVEn 2.0
*/

// Include MOLECULES.h only once.
#ifndef MOLECULES
#define MOLECULES

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


class Molecule;


class FileFormat
{
public:
	// Load a Molecule from a file.
	virtual bool parse(string) {return false;}
	// Write a Molecule to a file.
	//virtual void write(Molecule) {}
};


class PDB : public FileFormat
{
public:
	// Parse a PDB file.
	bool parse(string file_name, vector<string>& vch, vector<int>& vai, vector<string>& van, vector<int>& vri, vector<string>& vrn, vector<float> vc[3])
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
				vch.push_back(chain); vai.push_back(ai); van.push_back(an); vri.push_back(ri); vrn.push_back(rn);
				for (i = 0; i < 3; i++)
					vc[i].push_back(c[i]);
			}
			getline(in, line);  // Read the next line.
		}
	}
	// Write a PDB file.
	void write(Molecule&, string);
};


// The Molecule class stores atom information as vectors.
class Molecule
{
private:
	// Ensure that the data vectors have the same length.
	bool checkXYZ()
	{
		return ((c[0].size() == c[1].size()) && (c[1].size() == c[2].size()));
	}
public:
	// Declare containers for the following information for the Molecule:
	vector<float> c[3];  // atomic coordinates
	vector<float> r, q;  // numeric properties (r: radius, q: partial charge)
	vector<int> ai, ri;  // atomic (ai) and residue (ri) indices
	vector<string> chain, an, rn;  // chain, atomic (an) and residue (rn) names
	// Get the number of atoms in the Molecule.
	int size()
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
	bool load(string file_name, string format)
	{
		if (format == "PDB")
		{
			// Load from a PDB file.
			PDB parser;
			return parser.parse(file_name, chain, ai, an, ri, rn, c);
		}
		else
		{
			// Cannot load from an unrecognized file format.
			return false;
		}
	}
	// Rotate the molecule using a rotation matrix.
	void rotate(float mat[3][3])
	{
		int i, j, k;  // Declare index variables for the matrix multiplication.
		// Declare an array to hold the old coordinates of the point.
		float old[3];
		// Change all points in the vectors.
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
			/*
			// Dot the coordinates of this point with the columns of the rotation matrix. This operation is hard-coded, not looped, to speed it up.
			x_old = x.at(i); y_old = y.at(i);
			x.at(i) = x_old * mat[0][0] + y_old * mat[1][0] + z.at(i) * mat[2][0];
			y.at(i) = x_old * mat[0][1] + y_old * mat[1][1] + z.at(i) * mat[2][1];
			z.at(i) = x_old * mat[0][2] + y_old * mat[1][2] + z.at(i) * mat[2][2];
			*/
		}
	}
	// Print the molecular data.
	bool write(string file_name, string format)
	{
		// Ensure the coordinates are consistent.
		if (!checkXYZ())
			return false;
		if (format == "PDB")
		{
			// Write a PDB file.
			PDB writer;
			writer.write(*this, file_name);
		}
		else
		{
			// Cannot write to an unrecognized file format.
			return false;
		}
		return true;
	}
	// Default constructor.
	Molecule() {};
};


// Write a PDB file.
void PDB::write(Molecule& mol, string file_name)
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
