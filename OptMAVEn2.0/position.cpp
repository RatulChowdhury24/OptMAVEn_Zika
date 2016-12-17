/* This file contains the source code for creating the antigen positions. It should be run first in OptMAVEn. 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "MOLECULES2.cpp"
#include "UTILITIES.h"
using namespace std;


/* This program begins by reading the parameter file, which should contain all of the following information:

Angle       angle_begin  angle_end    angle_step
Antigen     /path/to/structure/of/antigen structure_file_format
AntibodyH   /path/to/structure/of/antibody_H_chain structure_file_format
AntibodyK   /path/to/structure/of/antibody_K_chain structure_file_format
*/

/* USAGE of argv:
1: path to file containing parameters
*/
int main(int argc, char* argv[])
{
	string parameter_file_path = argv[1];  // Get the location of the parameter file.
	ifstream in(parameter_file_path.c_str());  // Open the parameter file.
	string line;
	string parameter;
	float angle_begin;
	float angle_end;
	float angle_step;
	string antigen_file_path;
	string antigen_file_format;
	string chain_H_file_path;
	string chain_H_file_format;
	string chain_K_file_path;
	string chain_K_file_format;
	int residue;
	vector<int> epitope;
	getline(in, line);
	// Read until the end of the file.
	while (!in.eof())
	{
		stringstream info(line);  // Stream in the line.
		info >> parameter;  // Read the parameter name.
		// Determine where to put the value.
		if (parameter == "Angle")
			info >> angle_begin >> angle_end >> angle_step;
		else if (parameter == "Antigen")
			info >> antigen_file_path >> antigen_file_format;
		else if (parameter == "AntibodyH")
			info >> chain_H_file_path >> chain_H_file_format;
		else if (parameter == "AntibodyK")
			info >> chain_K_file_path >> chain_K_file_format;
		else if (parameter == "Epitope")
		{
			while (!info.eof())
			{
				info >> residue;
				epitope.push_back(residue);
			}
		}		
		getline(in, line);
	}
	// Make a molecule of the antigen.
	Antigen antigen;
	cout << "Loading antigen" << endl;
	antigen.load(antigen_file_path, antigen_file_format);
	cout << "Setting epitope" << endl;
	antigen.set_epitope(epitope);  // Tell the antigen which residues are part of the epitope.	
	cout << "Rotating zmin" << endl;
	antigen.epitope_zmin();  // Rotate the antigen so as to minimize the z coordinates of its epitope.
	cout << "Writing PDB" << endl;
	antigen.write("zneg.pdb", "PDB");  //FIXME
	return 0;
}
