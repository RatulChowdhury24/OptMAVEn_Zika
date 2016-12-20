/* This file contains the source code for creating the antigen positions. It should be run first in OptMAVEn. 
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "EXPERIMENT2.h"
#include "MOLECULES2.cpp"
#include "STANDARDS2.h"
#include "UTILITIES2.h"


/* USAGE of argv:
1: path to file containing parameters	
*/
int main(int argc, char* argv[])
{
	std::string parameter_file_path = argv[1];  // Get the location of the parameter file.	
	try
	{
		// Load the Experiment from the parameters.
		Experiment ex;
		ex.load(parameter_file_path);
		// Make a Molecule for the antigen.
		Antigen antigen;
		std::cout << "Loading Antigen" << std::endl;
		antigen.load(ex.antigen_file_path, ex.antigen_file_format);
		antigen.write("scanf.pdb", "PDB");
		std::cout << "Loading Epitope" << std::endl;
		antigen.define_epitope(ex.epitope);  // Tell the antigen which residues are part of the epitope.
		std::cout << "Mounting" << std::endl;
		antigen.mount_epitope();  // Move the antigen so its epitope is centered at the origin and points in the negative z direction.
		std::cout << "Writing" << std::endl;
		antigen.write(AgInit(), MFF());  // Save a file of the antigen in its new position.
		for (int i = 0; i < 3; i++)
		{
			std::cout << antigen.center_of_geometry[i] << ", ";
		}
		for (int i = 0; i < 3; i++)
		{
			std::cout << antigen.get_epicenter()[i] << ", ";
		}
	}
	catch (std::string msg)
	{
		std::cout << msg << std::endl;
		exit(EXIT_FAILURE);
	}	
		
	// FIXME: let the experiment load these values.
	/*
	// Make molecules for the antibodies.
	Molecule IgH; IgH.load("MoleculeH.pdb", "PDB");
	// Find positions in which the antibodies and antigen do not clash.
	std::cout << "Valid" << std::endl;
	for (float zt = 0.0; zt < 360.0; zt += 20.0)
	{
		// Simulate what it would be like to need to load the structure for every time.
		antigen.load(AgInit(), MFF());
		antigen.define_epitope(ex.epitope);  // Tell the antigen which residues are part of the epitope.
		array<float, 3> negz = {0, 0, -1};
		antigen.rotate(negz, zt);
		char name[81];
		snprintf(name, 81, "Agt%f", zt);
		antigen.write(std::string(name), "PDB");
		if (!clashing(antigen, IgH))
			std::cout << zt << std::endl;
	}
	*/
	return 0;
}
