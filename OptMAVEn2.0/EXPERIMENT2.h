/* This file sets up an Experiment, which loads, stores, and saves information about an OptMAVEn experiment.
*/

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <fstream>
#include <string>
#include <vector>
#include "STANDARDS2.h"


/* A parameter file should contain all of the following information:

AntibodyH   /path/to/structure/of/antibody_H_chain structure_file_format
AntibodyK   /path/to/structure/of/antibody_K_chain structure_file_format
Antigen     /path/to/structure/of/antigen structure_file_format
ZAngles     angle1 angle2 ... angleN

*/


class Experiment
{
public:
	// Attributes to be loaded.
	std::string antigen_file_path;
	std::string antigen_file_format;
	std::string chain_H_file_path;
	std::string chain_H_file_format;
	std::string chain_K_file_path;
	std::string chain_K_file_format;
	std::vector<int> epitope;
	std::vector<float> z_rots;
	std::vector<float> x_trans;
	std::vector<float> y_trans;
	std::vector<float> z_trans;
	int max_parallel;
	// Load a single line of floats into a std::vector.
	std::vector<float> load_line_float(std::stringstream& line, std::string what)
	{
		std::vector<float> values;
		float value;
		while (!line.eof())
		{
			line >> value;
			if (line.fail())
			{
				char message[128];
				snprintf(message, 128, "Bad value for %s.", what.c_str());
				throw std::string(message);
			}
			values.push_back(value);
		}
		return values;
	}
	// Load a single line of ints into a std::vector.
	std::vector<int> load_line_int(std::stringstream& line, std::string what)
	{
		std::cout << "Loading epitope:";
		std::vector<int> values;
		int value;
		while (!line.eof())
		{
			line >> value;
			if (line.fail())
			{
				char message[128];
				snprintf(message, 128, "Bad value for %s.", what.c_str());
				throw std::string(message);
			}
			std::cout << " " << value;
			values.push_back(value);
		}
		std::cout << std::endl;
		return values;
	}
	// Load an Experiment from a parameter file.
	void load(std::string parameter_file_path)
	{
		std::ifstream in(parameter_file_path.c_str());  // Open the file.
		// Abort if the parameter file cannot be opened.
		if (!in.is_open())
		{
			char message[128];
			snprintf(message, 128, "Cannot open file: %s", parameter_file_path.c_str());
			throw std::string(message);
		}
		// Make variables to track the loading.
		std::string parameter;
		std::string line;
		int int_in;
		float float_in;
		std::string bad_in;
		// Read until the end of the file.
		getline(in, line);
		while (!in.eof())
		{
			std::stringstream info(line);  // Stream in the line.
			info >> parameter;  // Read the parameter name.
			// Determine where to put the value.
			if (parameter == parameter_Ag())
				info >> antigen_file_path >> antigen_file_format;
			else if (parameter == parameter_IgH())
				info >> chain_H_file_path >> chain_H_file_format;
			else if (parameter == parameter_IgK())
				info >> chain_K_file_path >> chain_K_file_format;
			else if (parameter == parameter_epitope())
			{
				epitope = load_line_int(info, std::string("epitope residue index"));
			}
			else if (parameter == parameter_z_rots())
			{
				z_rots = load_line_float(info, std::string("z rotation angle"));
			}
			else
			{
				// A bad parameter indicates a malformatted file.
				char message[128];
				snprintf(message, 128, "Bad experimental parameter: %s", parameter.c_str());
				throw std::string(message);
			}
			getline(in, line);
		}
	}
	// Calculate the total number of combinations of translations and rotation angles, not considering that some combinations may involve clashes.
	int total_combinations()
	{
		return z_rots.size() * x_trans.size() * y_trans.size() * z_trans.size();
	}
	// Make a file that says which combinations of rotations and translations will be assigned to each cull_clashes process.
	void write_cull_clashes_jobs()
	{
		int pid = 0;
		
	}
	Experiment() {};
};

#endif
