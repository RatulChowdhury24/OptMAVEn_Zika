/* This file sets up an Experiment, which loads, stores, and saves information about an OptMAVEn experiment.
*/

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <fstream>
#include <string>
#include <vector>
using namespace std;


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
	string antigen_file_path;
	string antigen_file_format;
	string chain_H_file_path;
	string chain_H_file_format;
	string chain_K_file_path;
	string chain_K_file_format;
	vector<int> epitope;
	vector<float> z_rots;
	vector<float> x_trans;
	vector<float> y_trans;
	vector<float> z_trans;
	int max_parallel;
	// Load a single line of floats into a vector.
	vector<float> load_line_float(stringstream line, string what)
	{
		vector<float> values;
		float value;
		while (!line.eof())
		{
			line >> value;
			if (line.fail())
			{
				char message[128];
				line >> bad_in;
				snprintf(message, 128, "Bad value for %s: %s", what, bad_in.c_str());
				throw string(message);
			}
			values.push_back(value);
		}
	}
	// Load a single line of ints into a vector.
	vector<int> load_line_int(stringstream line, string what)
	{
		vector<int> values;
		int value;
		while (!line.eof())
		{
			line >> value;
			if (line.fail())
			{
				char message[128];
				line >> bad_in;
				snprintf(message, 128, "Bad value for %s: %s", what, bad_in.c_str());
				throw string(message);
			}
			values.push_back(value);
		}
	}
	// Load an Experiment from a parameter file.
	void load(string parameter_file_path)
	{
		ifstream in(parameter_file_path.c_str());  // Open the file.
		// Abort if the parameter file cannot be opened.
		if (!in.is_open())
		{
			char message[128];
			snprintf(message, 128, "Cannot open file: %s", parameter_file_path.c_str());
			throw string(message);
		}
		// Make variables to track the loading.
		string parameter;
		string line;
		int int_in;
		float float_in;
		string bad_in;
		// Read until the end of the file.
		getline(in, line);
		while (!in.eof())
		{
			stringstream info(line);  // Stream in the line.
			info >> parameter;  // Read the parameter name.
			// Determine where to put the value.
			if (parameter == parameter_Ag)
				info >> antigen_file_path >> antigen_file_format;
			else if (parameter == parameter_IgH)
				info >> chain_H_file_path >> chain_H_file_format;
			else if (parameter == parameter_IgK)
				info >> chain_K_file_path >> chain_K_file_format;
			else if (parameter == parameter_epitope)
			{
				epitope = load_line_int(info, "epitope residue index");
			}
			else if (parameter == parameter_z_rots)
			{
				z_rots = load_line(info, "z rotation angle");
			}
			else
			{
				// A bad parameter indicates a malformatted file.
				char message[128];
				snprintf(message, 128, "Bad experimental parameter: %s", parameter.c_str());
				throw string(message);
			}
			getline(in, line);
		}
	}
	// Calculate the total number of combinations of translations and rotation angles, not considering that some combinations may involve clashes.
	int total_combinations()
		return z_rots.size() * x_trans.size() * y_trans.size() * z_trans.size();
	// Make a file that says which combinations of rotations and translations will be assigned to each cull_clashes process.
	void write_cull_clashes_jobs()
	{
		int pid = 0;
		
	}
	Experiment() {};
};

#endif
