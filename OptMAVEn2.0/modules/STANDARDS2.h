/* This file defines standards that are used internally throughout OptMAVEn 2.0.
*/

#ifndef STANDARDS_H
#define STANDARDS_H

#include <string>


// The molecule file format (MFF) used withing OptMAVEn.
std::string MFF()
{
	return "PDB";
}

// The name of the file of the antigen structure once it has been rotated to minimize the z coordinates of its epitope and translated so that its epitope center lies at the origin.
std::string AgInit()
{
	return "AgInit.pdb";
}

// The name of the file that tells each cull_clashes process which orientations to test.
std::string cull_clashes_jobs()
{
	return "cull_clashes_jobs.dat";
}

// Names of parameters found in the parameters file.
std::string parameter_Ag()
{
	return "Antigen:";
}

std::string parameter_IgH()
{
	return "AntibodyH:";
}

std::string parameter_IgK()
{
	return "AntibodyK:";
}

std::string parameter_epitope()
{
	return "Epitope:";
}

std::string parameter_z_rots()
{
	return "ZAngles:";
}

std::string parameter_x_trans()
{
	return "XTrans:";
}

std::string parameter_y_trans()
{
	return "YTrans:";
}

std::string parameter_z_trans()
{
	return "ZTrans:";
}

#endif
