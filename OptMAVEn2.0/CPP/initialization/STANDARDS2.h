/* This file defines standards that are used internally throughout OptMAVEn 2.0.
*/

#ifndef STANDARDS_H
#define STANDARDS_H

#include <string>


// The molecule file format (MFF) used withing OptMAVEn.
string MFF()
	return "PDB";

// The name of the file of the antigen structure once it has been rotated to minimize the z coordinates of its epitope and translated so that its epitope center lies at the origin.
string AgInit()
	return "AgInit.pdb";

// The name of the file that tells each cull_clashes process which orientations to test.
string cull_clashes_jobs()
	return "cull_clashes_jobs.dat";

// Names of parameters found in the parameters file.
string parameter_Ag()
	return "Antigen:";

string parameter_IgH()
	return "AntibodyH:";

string parameter_IgK()
	return "AntibodyK:";

string parameter_epitope()
	return "Epitope:";

string parameter_z_rots()
	return "ZAngles:";

string parameter_x_trans()
	return "XTrans:";

string parameter_y_trans()
	return "YTrans:";

string parameter_z_trans()
	return "ZTrans:";

#endif
