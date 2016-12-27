""" This script takes the raw structure of the antigen, adds missing atoms, generates a PSF, performs an energy minimization in NAMD, and orients the minimized antigen with its epitope centered at the origin and pointing along the negative z axis.
"""

import os

# Add missing atoms and generate the PSF.
os.system("vmd -e test_vmd_prepare_minimization.tcl -dispdev none")

# Generate the NAMD configuration file.
config_file = "namd_minimization.conf"
minimized_prefix = "AgMin"
minimized_coords = "{}.coor".format(minimized_prefix)
minimized_pdb = "{}.pdb".format(minimized_prefix)
lines = open("namd_minimization_base.conf").readlines()
# Define the parameters to write into the configuration file.
parameters = {
	"structure": "Ag.psf",
	"coordinates": "AgAddedAtoms.pdb",
	"parameters": "par_all27_prot_na.prm",
	"set outputname": minimized_prefix
}

# Edit the base file line by line.
with open(config_file, "w") as cf:
	for line in lines:
		if line.strip() != "":
			# If the line begins with a parameter in the dict of parameters, add the value of that parameter to the configuration file.
			parameter = line[0: min(20, len(line))].strip()
			if parameter in parameters:
				line = "{:<20s}{}\n".format(parameter, parameters[parameter])
		cf.write(line)

# Run the minimization in NAMD.
os.system("namd {}".format(config_file))

# Rename the output coordinates to clearly be a PDB.
os.rename(minimized_coords, minimized_pdb)

# Position the minimized antigen.
os.system("vmd -e test_vmd_initial_position.tcl -f {} -dispdev none".format(minimized_pdb))
