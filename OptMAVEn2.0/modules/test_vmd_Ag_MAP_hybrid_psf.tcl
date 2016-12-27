# Merge an antigen and a MAPs part into a single molecule and generate a PDB and PSF.
# The following command-line arguments are required:
# -e: the name of this script
# -f: the coordinates of the antigen
# -args: 0: the experiment's directory
#        1: the topology file

# The psfgen package is required.
package require psfgen

# Read the topology file and alias the atoms.
topology top_all27_prot_na.rtf
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# Build "segment A" for the antigen.
set AgPDB "AgInit.pdb"
segment A {
	pdb $AgPDB
}

# Read antigen coordinates into segment A.
coordpdb $AgPDB A

# Build "segment I" for the immunoglobulin part.
set IgPDB "../databases/MAPs/HCDR3/HCDR3_1.pdb"
segment I {
	pdb $IgPDB
}

# Read the immunoglobulin coordinates into segment I.
coordpdb $IgPDB I

# Create the PDB and PSF.
writepdb "AgIg.pdb"
writepsf "AgIg.psf"

exit
