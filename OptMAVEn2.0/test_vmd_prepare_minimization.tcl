# This VMD script adds missing atoms to the initial structure of an antigen and generates the PSF that is necessary to perform energy calculations.

# Load required external modules.
package require psfgen

# Generate a PSF of the antigen for future NAMD energy calculations.
# Load the topology file.
topology top_all27_prot_na.rtf
# Alias atoms and residues to make the structure and the topology file compatible.
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
# Generate a segment labeled "A" (for antigen) that has a place for every atom that is present in the topology file.
set AgPDB "off_center_antigen.pdb"
segment A {pdb $AgPDB}
# Load the coordinates of the atoms that are present in the structure file.
coordpdb $AgPDB A
# Guess coordinates for missing atoms.
guesscoord
# Write a PDB and PSF for the antigen with newly added atoms.
set addedAtomsPDB "AgAddedAtoms.pdb"
set AgPSF "Ag.psf"
writepdb $addedAtomsPDB
writepsf $AgPSF

exit
