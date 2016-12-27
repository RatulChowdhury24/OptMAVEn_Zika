# This VMD script adds missing atoms to the initial structure of an antigen and generates the PSF that is necessary to perform energy calculations.
# The following command-line arguments are required:
# -e: the name of this script
# -args: 0: the input coordinate file (PDB) of the antigen
#        1: the output coordinate file (PDB)
#        2: the output structure file (PSF)
#        3 ... : the topology file(s)

# Load required external modules.
package require psfgen

# Load the arguments.
if {[llength $argv] < 4} {
	puts "Usage: -args input/coordinates.pdb output/coordinates.pdb output/structure.psf topology1.rtf topology2.rtf (optional) ..."
	exit 1
}

set inputCoords [lindex $argv 0]
set outputCoords [lindex $argv 1]
set outputStruct [lindex $argv 2]
for {set i 3} {$i < [llength $argv]} {incr i} {
	topology [lindex $argv $i]
}

# Generate a PSF of the antigen for future NAMD energy calculations.
# Alias atoms and residues to make the structure and the topology file(s) compatible.
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
# Generate a segment labeled "A" (for antigen) that has a place for every atom that is present in the topology file(s).
set AgPDB $inputCoords
segment A {pdb $AgPDB}
# Load the coordinates of the atoms that are present in the coordinate file.
coordpdb $AgPDB A
# Guess coordinates for missing atoms.
guesscoord
# Write a PDB and PSF for the antigen with newly added atoms.
writepdb $outputCoords
writepsf $outputStruct

exit 0
