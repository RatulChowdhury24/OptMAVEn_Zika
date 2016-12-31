# Merge an antigen and a MAPs part into a single molecule and generate a PDB and PSF.
# The following command-line arguments are required:
# -e: the name of this script
# -args: 0: the experiment's directory
#        1: the coordinates of the antigen
#        2: the coordinates of the MAPs part
#        3: the prefix of the combined structure
#        4: the first topology file
#        5: (optional) the second topology file
#        etc.

# The psfgen package is required.
package require psfgen

# Load the VMD functions.
set InstallFolder "/gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0"
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] < 5} {
	puts "Usage: -args experiment/directory antigen/coordinates.pdb MAPs/part/coordinates.pdb combined/structure/prefix topology/file1.rtf (optional) topology/file2.rtf ..."
	exit 1
}
set expDir [lindex $argv 0]
set AgPDB [lindex $argv 1]
set MAPsPDB [lindex $argv 2]
set combined [lindex $argv 3]

# Read the topology file and alias the atoms.
for {set i 4} {$i < [llength $argv]} {incr i} {
	topology [lindex $argv $i]
}
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# Build "segment A" for the antigen.
segment A {
	pdb $AgPDB
}

# Read antigen coordinates into segment A.
coordpdb $AgPDB A

# Build "segment I" for the immunoglobulin part.
segment I {
	pdb $MAPsPDB
}

# Read the immunoglobulin coordinates into segment I.
coordpdb $MAPsPDB I

# Create the PDB and PSF.
writepdb "$combined.pdb"
writepsf "$combined.psf"

exit 0
