# Merge an antigen and a MAPs part into a single molecule and generate a PDB and PSF.
# The following command-line arguments are required:
# -e: the name of this script
# -args: 0: the experiment's directory
#        1: the coordinates of the antigen
#        2: the coordinates of the MAPs part
#        3: the prefix of the combined structure
#        4: the segment of the antigen
#        5: the segment of the MAPs part
#        6: the first topology file
#        7: (optional) the second topology file
#        etc.

# The psfgen package is required.
package require psfgen

# Load the VMD functions.
set InstallFolder "/home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0"
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] < 7} {
	puts "Usage: -args experiment/directory antigen/coordinates.pdb MAPs/part/coordinates.pdb combined/structure/prefix antigenSegment MAPsSegment topology/file1.rtf (optional) topology/file2.rtf ..."
	exit 1
}
set expDir [lindex $argv 0]
set AgPDB [lindex $argv 1]
set MAPsPDB [lindex $argv 2]
set combined [lindex $argv 3]
set antigenSegment [lindex $argv 4]
set MAPsSegment [lindex $argv 5]

# Read the topology file and alias the atoms.
for {set i 6} {$i < [llength $argv]} {incr i} {
	topology [lindex $argv $i]
}
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# Build a segment for the antigen.
segment $antigenSegment {
	pdb $AgPDB
}

# Read antigen coordinates into the antigen segment.
coordpdb $AgPDB $antigenSegment

# Build a segment for the MAPs part. Do not patch the terminal residues.
segment $MAPsSegment {
    first none
    last none
	pdb $MAPsPDB
}

# Read the immunoglobulin coordinates into the MAPs segment.
coordpdb $MAPsPDB $MAPsSegment

# Create the PDB and PSF.
writepdb "$combined.pdb"
writepsf "$combined.psf"

exit 0
