# Assemble an antigen structure (with a provided position) and six MAPs parts into an antibody/antigen complex. Generate a PDB and a PSF.
# The following command-line arguments are required:
# -e: the name of this script
# -args: 0: the experiment's directory
#        1: the experiment details file
#        2: the antigen zRot
#        3: the antigen x position
#        4: the antigen y position
#        5: the antigen z position
#        6: the coordinates of the antigen
#        7: the coordinates of the HV part
#        8: the coordinates of the HJ part
#        9: the coordinates of the HCDR3 part
#       10: the coordinates of the LV part
#       11: the coordinates of the LJ part
#       12: the coordinates of the LCDR3 part
#       13: the prefix of the assembled structure
#       14: the segment of the antigen
#       15: the segment of the antibody heavy chain
#       16: the segment of the antibody light chain
#       17: the first topology file
#       18: (optional) the second topology file
#        etc.

# The psfgen package is required.
package require psfgen

# Load the VMD functions.
set InstallFolder "/gpfs/work/m/mfa5147/OptMAVEn2.0"
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] < 17} {
	puts "Usage: -args experiment/directory antigen/coordinates.pdb MAPs/part/coordinates.pdb combined/structure/prefix antigenSegment MAPsSegment topology/file1.rtf (optional) topology/file2.rtf ..."
	exit 1
}
set expDir [lindex $argv 0]
set detailsFile [lindex $argv 1]
set zRot [lindex $argv 2]
set x [lindex $argv 3]
set y [lindex $argv 4]
set z [lindex $argv 5]
set AgPDB [lindex $argv 6]
set HV [lindex $argv 7]
set HJ [lindex $argv 8]
set HCDR3 [lindex $argv 9]
set LV [lindex $argv 10]
set LJ [lindex $argv 11]
set LCDR3 [lindex $argv 12]
set assembly [lindex $argv 13]
set antigenSegment [lindex $argv 14]
set AbHSegment [lindex $argv 15]
set AbLSegment [lindex $argv 16]
# Read the topology file and alias the atoms.
for {set i 17} {$i < [llength $argv]} {incr i} {
	topology [lindex $argv $i]
}
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD


### First, transform the antigen to position it in the correct place.
mol new $AgPDB
# Read the position.
set antigenMolID 0
# Select the entire antigen.
set Ag [atomselect $antigenMolID "all"]
# Select the epitope in the antigen.
#set epitopeResidues [readAttributeList $detailsFile "Epitope Position" 1]
#set epitope [atomselect $antigenMolID "resid $epitopeResidues"]
#if {[llength [$epitope get name]] == 0} {
#	puts "Failed to load epitope from $detailsFile."
#	exit 1
#}
# Position the antigen and save the new coordinates to a temporary file.
set faceAngle 0
set pos "0 0 0"
set P "$zRot $x $y $z"
repositionAntigen $Ag $P $faceAngle $pos
#puts "EPITOPE:"
#puts $epitopeResidues
set transformedAgPDB "$AgPDB.temp"
puts [getAntigenPosition $Ag $epitope $antigenMolID $antigenSegment]
$Ag writepdb $transformedAgPDB


# Build a segment for the antigen.
segment $antigenSegment {
	pdb $transformedAgPDB
}

# Read the transformed antigen coordinates into the antigen segment and delete the temporary file.
coordpdb $transformedAgPDB $antigenSegment
#FIXME file delete $transformedAgPDB

# Build a segment for the antibody heavy chain.
segment $AbHSegment {
    pdb $HV
    pdb $HCDR3
    pdb $HJ
}

# Read the MAPs parts coordinates into the antibody.
coordpdb $HV $AbHSegment
coordpdb $HCDR3 $AbHSegment
coordpdb $HJ $AbHSegment

# Build a segment for the antibody light chain.
segment $AbLSegment {
    pdb $LV
    pdb $LCDR3
    pdb $LJ
}

# Read the MAPs parts coordinates into the antibody.
coordpdb $LV $AbLSegment
coordpdb $LCDR3 $AbLSegment
coordpdb $LJ $AbLSegment

# Create the PDB and PSF.
writepdb "$assembly.pdb"
writepsf "$assembly.psf"

exit 0
