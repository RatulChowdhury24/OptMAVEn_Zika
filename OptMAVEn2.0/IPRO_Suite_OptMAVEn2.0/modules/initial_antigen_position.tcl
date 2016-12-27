# This script rotates the relaxed antigen so as to minimize the z coordinates of its epitope, centers the epitope at the origin, and rotates the antigen so it points in the zero degree direction.
# The following arguments are required:
# -e: the name of this script
# -f: the file name of the antigen
# -args: 0: the file name of the output antigen
#        1: the name of the experiment details file

# Load the VMD functions.
set InstallFolder "/home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0"
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Get the arguments.
if {[llength $argv] != 2} {
	puts "Usage: -f input/coordinates.pdb -args output/coordinates.pdb experiment/details.txt"
	exit 1
}
set outFile [lindex $argv 0]
set epiRes [readAttributeList [lindex $argv 1] "Epitope Position" 1]

# Select the entire antigen and the epitope.
set Ag [atomselect 0 "all"]
set epitope [atomselect 0 "resid $epiRes"]

# Move the antigen to the initial position, with the epitope's z coordinates minimized, the epitope centered at the origin, and the antigen facing the zero degree direction. Save a structure of the initial antigen position.
mountAntigen $Ag $epitope 0
$Ag writepdb $outFile

exit 0
