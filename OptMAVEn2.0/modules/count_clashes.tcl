# This VMD script counts steric clashes between two molecules. The following arguments are required:
# -e: the name of this script.
# -m: 0: molecule1.pdb
#     1: molecule2.pdb
# -args: 0: output/file.dat
#        1: clashCutoff

set InstallFolder "/gpfs/work/m/mfa5147/OptMAVEn2.0"

# Load the VMD functions.
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] != 2} {
	puts "Usage: -args output/file.dat clashCutoff"
	exit 1
}
set outputFile [lindex $argv 0]
set clashCutoff [lindex $argv 1]

# Select the atoms in the molecules.
set mol1 [atomselect 0 "all"]
set mol2 [atomselect 1 "all"]

# Count clashes between the molecules.
set clashes [llength [lindex [measure contacts $clashCutoff $mol1 $mol2] 0]]

# Open a file to write the number of clashes.
set results [open $outputFile "w"]
puts $results $clashes
close $results

exit 0
