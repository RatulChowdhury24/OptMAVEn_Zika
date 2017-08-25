# This VMD script calculates the interaction energy between two segments of a molecule.
# The following arguments are required:
# -e: the name of this script
# -args: 0: the PSF
#        1: the PDB
#        2: the first segment
#        3: the second segment
#        4: the name of the output file
#        5: the first parameter file
#        6: (optional) the second parameter file
#        etc.

# Load the VMD functions.
package require namdenergy
set InstallFolder "/gpfs/work/m/mfa5147/OptMAVEn2.0"
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] < 6} {
	puts "Usage: structure.psf coordinates.pdb segment1 segment2 output/energies.dat parameter/file1.prm parameter/file2.prm (optional) ..."
	exit 1
}
set inStruct [lindex $argv 0]
set inCoords [lindex $argv 1]
set seg2 [lindex $argv 2]
set seg1 [lindex $argv 3]
set outEnergy [lindex $argv 4]
set parameterFiles []
for {set i 5} {$i < [llength $argv]} {incr i} {
	lappend parameterFiles [lindex $argv $i]
}

# Load the molecule.
mol new $inStruct
mol addfile $inCoords

# Select every atom in each segment. These segments were named by merge_antigen_part.tcl.
set sel1 [atomselect 0 "segname $seg1"]
set sel2 [atomselect 0 "segname $seg2"]

# Calculate the energy at all positiions.
set energyList [lindex [namdenergy -sel $sel1 $sel2 -elec -vdw -par $parameterFiles] 0]
# Read the energy from the fourth position in the list.
set energy [lindex [split $energyList] 4]
puts "Energy: $energy"
# Write the energy to the file of all energies.
set eFile [open $outEnergy "w"]
puts $eFile "$energy"
close $eFile

exit 0
