# This VMD script tests all of the combinations of levels in the OptMAVEn grid and outputs the combinations that do not cause clashes between the antigen and both prototype antibodies.
# The following arguments are required:
# -e: the name of this script.
# -m: 0: antigen/coords.pdb 
#     1: MoleculeH/coords.pdb
#     2: MoleculeK/coords.pdb
# -args: 0: experiment/details.txt
#        1: output/file.dat
#        2: clashCutoff
#        3: clashesPermitted

set InstallFolder "/gpfs/work/m/mfa5147/OptMAVEn2.0"

# Load the VMD functions.
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Select every atom in each molecule.
set Ag [atomselect 0 "all"]
set IgH [atomselect 1 "all"]
set IgK [atomselect 2 "all"]

# Make sure each molecule exists.
if {[llength [$Ag get name]] == 0 || [llength [$IgH get name]] == 0 || [llength [$IgK get name]] == 0} {
	puts "One or more molecules failed to load."
	exit 1
}

# Load the arguments.
if {[llength $argv] != 4} {
	puts "Usage: -args experiment/details.txt output/file.dat clashCutoff"
	exit 1
}
set detailsFile [lindex $argv 0]
set positionsFile [lindex $argv 1]
set clashCutoff [lindex $argv 2]
set clashesPermitted [lindex $argv 3]

# Select the epitope in the antigen.
set epitopeResidues [readAttributeList $detailsFile "Epitope Position" 1]
set epitope [atomselect 0 "resid $epitopeResidues"]
if {[llength [$epitope get name]] == 0} {
	puts "Failed to load epitope from $detailsFile."
	exit 1
}

# Load z rotations and x, y, and z translations.
set levels []
for {set i 1} {$i < 5} {incr i} {
	set dimlevels [split [readAttribute $detailsFile "Optmaven Grid" $i]]
	puts dimlevels
	puts $dimlevels
	if {[llength $dimlevels] == 0} {
		puts "Experiment details file is missing one or more grid levels."
		exit 1
	}
	lappend levels $dimlevels
}

# Calculate the center of geometry of the antigen and epitope.
set AgCenter [measure center $Ag]
set epiCenter [measure center $epitope]
set faceAngle 0

# Determine the maximum z coordinate of any atom in either antibody.
set maxZH [lindex [measure minmax $IgH] {1 2}]
set maxZK [lindex [measure minmax $IgK] {1 2}]
set maxZIg [expr max($maxZH, $maxZK)]

# Determine the minimum z coordinate of any atom in the antigen.
set minZAg [lindex [measure minmax $IgH] {0 2}]

# Test all positions and record those that do not clash.
set positions [open $positionsFile "w"]
foreach zRot [lindex $levels 0] {
	# Rotate around the z axis to point the antigen in the direction of zRot.
	rotateZAxis $Ag $AgCenter [expr $zRot - $faceAngle]
	# Update the facing angle.
	set faceAngle $zRot
	foreach z [lindex $levels 3] {
		# Select all of the atoms in the antigen with z coordinates within the cutoff z coordinate of either antibody.
		set AgNear [atomselect 0 "z < $maxZIg + $clashCutoff"]
		# Select the atoms in the antibody with z coordinates within the cutoff z coordinates of the antigen.
		set IgHNear [atomselect 1 "z > $minZAg - $clashCutoff"]
		set IgKNear [atomselect 2 "z > $minZAg - $clashCutoff"]
		foreach y [lindex $levels 2] {
			foreach x [lindex $levels 1] {
				# Define the point P to which to move the center of the epitope.
				set P "$x $y $z"
				puts "$zRot $P"
				# Calculate the translation needed to move the antigen so that the center of geometry of its epitope lies at P.
				set translate [vecsub $P $epiCenter]
				# Move the antigen.
				$Ag moveby $translate
				# Update the centers.
				set AgCenter [vecadd $AgCenter $translate]
				set epiCenter [vecadd $epiCenter $translate]
				# Update the minimum z coordinate of the antigen.
				set translateZ [lindex $translate 2]
				set minZAg [expr $minZAg + $translateZ]
				# Count the number of clashes between the antigen and the antibody.
				set clashes [llength [lindex [measure contacts $clashCutoff $AgNear $IgHNear] 0]]
				if {$clashes <= $clashesPermitted} {
					set clashes [llength [lindex [measure contacts $clashCutoff $AgNear $IgKNear] 0]]
					if {$clashes <= $clashesPermitted} {
						# Write non-clashing positions to the file.
						puts $positions "$zRot $P"
						# $Ag writepdb "zr$zRot.x$x.y$y.z$z.pdb"  ;# To save disk space, only write PDBs while benchmarking.
					}
				}
			}
		}
	}
}
close $positions

exit 0
