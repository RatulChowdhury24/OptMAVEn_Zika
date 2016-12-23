# Load the VMD functions.
source VMD_FUNCTIONS.tcl

# Load the molecules.
mol new AgInit.pdb;  # molecule 0
mol new MoleculeH.pdb;  # molecule 1
mol new MoleculeK.pdb;  # molecule 2

# Select every atom in each molecule.
set Ag [atomselect 0 "all"]
set IgH [atomselect 1 "all"]
set IgK [atomselect 2 "all"]

# Select the epitope in the antigen.
set epitopeResidues [gets [open epitope.dat "r"]]
set epitope [atomselect 0 "resid $epitopeResidues"]

# Calculate the center of geometry of the antigen, epitope, and the alpha carbon of the antigen's first residue.
set AgCenter [measure center $Ag]
set epiCenter [measure center $epitope]
set CA1Center [measure center [atomselect 0 "resid 1 and name CA"]]
set faceAngle 0

# Load z rotations and x, y, and z translations.
set levels []
set dimensions [split [read [open "grid_levels.dat"]] "\n"]
for {set i 0} {$i < 4} {incr i} {
	lappend levels [split [lindex $dimensions $i]]
}

# Define the distance within which atoms are considered to clash.
set clashCutoff 1

# Determine the maximum z coordinate of any atom in either antibody.
set maxZH [lindex [measure minmax $IgH] {1 2}]
set maxZK [lindex [measure minmax $IgK] {1 2}]
set maxZIg [expr max($maxZH, $maxZK)]

# Determine the minimum z coordinate of any atom in the antigen.
set minZAg [lindex [measure minmax $IgH] {0 2}]

# Test all positions and record those that do not clash.
set positionsFile [open "energy_positions.dat" "w"]
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
				if {$clashes == 0} {
					set clashes [llength [lindex [measure contacts $clashCutoff $AgNear $IgKNear] 0]]
					if {$clashes == 0} {
						# Write non-clashing positions to the file.
						$Ag writepdb "clashtest$zRot$x$y$z.pdb"
						puts $positionsFile "$zRot $P"
					}
				}
			}
		}
	}
}
close $positionsFile

exit
