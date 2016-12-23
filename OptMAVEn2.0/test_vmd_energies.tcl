# Load the VMD functions.
package require namdenergy
source VMD_FUNCTIONS.tcl

# Load the antigen and the MAPs part.
mol new AgIg.psf
mol addfile AgIg.pdb

# Select every atom in each segment.
set Ag [atomselect 0 "segname A"]
set MAPsPart [atomselect 0 "segname I"]

# Select the epitope in the antigen.
set epitopeResidues [gets [open epitope.dat "r"]]
set epitope [atomselect 0 "resid $epitopeResidues and segname A"]

# Define the initial position of the antigen, which is zero degrees and point (0, 0, 0).
set faceAngle 0
set pos "0 0 0"

# Calculate the energy at all positiions.
set positions [open "energy_positions.dat"]
while {[gets $positions P] >= 0} {
	# Get the difference between the needed and current coordinates.
	set dP [vecsub $P "$faceAngle $pos"]
	set dzRot [lindex $dP 0]
	set dx [lindex $dP 1]
	set dy [lindex $dP 2]
	set dz [lindex $dP 3]
	set dxyz "$dx $dy $dz"
	# Rotate the antigen to the specified direction.
	if {$dzRot != 0} {
		rotateZAxis $Ag $pos $dzRot
		set faceAngle [expr $faceAngle + $dzRot]
	}
	# Translate the antigen by the amount needed.
	$Ag moveby $dxyz
	set pos [vecadd $pos $dxyz]
	# Ensure correct position.
	puts Position
	puts $pos
	puts [getAntigenPosition $Ag $epitope 0]
	# Calculate the interaction energy (electrostatic and van der Waals) between the antigen and the MAPs part.
	namdenergy -sel $Ag $MAPsPart -elec -vdw -par par_all27_prot_na.prm -ofile "namdenergy$faceAngle$dx$dy$dz.dat"
}
close $positions

exit
