# Generate a range of numbers.
proc range {start end step} {
	set out {}
	set iMax [expr {1 + ((abs($end - $start) - 1) / abs($step))}]
	for {set i 0} {$i < $iMax} {incr i} {
		lappend out [expr {$start + ($i * $step)}]
	}
	return $out
}

# Set working directory.
cd /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/modules

# Load the antigen.
mol new off_center_antigen.pdb;  # molecule 0

# Select the entire antigen and the epitope.
set Ag [atomselect 0 "all"]
set epitope [atomselect 0 "resid 306 307 342 343 344 350 353 355 388 391 392 393 395"]

# Move the center of geometry of the antigen to the origin.
set AgCenter [measure center $Ag]
$Ag moveby [vecscale -1 $AgCenter]

# Calculate the center of geometry of the epitope.
set epiCenter [measure center $epitope]

# Rotate the antigen so that the vector from the center of the antigen to the epitope points in the negative z direction (0, 0, -1), which minimizes the z coordinates of the epitope.
set rotAxis [veccross "0 0 -1" $epiCenter];  # The axis of rotation is the cross product of the final vector (0, 0, -1) with the initial vector (epiVect).
set normAxis [veclength $rotAxis]
set normEpi [veclength $epiCenter]
set rotAngle [expr -1 * asin($normAxis / $normEpi)];  # The angle (radians) between the vectors is computed using sin(angle) = (|a x b|) / (|a||b|), where a = (0, 0, -1), |a| = 1, b = epiVect, and a x b = axis.
# Make the rotation matrix.
set transfrm [trans axis $rotAxis $rotAngle rad]
$Ag move $transfrm

# Move the center of geometry of the epitope to the origin.
set epiCenter [measure center $epitope]
$Ag moveby [vecscale -1 $epiCenter]

# Save a structure of the mounted antigen.
$Ag writepdb "AgMount.pdb"
exit
