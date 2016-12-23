proc pi {} {
	return 3.14159265359
}

# Convert radians to degrees.
proc rad2deg {rad} {
	return [expr $rad * 180 / [pi]]
}

# Convert degrees to radians.
proc deg2rad {deg} {
	return [expr $deg * [pi] / 180]
}

# Calculate the angle (in radians, from 0 to pi) between two vectors (as of now, there is no single function to do this in VMD).
proc vectorAngle {v1 v2} {
	# Use the formula angle = acos((v1 dot v2) / (|v1||v2|)).
	return [expr acos([vecdot $v1 $v2] / ([veclength $v1] * [veclength $v2]))]
}

# Compute the angle (in radians, from 0 to 2pi) between the x axis and the projection of a vector onto the x axis.
proc angleToX {vec} {
	# Obtain the x and y coordinates of the vector.
	set xCoor [lindex $vec 0]
	set yCoor [lindex $vec 1]
	# Calculate the angle between the vector's projection onto the x-y plane and the x axis.
	set rotAngle [expr acos($xCoor / sqrt($xCoor * $xCoor + $yCoor * $yCoor))]
	if {$yCoor < 0} {
		set rotAngle [expr (2 * [pi]) - $rotAngle];  # If the y coordinate is negative, subtract the angle from 2pi radians.
	}
	return $rotAngle
}

# Calculate the transformation matrix needed to minimize the z coordinates of the epitope.
proc transformationZMin {antigen epitope} {
	# Calculate the "arm" vector from the center of geometry of the antigen to that of the epitope.
	set AgCenter [measure center $antigen]
	set arm [vecsub [measure center $epitope] $AgCenter]
	# Create the matrix needed to rotate around the axis formed by the cross product of the arm vector and the desired direction vector (0, 0, -1) by the angle between those two vectors.
	return [trans origin $AgCenter axis [veccross $arm "0 0 -1"] [vectorAngle $arm "0 0 -1"] rad]
}

# Calculate the z angle (in degrees) at which an antigen is facing. This angle is defined arbitrarily as the angle between the positive x axis and the projection onto the x-y plane of the vector from the antigen's center of mass to the C-alpha atom of its first residue, when the antigen has been rotated such that the z coordinates of its epitope are minimized. This definition is in place to provide a standard way to calculate the z rotation of any antigen structure.
proc getZAngle {antigen epitope antigenMolID} {
	# Transform the coordinates of the C-alpha atom of the first residue with the transformation matrix needed to minimize the z coordinates of the epitope and center the antigen at the origin. Then find the angle its x-y projection makes with the x axis.
	return [rad2deg [angleToX [coordtrans [transformationZMin $antigen $epitope] [measure center [atomselect $antigenMolID "resid 1 and name CA"]]]]]
}

# Rotate the antigen around the z axis by zRot degrees, ASSUMING that the epitope has had its z coordinates minimized.
proc rotateZAxis {antigen AgCenter zRot} {
	$antigen moveby [vecscale -1 $AgCenter]
        $antigen move [transaxis z $zRot deg]
        $antigen moveby $AgCenter
}

# Rotate the antigen's epitope to minimize its z coordinates.
proc rotateZMin {antigen epitope} {
	$antigen move [transformationZMin $antigen $epitope]
}

# Calculate the z rotation angle and (x, y, z) position of an antigen's epitope. NOTE: this function should only be used for antigens whose epitopes are already aligned with the negative z axis.
proc getAntigenPosition {antigen epitope antigenMolID} { 
	return "[getZAngle $antigen $epitope $antigenMolID] [measure center $epitope]"
}

# Move the antigen to the position with the epitope centered at the given (x, y, z) point and aligned with the negative z axis, and with the antigen facing the direction (in degrees) given by zRot.
proc positionAntigen {antigen epitope antigenMolID zRot x y z} {
	# Rotate the epitope to minimize its z coordinates.
	rotateZMin $antigen $epitope
	# Translate the antigen so that its epitope is centered at (x, y, z).
	$antigen moveby [vecsub "$x $y $z" [measure center $epitope]]
	# Rotate the antigen around the z axis so that it points in the direction given by zRot.
	$antigen move [transaxis z [expr $zRot - [getZAngle $antigen $epitope $antigenMolID]] deg]
}

# Move the antigen to the position with the epitope centered at the origin and pointing towards the negative z axis and with the antigen facing zero degrees.
proc mountAntigen {antigen epitope antigenMolID} {
	positionAntigen $antigen $epitope $antigenMolID 0 0 0 0
}

# Generate a range of numbers.
proc range {start end step} {
	set out {}
	set iMax [expr {1 + ((abs($end - $start) - 1) / abs($step))}]
	for {set i 0} {$i < $iMax} {incr i} {
		lappend out [expr {$start + ($i * $step)}]
	}
	return $out
}
