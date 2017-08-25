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
        #FIXME
        puts "ANGLE TO X"
	# Obtain the x and y coordinates of the vector.
	set xCoor [lindex $vec 0]
	set yCoor [lindex $vec 1]
	puts "COORDS"
	puts $xCoor
	puts $yCoor
	# Calculate the angle between the vector's projection onto the x-y plane and the x axis.
	set rotAngle [expr acos($xCoor / sqrt($xCoor * $xCoor + $yCoor * $yCoor))]
	# Convert to degrees.
	set rotAngle [expr {$rotAngle * 180 / [pi]}]
	puts "ROT ANGLE"
	puts $rotAngle
	if {$yCoor < 0} {
		puts "subtracting"
		set rotAngle [expr {360 - $rotAngle}];  # If the y coordinate is negative, subtract the angle from a full rotation.
	}
	puts $rotAngle
	puts "FINISHED"
	return $rotAngle
}

# Calculate the transformation matrix needed to minimize the z coordinates of the epitope.
proc transformationZMin {antigen epitope} {
        puts "Z min transformation"
	# Calculate the "arm" vector from the center of geometry of the antigen to that of the epitope.
        puts "Measuring centers"
	set AgCenter [measure center $antigen]
	set arm [vecsub [measure center $epitope] $AgCenter]
	# Create the matrix needed to rotate around the axis formed by the cross product of the arm vector and the desired direction vector (0, 0, -1) by the angle between those two vectors.
	puts "Calculating matrices"
	if {[veclength $arm] > 0.001} {
		#FIXME
		puts "Returned matrix"
		puts [trans center $AgCenter offset $AgCenter axis [veccross $arm "0 0 -1"] [vectorAngle $arm "0 0 -1"] rad]
		return [trans center $AgCenter offset $AgCenter axis [veccross $arm "0 0 -1"] [vectorAngle $arm "0 0 -1"] rad]
	} else {
		puts "Return identity matrix"
		return { {1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1} }
	}
}

# Calculate the z angle (in degrees) at which an antigen is facing. This angle is defined arbitrarily as the angle between the positive x axis and the projection onto the x-y plane of the vector from the antigen's center of mass to the C-alpha atom of its first residue, when the antigen has been rotated such that the z coordinates of its epitope are minimized. This definition is in place to provide a standard way to calculate the z rotation of any antigen structure.
proc getZAngle {antigen epitope antigenMolID antigenSegment} {
	# Transform the coordinates of the C-alpha atom of the first residue with the transformation matrix needed to minimize the z coordinates of the epitope and center the antigen at the origin. Then find the angle its x-y projection makes with the x axis.
	puts "get Z angle"
	puts "Transformation"
	puts  [transformationZMin $antigen $epitope]
	puts "Untransformed"
	puts   [measure center [atomselect $antigenMolID "segid $antigenSegment and residue 0 and name CA"]]
	puts "Transformed"
	puts [coordtrans [transformationZMin $antigen $epitope] [measure center [atomselect $antigenMolID "segid $antigenSegment and residue 0 and name CA"]]]
	puts "angle"
	set transformation [transformationZMin $antigen $epitope]
	set transformedAntigen [coordtrans $transformation [measure center $antigen]]
	set transformedCAlpha [coordtrans $transformation [measure center [atomselect $antigenMolID "segid $antigenSegment and residue 0 and name CA"]]]
	#puts [angleToX [coordtrans [transformationZMin $antigen $epitope] [measure center [atomselect $antigenMolID "segid $antigenSegment and residue 0 and name CA"]]]]
	puts [angleToX [vecsub $transformedCAlpha $transformedAntigen]]
	puts "ZQMLES"
	#FIXME
	#return [rad2deg [angleToX [coordtrans [transformationZMin $antigen $epitope] [measure center [atomselect $antigenMolID "segid $antigenSegment and residue 0 and name CA"]]]]]
	return [angleToX [vecsub $transformedCAlpha $transformedAntigen]]
}

# Rotate the antigen around the z axis by zRot degrees, ASSUMING that the epitope has had its z coordinates minimized.
proc rotateZAxis {antigen AgCenter zRot} {
	#$antigen moveby [vecscale -1 $AgCenter]
	#$antigen move [transaxis z $zRot deg]
	#$antigen moveby $AgCenter
	$antigen move [trans center $AgCenter offset $AgCenter axis z $zRot deg]

}

# Rotate the antigen's epitope to minimize its z coordinates.
proc rotateZMin {antigen epitope} {
	$antigen move [transformationZMin $antigen $epitope]
}

# Calculate the z rotation angle and (x, y, z) position of an antigen's epitope. NOTE: this function should only be used for antigens whose epitopes are already aligned with the negative z axis.
proc getAntigenPosition {antigen epitope antigenMolID antigenSegment} { 
	return "[getZAngle $antigen $epitope $antigenMolID $antigenSegment] [measure center $epitope]"
}

# Move the antigen to the position with the epitope centered at the given (x, y, z) point and aligned with the negative z axis, and with the antigen facing the direction (in degrees) given by zRot.
proc positionAntigen {antigen epitope antigenMolID antigenSegment zRot x y z} {
	#FIXME
        #$antigen writepdb "positionAntigen0.pdb"
	# Rotate the epitope to minimize its z coordinates.
	puts "Rotate Z min"
	rotateZMin $antigen $epitope
	#$antigen writepdb "positionAntigen1.pdb"
	# Rotate the antigen around the z axis so that it points in the direction given by zRot.
	puts "Rotate Z axis"
        set AgCenter [measure center $antigen]
	#$antigen move [trans center $AgCenter offset $AgCenter axis z [expr $zRot - [getZAngle $antigen $epitope $antigenMolID $antigenSegment]] deg]
	rotateZAxis $antigen $AgCenter [expr $zRot - [getZAngle $antigen $epitope $antigenMolID $antigenSegment]] 
        #$antigen writepdb "positionAntigen2.pdb"
	# Translate the antigen so that its epitope is centered at (x, y, z).
	puts "Translate into position"
	$antigen moveby [vecsub "$x $y $z" [measure center $epitope]]
	#$antigen writepdb "positionAntigen3.pdb"
        puts [measure center $antigen]
}

# Move the antigen to a new position P, given its current position (faceAngle and pos). Assume that the epitope Z coordinates have already been minimized, so do not perform this operation.
proc repositionAntigen {Ag P faceAngle pos} {
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
}

# Move the antigen to the position with the epitope centered at the origin and pointing towards the negative z axis and with the antigen facing zero degrees.
proc mountAntigen {antigen epitope antigenMolID antigenSegment} {
	positionAntigen $antigen $epitope $antigenMolID $antigenSegment 0 0 0 0
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

# Read an attribute from an experiment.
proc readAttribute {fileName attribute number} {
	set f [open $fileName]
	set getting 0
	set count 0
	# Read each line of the Experiment Details file.
	while {$count < $number && $getting != -1} {
		set getting [gets $f line]
		# If the line has a colon, then it contains an attribute.
		set div [split $line ":"]
		if {[llength $div] > 1} {
			set attr [lindex $div 0]
			# Count the number of times the attribute has been found.
			if {$attr == $attribute} {
				incr count
			}
		}
	}
	# Close the file.
	close $f
	# Make sure the attribute was found.
	if {$count < $number} {
		puts "$attribute not found."
		exit 1
	}
	# If so, return the associated value.
	return [string trim [lindex $div 1]]
}

# Read an attribute list from an experiment.
proc readAttributeList {fileName attribute position} {
	set f [open $fileName]
	set getting 0
	set attrList []
	# Read each line of the Experiment Details file.
	while {$getting != -1} {
		set getting [gets $f line]
		# If the line has a colon, then it contains an attribute.
		set div [split $line ":"]
		if {[llength $div] > 1} {
			set attr [lindex $div 0]
			# Add the value to the list of attributes.
			if {$attr == $attribute} {
				# Use the value at the specified position.
				if {$position >= 0} {
					lappend attrList [lindex [split [string trim [lindex $div 1]]] $position]
				} else {
					# If no position was specified, use the whole value.
					lappend attrList [lindex $div 1]
				}
			}
		}
	}
	# Close the file.
	close $f
	# Ensure there are items in the list.
	if {[llength $attrList] == 0} {
		puts "$attribute not found."
		exit 1
	}
	# Return the values.
	return $attrList
}
