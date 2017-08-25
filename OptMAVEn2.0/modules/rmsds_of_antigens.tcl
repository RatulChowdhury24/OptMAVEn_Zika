# This VMD script calculates the RMSD between the atoms in the experimental antigen position and those in each position found during the clash culling.
# The following arguments are required:
# -e: the name of this script.
# -m: 0: experimental/antigen/coords.pdb
#     1: mounted/antigen/coords.pdb
# -args: 0: experiment/details.txt
#        1: positions/file.dat
#        2: output/file.dat

set InstallFolder "/gpfs/work/m/mfa5147/OptMAVEn2.0"

# Load the VMD functions.
source $InstallFolder/modules/VMD_FUNCTIONS.tcl

# Load the arguments.
if {[llength $argv] < 3} {
	puts "Usage: structure.psf coordinates.pdb antigenSegment MAPsSegment positions/file.dat output/energies.dat experiment/details.txt parameter/file1.prm parameter/file2.prm (optional) ..."
	exit 1
}
set detailsFile [lindex $argv 0]
set positionsFile [lindex $argv 1]
set outputFile [lindex $argv 2]

# Select the epitope in the antigen.
set epitopeResidues [readAttributeList $detailsFile "Epitope Position" 1]
set epitope [atomselect 0 "resid $epitopeResidues"]
if {[llength [$epitope get name]] == 0} {
	puts "Failed to load epitope from $detailsFile."
	exit 1
}

# Define the initial position of the antigen, which is zero degrees and point (0, 0, 0).
set faceAngle 0
set pos "0 0 0"

# Make a file of rmsds.
set rmsdFile [open $outEnergy "w"]
close $rmsdFile

# Calculate the RMSD at all positiions.
set positions [open $positionsFile]
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
	# FIXME
	puts " "
	puts $pos
	puts [getAntigenPosition $Ag $epitope 0 $antigenSegment]
	# Calculate the RMSD between the antigen's position and the experimental position.
	set rmsd 
	# Read the energy from the fourth position in the list.
	set energy [lindex $energyList 4]
	# Write the energy to the file of all energies.
	set eFile [open $outEnergy "a"]
	puts $eFile "$P $energy"
	close $eFile
}
close $positions

# Show that the calculations are finished.
set fin [open "finished" "w"]
close $fin

# Remove the PDB and PSF to save disk space.
file delete $inStruct
file delete $inCoords

exit 0
