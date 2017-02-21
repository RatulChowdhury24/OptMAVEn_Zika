# This script extracts specific chains from a coordinate file and saves them to a new coordinate file.
# The required arguments are as follows:
# -e: the name of this script
# -f: the name of the coordinate file
# -argv: 0: the name of the output coordinate file
#        1: the first chain to extract
#        2: (optional) the second chain to extract
#        etc.

# Load the arguments.
if {[llength $argv] < 2} {
	puts "Usage: -f input/coordinates.pdb -args output/coordinates.pdb chain1 chain2 (optional) ..."
	exit 1
}

set outputCoords [lindex $argv 0]

set chains []
for {set i 1} {$i < [llength $argv]} {incr i} {
	lappend chains [lindex $argv $i]
}

if {[llength $chains] < 1} {
	puts "At least one chain must be specified."
	exit 1
}

# Select all atoms in the specified chain(s).
set selection [atomselect 0 "chain $chains"]

# Make sure the selection contains atoms.
if {[llength [$selection get chain]] == 0} {
	puts "The selection of chain(s) $chains contains no atoms."
	exit 1
}

# Write the selection to a PDB.
$selection writepdb $outputCoords

exit 0
