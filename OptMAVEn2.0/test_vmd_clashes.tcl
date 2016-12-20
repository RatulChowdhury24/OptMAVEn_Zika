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

# Load the molecules.
mol new AgInit.pdb;  # molecule 0
mol new MoleculeH.pdb;  # molecule 1

# Select every atom in each molecule.
set Ag [atomselect 0 "all"]
set IgH [atomselect 1 "all"]

# Select the epitope.
set epitope [atomselect 0 "resid 306 307 342 343 344 350 353 355 388 391 392 393 395"]

# Calculate the center of geometry of the epitope.
set epiCenter [measure center $epitope]

# Define x, y, and z translations.
set xt [range -5 10 1.25]
set yt [range -10 5 1.25]
set zt [range -16 -3 1]

# Test all positions.
set count 0;
foreach x $xt {
	foreach y $yt {
		foreach z $zt {
			# Define the point to which to move the center of the epitope.
			set P "$x $y $z"; #lappend $P $x; lappend $P $y; lappend $P $z;
			# Calculate the translation needed to move the antigen so that the center of geometry of its epitope lies at P.
			set trans [vecsub $P $epiCenter]
			# Move the antigen.
			$Ag moveby $trans
			# Update the center of the epitope.
			set epiCenter [vecadd $epiCenter $trans]
			# Count the number of clashes (atoms closer than 1.5 A) between the antigen and the antibody.
			#set clashes [llength [lindex [measure contacts 1.5 $Ag $IgH] 0]]
			set clashes [llength [exwithin 1 of $IgH]]
			if {$clashes == 0} {
				puts $P
				$Ag writepdb "noclash$count.pdb"
			}
			incr count;
		}
	}
}

#namdenergy

exit
