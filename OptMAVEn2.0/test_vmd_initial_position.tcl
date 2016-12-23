# This script rotates the relaxed antigen so as to minimize the z coordinates of its epitope, centers the epitope at the origin, and rotates the antigen so it points in the zero degree direction.

# Load the VMD functions.
source VMD_FUNCTIONS.tcl

# Select the entire antigen and the epitope.
set Ag [atomselect 0 "all"]
set epitope [atomselect 0 "resid 306 307 342 343 344 350 353 355 388 391 392 393 395"]

# Move the antigen to the initial position, with the epitope's z coordinates minimized, the epitope centered at the origin, and the antigen facing the zero degree direction. Save a structure of the initial antigen position.
mountAntigen $Ag $epitope 0
$Ag writepdb "AgInit.pdb"

exit
