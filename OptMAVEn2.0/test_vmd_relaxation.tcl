package require psfgen

# Open the configuration file that describes how to perform the energy minimization.
set configuration [open "relaxation.dat"]

# Read the lines of the configuration file.
set inPDB [gets $configuration]  # The first line is the name of the PDB file to minimize.
set selectionText [gets $configuration]  # The second line is how much of the structure to select (e.g. "all", "protein", etc.)
set outPDB [gets $configuration]  # The third line is the name of the output PDB file.
set outPSF [gets $configuration]  # The fourth line is the name of the output PSF file.
set chainName [gets $configuration]  # The fifth line is the name of the chain as it will appear in the output PDB file.
set topolFile [gets $configuration]  # The sixth line is the name of the topology file.

# Load the PDB and select the specified atoms.
mol new $inPDB
set selection [atomselect 0 $selectionText]  

# Create a PDB of the selection. This PDB will be overwritten after the energy minimization.
$selection writepdb $outPDB

# Load the topology file and make residue and atom name aliases to ensure compatibility between the topology and PDB files.
topology $topolFile
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

# Create a virtual segment that has all of the atoms but lacks coordinates.
segment $chainName {pdb $outPDB}

# Load the coordinates of the atoms that exist and guess the coordinates of those that do not.
coordpdb $outPDB $chainName
guesscoord

# Write the PDB (with the newly guessed coordinates) and PSF files.
writepdb $outPDB
writepsf $outPSF


