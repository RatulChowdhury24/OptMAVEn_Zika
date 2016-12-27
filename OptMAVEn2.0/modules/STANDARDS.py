#!/usr/bin/env python

# The name of the file
__name__ = "IPRO Suite Standards Module"
# The documenation string for the module
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains standards and defaults used throughout the IPRO Suite of 
Programs."""

import os
import sys

# Make sure the modules know where they are installed
#InstallFolder = "/home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0"
InstallFolder = "/home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0"

# Where the CPLEX solver is located
CPLEXFolder = "/usr/global/ilog/CPLEX_Studio124/cplex/python/x86-64_sles10_4.1/"

# The standard user of this installation of the IPRO suite
defaultUser = "rzc158"

# The programs included in the IPRO suite
supportedPrograms = ['IPRO', 'OptGraft', 'OptCDR', 'OptZyme', "Mutator", "OptMAVEn"]

# The supported force fields
supportedFields = ["CHARMM"]
# The default force field
defaultField = "CHARMM"

# The list of all molecule file formats that are supported by the IPRO suite
supportedFormats = ['PDB']
# The default format
defaultFormat = "PDB"

# The default Lennard-Jones energy potential softening term
defaultUseSoftening = True 
defaultLjPhi = 0.005 

# The names of amino acids, organized by file format. Note that for PDB formats,
# there are both HIS and HSD for histidine, as HSD is the preferred protonation
# state in CHARMM
aminoAcids = {
"PDB":['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'HSD', 'ILE', 'LYS', \
'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
}
# Conversions between 1 and 3 letter amino acid names, in alphabetical order
convertAA = {
"PDB":{"ALA":"A", "A":"ALA", "CYS":"C", "C":"CYS", "ASP":"D", "D":"ASP", \
       "GLU":"E", "E":"GLU", "PHE":"F", "F":"PHE", "GLY":"G", "G":"GLY", \
"HIS":"H", "HSD":"H", "H":"HIS", "ILE":"I", "I":"ILE", "LYS":"K", "K":"LYS", \
       "LEU":"L", "L":"LEU", "MET":"M", "M":"MET", "ASN":"N", "N":"ASN", \
       "PRO":"P", "P":"PRO", "GLN":"Q", "Q":"GLN", "ARG":"R", "R":"ARG", \
       "SER":"S", "S":"SER", "THR":"T", "T":"THR", "VAL":"V", "V":"VAL", \
       "TRP":"W", "W":"TRP", "TYR":"Y", "Y":"TYR"}
}
# Group the amino acids by chemical properties
propertiesAA = {
"PDB":{"Aliphatic":['ALA', 'GLY', 'ILE', 'LEU', 'MET', 'PRO', 'VAL'], \
       "Aromatic":['PHE', 'TRP', 'TYR'], \
       "Charged":['ASP', 'GLU', 'HIS', 'LYS', 'ARG'], \
       "omp":['ALA', 'GLY', 'ILE', 'LEU', 'MET', 'PRO', 'VAL','PHE', 'TRP']}
}
# The names of atoms that can be in the backbone of amino acids, organized by
# file format
backboneAtoms = {
"PDB":['N', 'HN', 'HN1', 'HN2', 'HT1', 'HT2', 'HT3', 'CA', 'HA', 'HA2', 'C', \
'O', 'OT1', 'OT2']
}

# Defaults associated with the CHARMM force field
# The maximum number of iterations each energy minimization will run for
defaultCHARMMIterations = 5000
# The energy terms CHARMM calculates
defaultCHARMMEnergies = ['angl', 'bond', 'dihe', 'elec', 'impr', 'urey', 'vdw']
# The topology files
defaultCHARMMTopologies = ['top_all27_prot_na.rtf']
# Parameter files
defaultCHARMMParameters = ['par_all27_prot_na.prm']
# The command to access CHARMM
CHARMMCommand = "/gpfs/group/cdm/c34b1.xj.gnu/exec/gnu/charmm.serial.xlarge"

# The command to run GAMS
GAMSCommand = "gams"

# Information about how to use implicit solvation
# Whether or not to use solvation by default
defaultUseSolvation = True
# The types of solvation that are supported for non-bonded energy calculations
# NOT in a force field
supportedSolvations = ["Lazaridis-Karplus"]
# The default type of implicit solvation to use for these calculations
defaultSolvationType = "Lazaridis-Karplus"
# The default list of LK solvation files
defaultLKSolvationFiles = ["solvation.dat"]

# The location of the rotamer library
RotamerLibraryPath = "/gpfs/group/cdm/rotlib/"
# Because of memory constraints there is a maximum number of rotamers that an
# Optimal Rotamer Selection protocol can choose between.
defaultMaxRotamers = 1300
# To help in making sure that the maximum rotamer number is not exceeded, only
# rotamers that are within "window" kcal/mol of the lowest energy rotamer at
# their position are permitted.
defaultRotamerWindow = 40.0

# Docking information
# How often docking should be run
defaultDockingFrequency = 3
# How long to run docking for
defaultDockingIterations = 500
# The standard deviation of moving molecules during docking
defaultSDMove = 0.2
# And of rotating molecules (in degrees)
defaultSDRotate = 2.0
# The simulated annealing temperature at the start of docking. See the next
# block of standards for a definition of what this means. This temperature
# retains 25% of molecule postions within 10 kcal/mol of the best.
defaultDockTempStart = 3640.0
# And at the end of docking (10% within 10 kcal/mol of the best)
defaultDockTempEnd = 2190.0

# Information related to running IPRO in any IPRO suite program
# The number of IPRO iterations to run
defaultIterations = 3000
# The gas constant, R, in kcal/mol
GasConstant = 0.001986
# The temperature parameter for simulated annealing, which has units of Kelvin,
# is calculated as follows:
# T = -(delta energy) / (GasConstant * LN(% of designs retained))
# 3640 = (-10.0)/(0.001986 * LN(0.25))
defaultAnnealingTemp = 3640.0
# If multiple processors are working on this job simultaneously, should they
# share simulated annealing results?
defaultShareAnnealing = True
# The default type of energy calculations to use in IPRO
defaultEnergyCalc = 'Interaction'

# Information related to which Residues should be assigned rotamers after
# backbone perturbations
# Should Residues be assigned rotamers if they are close in linear sequence
# (this is ALWAYS done to accomodate the perturbations) or based on having close
# distances in 3D space
defaultPackingMethod = "Distance"
# If repacking is being done based on Distance, do other Residues need to be
# close to the "RESTRAINED" residues or either restrained and FREE Residues?
defaultPackingSelection = "RESTRAINED"
# What is close? If a Residue has any non-hydrogen atom within this distance of
# any non-hydrogen atom in the packing selection, it will be repacked
defaultPackingCutoff = 4.5

# Whether or not to do refinements by default
defaultDoRefinement = True
# And how many iterations they should run for
defaultRefinementIterations = 25
# The number of ensembles to generate while doing them
defaultEnsembleNumber = 10

# OptMAVEn information
# Default human 9mer sequences file for deimmunization 
defaultHumanSequencesFile = "Human_9mer_Sequences.txt"
# Default Maps database integer cuts file for clashing
defaultIntegerCutsFile = "MAPs_Integer_Cuts.txt"

# Default OptZyme information
defaultWeight = 0.03957

# Spacing for printing things to a screen
SPACES = 80

def list_items(items):
    """Generate a grammatically correct list of strings."""
    # Store the formatted text here
    text = ''
    # If there are no items, say that
    if len(items) == 0:
        text += " NONE"
    # If there is only one item, just include it
    elif len(items) == 1:
        text += " " + items[0]
    # Use an 'and' between items when there are two
    elif len(items) == 2:
        text += " " + items[0] + " and " + items[1]
    # Otherwise use commas
    else:
        for i in range(len(items)):
            if i != len(items) - 1:
                text += " " + items[i] + ","
            else:
                text += " and " + items[i]
    return text

def screen_formatting(text, length = SPACES):
    """Format a string of text so that it conveniently fits on a screen."""
    # Store the output of this function in this string
    output = ''
    # Make sure the text is a string
    if not isinstance(text, str):
        text = str(text)
    # Split the text on end line characters so they are preserved in the output
    lines = text.split("\n")
    # Loop through lines
    for line in lines:
        # If it is a string of an Atom, include it as is even if it exceeds the
        # length input
        if line.startswith(("ATOM", "HETATM")):
            output += line + "\n"
        # Otherwise, include the contents of the line while obeying the length
        # requirement
        else:
            # Keep track of where the last white space (or the beginning of the
            # line) was encountered
            lastIndex = 0
            # Store the phrase being assembled here
            phrase = ''
            # Go through the characters in the line
            for i in range(len(line)):
                # If this is the last character or the next character is white
                # space
                if i == len(line) - 1 or line[i+1] in [' ', '\t']:
                    # Get the word that is being included
                    word = line[lastIndex:i+1]
                    # If adding this to the phrase makes the phrase too long
                    if len(word) + len(phrase) > length:
                        output += phrase + "\n"
                        phrase = ''
                        # Increment the last index by 1 (to remove the extra
                        # space) and restate what the word is
                        lastIndex += 1
                        word = line[lastIndex:i+1]
                    # Add the word to the phrase
                    phrase += word
                    # Set the last index appropriately
                    lastIndex = i+1
            # Add the phrase to the output
            output += phrase + "\n"
    return output[:-1]

class IPRO_Error(Exception):
    """An error for problems in the IPRO suite of programs."""
    def __init__(self, error = '', cancelJobs = False):
        """The initialization routine of the IPRO_Error class."""
        # Store the error
        if isinstance(error, str):
            self.error = error
        else:
            self.error = ''
        if cancelJobs:
            os.system("qstat > remove_jobs.out")
            jobs = []
            with open("remove_jobs.out") as file:
                for line in file:
                    items = line.split()
                    if items[0] == "Job" or items[0][0] == "-":
                        continue
                    elif items[2] == defaultUser:     
                        jobs.append(items[0].split(".")[0])
            for job in jobs:
                os.system("qdel " + job)
                             
    def __str__(self):
        """Create a formatted string describing the error."""
        # If there is no content in the error, return an empty string
        if not isinstance(self.error, str) or len(self.error) == 0:
            return ''
        # Otherwise use the screen_formatting function
        else:
            return "\n" + screen_formatting(self.error, SPACES)
