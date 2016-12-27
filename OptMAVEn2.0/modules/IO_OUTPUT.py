#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Experiment Details Output Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This file contains functions for outputting the details needed to run an IPRO
Suite Experiment. That information is stored in an Experiment_Details.txt
file"""

# Include standard PYTHON modules
import os
import sys
# Include other I/O functions
import IO_FUNCTIONS as FUNCTIONS
import IO_VALIDATE as VALIDATE

def standard_format(attribute, value):
    """A standard format for outputting attributes to Experiment Details"""
    # Store the created text here
    text = ''
    # If the value is a list, call this function recursively on each item
    if isinstance(value, list):
        for item in value:
            text += standard_format(attribute[:-1], item)
    else:
        # Create the LHS of the string
        lhs = attribute + ": "
        text += lhs.ljust(30)
        # And include options for how to format the value, based on what type of
        # variable it is
        if isinstance(value, float):
            text += format(value, '.3f')
        elif isinstance(value, bool):
            if value:
                text += "yes"
            else:
                text += "no"
        elif not isinstance(value, str):
            text += str(value)
        else:
            text += value
        # Add an end line character
        text += "\n"
    # Return the formatted text
    return text

def basic_info(experiment):
    """Create formatted text describing the basic Experiment information"""
    # Validate the basic information
    VALIDATE.basic_info(experiment)
    # Store the formatted text here
    text = 'Basic Experiment Information\n'
    # Loop through the attributes
    for name in ['User', 'Type', 'Name', 'File Format', 'Force Field','Folder']:
        text += standard_format(name, experiment[name])
    # Return this formatted text
    return text

def CHARMM_info(experiment):
    """Create formatted text describing how to use CHARMM."""
    # Store the formatted text here
    text = ''
    # Only do this if the Experiment is using the CHARMM force field
    if experiment["Force Field"] == "CHARMM":
        # Validate the force field information
        error = VALIDATE.CHARMM_info(experiment)
        if error != '':
            raise FUNCTIONS.IPRO_IOError(text)
        # Give the text an appropriate header
        text = "\nHow to use CHARMM\n"
        # Since the information is all correct, loop through the attributes
        for name in ['CHARMM Topology Files', 'CHARMM Parameter Files', \
                     'CHARMM Energy Terms', 'CHARMM Iterations', 'Lj Phi']:
            # Include each piece of information
            text += standard_format(name, experiment[name])
    return text

def docking_info(experiment):
    """Create text describing how to run docking"""
    # Validate the docking information
    error = VALIDATE.docking_info(experiment)
    # If there's a problem, raise an error
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # STore the text here
    text = "\nHow to run Docking\n"
    # Go through the attributes
    for attribute in ['Docking Frequency', 'Docking Iterations', \
    'Docking SD Movements', 'Docking SD Rotations', 'Docking Start Temp', \
    'Docking End Temp']:
        text += standard_format(attribute, experiment[attribute])
    return text

def IPRO_info(experiment):
    """Provide a description of how IPRO should be run"""
    # Validate the IPRO information
    error = VALIDATE.IPRO_info(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the formatted text here
    text = "\nHow to run IPRO\n"
    # Go through the attributes
    for name in ['IPRO Iterations', "IPRO Annealing Temperature", \
                 "Annealing Sharing", "Energy Calculation"]:
        text += standard_format(name, experiment[name])
    return text

def refinement_info(experiment):
    """Create text describing how to run structure refinements"""
    # Validate the refinement information
    error = VALIDATE.refinement_info(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the formatted text here
    text = "\nHow to do Structure Refinements\n"
    # Go through the attributes
    for name in ['Do Refinement', 'Refinement Iterations', \
                 'Ensemble Number']:
        # If this is a mutator experiment, skip the do refinement attribute
        if "Type" in experiment and experiment["Type"] == "Mutator" and name \
        == "Do Refinement":
            continue
        text += standard_format(name, experiment[name])
    return text

def rotamer_info(experiment):
    """Describe how to use rotamers"""
    # Validate the rotamer information
    error = VALIDATE.rotamer_info(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store formatted text here
    text = "\nHow to use rotamers\n"
    # Go through the attributes
    for name in ['Rotamer Library', 'Rotamer Window', 'Max Rotamer Number', \
                 'Packing Cutoff', 'Packing Method', 'Packing Selection']:
        text += standard_format(name, experiment[name])
    return text

def solvation_info(experiment):
    """Provide information about how to use implicit solvation"""
    # Validate the information
    error = VALIDATE.solvation_info(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the formatted text here
    text = "\nHow to use Implicit Solvation\n"
    # Go through the variables one at a time...
    text += standard_format("Use Solvation", experiment["Use Solvation"])
    # If solvation is being used, say the type
    if experiment["Use Solvation"]:
        text += standard_format("Solvation Type", \
                                    experiment["Solvation Type"])
        # If LK solvation is being used, list those file
        if experiment["Solvation Type"] == "Lazaridis-Karplus":
            text += standard_format("LK Solvation Files", \
                                        experiment["LK Solvation Files"])
    return text

def Dimers(experiment):
    """If there are Dimers in an experiment, list them"""
    # Validate the Dimer information
    error = VALIDATE.Dimers(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = ''
    # Go through the dimers
    if "Dimers" in experiment:
        for pair in experiment["Dimers"]:
            RHS = "Molecules " + pair[0] + " and " + pair[1]
            text += standard_format("Dimers", RHS)
    # If there is Dimer information, give a heading
    if text != '':
        text = "\nMolecules that are Dimers\n" + text
    return text

def Molecules(experiment):
    """Describe the Molecules being used in an experiment"""
    # Validate the Molecules
    error = VALIDATE.Molecules(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = "\nThe Molecules used in the experiment\n"
    # Go through the Molecules
    for data in experiment["Molecules"]:
        # Create the LHS of the text
        LHS = "Molecule"
        # And the RHS
        RHS = "Molecule "
        # If the Molecule's name is a blank, use '_'
        if data[1] == ' ':
            RHS += "_"
        else:
            RHS += data[1]
        RHS += " from file " + data[0] + " is "
        if data[2].design:
            RHS += "Design Molecule "
        else:
            RHS += "Target Molecule "
        RHS += data[2].name
        # Add this to the text
        text += standard_format(LHS, RHS)
    # Add text listing Dimers, if there are any
    text += Dimers(experiment)
    return text

def DesignPositions(experiment):
    """List what Residues are permitted to mutate"""
    # Validate the Design Positions
    error = VALIDATE.DesignPositions(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = ''
    # If this isn't an OptMAVEn or OptCDR experiment, list the Design Positions
    if "Type" in experiment and experiment["Type"] in ['OptMAVEn', 'OptCDR']:
        pass
    else:
        text += "\nResidues that are Permitted to Mutate\n"
        # Go through the Molecules
        for mn in experiment["Design Positions"]:
            # Go through the Residues
            for rn in experiment["Design Positions"][mn]:
                # Create the RHS
                RHS = "Residue " + rn + " in Molecule " + mn
                text += standard_format("Design Position", RHS)
    return text

def EpitopePositions(experiment):
    """List what Residues belongs to the antigen epitope"""
    # Validate the Design Positions
    error = VALIDATE.EpitopePositions(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = ''
    if "Type" in experiment and experiment["Type"] in ['OptMAVEn', 'OptCDR']:
        text += "\nAntigen epitope residues \n"
        # Go through the Molecules
        for mn in experiment["Epitope Positions"]:
            # Go through the Residues
            for rn in experiment["Epitope Positions"][mn]:
                # Create the RHS
                RHS = "Residue " + rn + " in Molecule " + mn
                text += standard_format("Epitope Position", RHS)
    return text

#mfa#
def OptmavenGrid(experiment):
    """ Print the grid levels for the OptMAVEn grid search. """
    return "\n".join(["{:<30}{}".format("Optmaven Grid:", " ".join(map(str, dimension))) for
            dimension in experiment["Optmaven Grid"]])

def DesignGroups(experiment):
    """Describe what Molecules are simultaneously present in the system"""
    # Validate the Design Groups
    error = VALIDATE.DesignGroups(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = "\nThe Design Groups\n"
    # Go through the groups
    for group in experiment["Design Groups"]:
        RHS = group[0].capitalize() + " binding to Molecules:"
        for i in range(1, len(group)):
            RHS += " " + group[i]
        text += standard_format("Design Group", RHS)
    return text

def OptzymeGroups(experiment):
    """Describe the objectives and Molecule present in the system"""
    # Validate the Optzyme Groups
    error = VALIDATE.OptzymeGroups(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = "\nThe Optzyme Groups\n"
    # Go through the groups
    for group in experiment["Optzyme Groups"]:
        RHS = group[1].capitalize() + " " + group[0].capitalize()
        if group[0].capitalize() == "Kcat":
            RHS += " with a weight of " + str(group[2])
        for i in range(len(group[3])):
            RHS += " using "
            if group[3][i][0] == "improve":
                RHS += "TSA Molecules:"
            else:
                RHS += "ground Molecules:"
            for mol in group[3][i][1:]:
                RHS += " " + mol
            if len(group[3]) == 2 and i == 0:
                RHS += " and"
        text += standard_format("Optzyme Group", RHS)
    return text

def PermittedKinds(experiment):
    """Specific limits on how Residues may mutate"""
    # Validate the permitted kinds
    error = VALIDATE.PermittedKinds(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the formatted text here
    text = ''
    if "Permitted Kinds" in experiment and experiment["Permitted Kinds"] != {}:
        # Create a header for the information
        text += "\nPermitted Kinds of Amino Acid Mutations\n"
        # Get the Molecule names
        mns = experiment["Permitted Kinds"].keys()
        mns.sort()
        # Go through the Molecules
        for mn in mns:
            # Get the Residue names
            rns = experiment["Permitted Kinds"][mn].keys()
            # Sort the Residue names manually
            # Ironically, the first step is to sort them automatically
            rns.sort()
            # Now sort them manually
            sorted = False
            while not sorted:
                sorted = True
                for i in range(len(rns) - 1):
                    try:
                        n1 = int(rns[i])
                    except ValueError:
                        n1 = int(rns[i][:-1])
                    try:
                        n2 = int(rns[i+1])
                    except ValueError:
                        n2 = int(rns[i+1][:-1])
                    # If appropriate, reverse the order of the entries
                    if n2 < n1:
                        temp = rns[i]
                        rns[i] = rns[i+1]
                        rns[i+1] = temp
                        sorted = False
            # Output the information about each Residue
            for rn in rns:
                RHS = 'Residue ' + rn + " in Molecule " + mn + ":"
                for aa in experiment["Permitted Kinds"][mn][rn]:
                    RHS += " " + aa
                text += standard_format("Permitted Kinds", RHS)
    return text

def Mutants(experiment):
    """List the mutations of each mutant"""
    # Validate the mutants
    error = VALIDATE.Mutants(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the text here
    text = "\nMutants\n"
    # Go through the mutants
    for i, mutant in enumerate(experiment["Mutants"]):
        # Go through the mutations
        for mutation in mutant:
            RHS = "Mutant " + str(i+1) + ": Mutate Residue "
            RHS += mutation[1] + " in Molecule " + mutation[0] + " to "
            RHS += mutation[2]
            text += standard_format("Mutation", RHS)
    return text

def fixedAtoms(experiment):
    """Output text about Atoms that may never move"""
    # Error checking is done in the restraints function. Store the text here
    text = ''
    if "Restraints" in experiment and "Fixed Atoms" in experiment["Restraints"]:
        # Get the data in a more accessible name
        data = experiment["Restraints"]["Fixed Atoms"]
        # Go through the Design Groups, Molecules, and Residues
        for gn in data:
            for mn in data[gn]:
                for rn in data[gn][mn]:
                    RHS = "In Residue " + rn + " of Molecule " + mn + " in "
                    if gn == 'all':
                        RHS += "all Design Groups:"
                    else:
                        RHS += "Design Group " + str(gn) + ":"
                    for an in data[gn][mn][rn]:
                        RHS += " " + an
                    # Include that text
                    text += standard_format("Fixed Atoms", RHS)
        # If there is text
        if text != '':
            text = "\nAtoms that may never move\n" + text
    return text

def position_restraints(experiment):
    """Describe position restraints"""
    # Error checking is done in the restraints function. Store the text here
    text = ''
    # If there are position restraints
    if "Restraints" in experiment and "Position" in experiment["Restraints"]:
        # go through the restraints
        for restraint in experiment["Restraints"]["Position"]:
            # Create the text
            RHS = ""
            if restraint[3] == 'all':
                RHS += "All Atoms in "
            else:
                RHS += "Atom " + restraint[3] + " in "
            if restraint[2] == "all":
                RHS += "all Residues in "
            else:
                RHS += "Residue " + restraint[2] + " in "
            if restraint[1] == 'all':
                RHS += "all Molecules in "
            else:
                RHS += "Molecule " + restraint[1] + " in "
            if restraint[0] == 'all':
                RHS += "all Design Groups"
            else:
                RHS += "Design Group " + str(restraint[0])
            # Have a different method for each force field
            if experiment["Force Field"] == "CHARMM":
                RHS += ", using a force constant of " 
                RHS += format(restraint[4], '.3f') + ", to "
                if len(restraint) == 6:
                    RHS += "Design Group " + str(restraint[5])
                elif restraint[3] == 'all':
                    RHS += "their initial positions"
                else:
                    RHS += "its initial position"
            else:
                text = "The I/O OUTPUT position restraints function does not "
                text += 'support the ' + str(experiment["Force Field"])
                text += " force field."
                raise FUNCTIONS.IPRO_IOError(text)
            # Store the information in the text
            text += standard_format("Position Restraint", RHS)
        # If there are position restraints, modify the text
        if text != '':
            text = "\nRestraints on Atom Positions\n" + text
    return text

def distance_restraints(experiment):
    """Create text describing distance restraints"""
    # Error checking is done in the restraints function. Store the distance
    # restraints here
    text = ''
    # If there are restraints
    if "Restraints" in experiment and "Distance" in experiment["Restraints"]:
        # Go through the restraints
        for restraint in experiment["Restraints"]["Distance"]:
            # Start creating the RHS
            RHS = "Between "
            # List the two atoms
            for i in range(1, 3):
                RHS += "Atom " + restraint[i][2] + " in Residue " + \
                       restraint[i][1] + " in Molecule " + restraint[i][0]
                if i == 1:
                    RHS += " and "
            # Say the Design Group
            if restraint[0] == 'all':
                RHS += " in all Design Groups"
            else:
                RHS += " in Design Group " + str(restraint[0])
            # Have a different method for the remaining values based on force
            # field
            if experiment["Force Field"] == "CHARMM":
                names = ['KMIN', 'RMIN', 'KMAX', 'RMAX', 'FMAX']
                for j in range(len(names)):
                    RHS += " " + names[j] + ": " + format(restraint[3+j], '.3f')
            else:
                text = "The I/O OUTPUT distance restraints function does not "
                text += "support the " +experiment["Force Field"]+" force field"
                raise FUNCTIONS.IPRO_IOError(text)
            # Include the text
            text += standard_format("Distance Restraint", RHS)
        # If text was created
        if text != '':
            text = "\nRestraints on Atom-Atom distances\n" + text
    return text

def dihedral_restraints(experiment):
    """Create text describing dihedral angle restraints"""
    # Error checking is done in the restraints function. Store the text here
    text = ''
    # If there are restraints
    if "Restraints" in experiment and "Dihedral" in experiment["Restraints"]:
        # Go through them
        for restraint in experiment["Restraints"]["Dihedral"]:
            # Start creating this restraint
            RHS = "Between "
            for i in range(1, 5):
                if i == 4:
                    RHS += "and "
                RHS += "Atom " + restraint[i][2] + " in Residue "
                RHS += restraint[i][1] + " in Molecule " + restraint[i][0]
                if i != 4:
                    RHS += ", "
            # Say the Design Group
            if restraint[0] == 'all':
                RHS += " in all Design Groups"
            else:
                RHS += " in Design Group " + str(restraint[0])
            # Have a different method for the remaining terms based on force
            # field
            if experiment["Force Field"] == "CHARMM":
                RHS += " using a force constant of "+format(restraint[5],'.3f')
                RHS += " to a minimum angle of " + format(restraint[6], '.3f')
            else:
                text = "The I/O OUTPUT dihedral restraints function does not "
                text += "support the "+experiment["Force Field"]+" force field."
                raise FUNCTIONS.IPRO_IOError(text)
            # Add the restraint to the text
            text += standard_format("Dihedral Restraint", RHS)
        # If there are restraints, give htem a header
        if text != '':
            text = "\nRestraints on Dihedral Angles\n" + text
    return text

def Restraints(experiment):
    """Create text describing the 4 types of structure restraints"""
    # Validate the restraints that exist
    error = VALIDATE.Restraints(experiment)
    if error != '':
        raise FUNCTIONS.IPRO_IOError(error)
    # Store the formatted text here
    text = ''
    # Add the text for each type of restraint
    text += fixedAtoms(experiment)
    text += position_restraints(experiment)
    text += distance_restraints(experiment)
    text += dihedral_restraints(experiment)
    return text
