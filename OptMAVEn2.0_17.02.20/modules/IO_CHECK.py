#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Input / Output Value Checking"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for checking values input by a user for use in an
IPRO Suite Experiment."""

# Include standard python modules
import os
import sys
# Include the FUNCTIONS module for I / O
import IO_FUNCTIONS as FUNCTIONS
# Include all contents of the STANDARDS module
from STANDARDS import *
# And include access to other needed modules
import MOLECULES

def has_whitespace(value, slash = False):
    """Check a string for whitespace"""
    # If it is not a string or has whitespace, that is a problem
    if not isinstance(value, str) or ' ' in value or '\t' in value or '\n' in \
    value:
        return True
    # Similarly, if the value should not have a slash in it (indicating a
    # different folder), that is a problem
    elif slash and '/' in value:
        return True
    # Otherwise there is no problem
    else:
        return False

def for_file(fileName, path = None, OPEN = False):
    """Check for a file in standard IPRO Suite locations"""
    # Validate the file's name
    if has_whitespace(fileName, True):
        text = "The name of files used in the IPRO Suite of Programs must not "
        text += "contain whitespace or '/', which would indicate that they are "
        text += "in a different folder."
        raise FUNCTIONS.IPRO_IOError(text)
    # If there is a problem with the path, reset to the current
    if has_whitespace(path):
	#mfa# path = os.getcwd() + "/"
	path = InstallFolder
    # There are three locations to check for the file
    places = [path, os.path.join(path, "input_files"), os.path.join(path,
            "structures")] #mfa#
    # Determine whether or not the file can be found
    found = False
    # Check each location
    for place in places:
        try:
            if fileName in os.listdir(place): #mfa#
                f = os.path.join(place, fileName)
                # If the file should be opened and returned, do so.
                if OPEN:
                    return open(f, "r"), f #mfa#
                # Otherwise, return the full file path.
                return f #mfa#
        # If there's an error, skip this folder
        except (OSError, IOError):
            pass
    # If the file wasn't found, raise an error
    if not found:
        text = "The " + fileName + " file could not be found in any expected "
        text += "folder."
        raise FUNCTIONS.IPRO_IOError(text)

# Have functions to check the attributes of various terms, based on what they
# are used for in the IPRO Suite of Programs
def basic_info(attribute, value):
    """Check the values of attributes that every Experiment needs"""
    # Who is running the experiment
    if attribute == "User":
        if not isinstance(value, str) or len(value) == 0:
            text = "The name of an IPRO Suite Experiment User must be a string "
            text += "of at least one character."
            raise FUNCTIONS.IPRO_IOError(text)
    # What type of Experiment it is
    elif attribute == "Type":
        # supportedPrograms is a list from STANDARDS
        if value not in supportedPrograms:
            text = "The Type of an IPRO Suite Experiment must be in the "
            text += "supported programs list of the STANDARDS module."
            raise FUNCTIONS.IPRO_IOError(text)
    # The name of an Experiment
    elif attribute == "Name":
        if has_whitespace(value, True) or len(value) == 0:
            text = "The Name of an IPRO Suite Experiment must be a string that "
            text += "does not contain whitespace or '/'"
            raise FUNCTIONS.IPRO_IOError(text)
        # Find out if this name is unique or not
        unique = True
        if "experiments" in os.listdir("./"):
            if value in os.listdir("experiments/"):
                unique = False
        elif value in os.listdir("./"):
            unique = True
        if not unique:
            text = value + " is an existing experiment name, so it may not be "
            text += "reused."
            raise FUNCTIONS.IPRO_IOError(text)
    # The type of file format
    elif attribute == "File Format":
        # Supported formats is from STANDARDS
        if value not in supportedFormats:
            text = str(value) + " is not a file format supported by the "
            text += "IPRO Suite of Programs."
            raise FUNCTIONS.IPRO_IOError(text)
    # The force field
    elif attribute == "Force Field":
        if value not in supportedFields:
            text = str(value) + " is not a force field supported by the IPRO "
            text += "Suite of Programs."
            raise FUNCTIONS.IPRO_IOError(text)
    # The folder of the Experiment
    elif attribute == "Folder":
        if has_whitespace(value):
            text="An IPRO Suite Experiment's Folder may not contain whitespace."
            raise FUNCTIONS.IPRO_IOError(text)
    # Otherwise, raise a generic error
    else:
        text = attribute + " is not a basic attribute of an IPRO Suite "
        text += "Experiment."
        raise FUNCTIONS.IPRO_IOError(text)

def CHARMM_info(attribute, value):
    """Check attributes associated with the use of CHARMM."""
    # If it is a list of CHARMM topology or parameter files, or the energy terms
    if attribute in ['CHARMM Topology Files', 'CHARMM Parameter Files', \
                     'CHARMM Energy Terms']:
        # If it isn't a list or is empty
        if not isinstance(value, list) or len(value) == 0:
            text = "The " + attribute + " must be a list of lowercase strings."
            raise FUNCTIONS.IPRO_IOError(text)
        # Check each File
        errors = ''
        for fileName in value:
            # Check that the expected file exists (otherwise an IPRO IO Error
            # will be raised by the for_file function)
            if attribute in ['CHARMM Topology Files', 'CHARMM Parameter Files']:
                try:
                    for_file(fileName)
                    if not fileName == fileName.lower():
                        errors += "\n" + fileName+" is not a lower case string."
                except FUNCTIONS.IPRO_IOError as error:
                    errors += str(error)
        if errors != '':
            raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Check the iterations
    elif attribute == "CHARMM Iterations":
        if not isinstance(value, int) or value < 50:
            text = "The number of " + attribute + " must be an integer greater "
            text += "than or equal to 50."
            raise FUNCTIONS.IPRO_IOError(text)
    # Check the softening term
    elif attribute == "Lj Phi":
        if not isinstance(value, float) or value > 1.0:
            text = "The value of " + attribute + " must be a floating point "
            text += "number between 0 and 1."
            raise FUNCTIONS.IPRO_IOError(text)
    # Otherwise, the attribute isn't for CHARMM
    else:
        text = attribute + " is not a CHARMM-related IPRO Suite Experiment "
        text += "attribute."
        raise FUNCTIONS.IPRO_IOError(text)

def docking_info(attribute, value, experiment = {}):
    """Check information associated with how to run docking"""
    # Positive integers
    if attribute in ['Docking Frequency', 'Docking Iterations']:
        if not isinstance(value, int) or value <= 0:
            text = "The " + attribute + " must be a positive integer."
            raise FUNCTIONS.IPRO_IOError(text)
    # Positive numbers
    elif attribute in ['Docking SD Movements', 'Docking SD Rotations', \
                       'Docking Start Temp', 'Docking End Temp']:
        if not isinstance(value, (float, int)) or value <= 0:
            text = "The " + attribute + " must be a positive number."
            raise FUNCTIONS.IPRO_IOError(text)
        # Docking temperatures must be compared with each other
        if attribute == "Docking Start Temp" and "Docking End Temp" in \
        experiment and value < experiment["Docking End Temp"]:
            text = "The " + attribute + " must be greater than the temperature "
            text += "at the end of docking."
            raise FUNCTIONS.IPRO_IOError(text)
        elif attribute == "Docking End Temp" and "Docking Start Temp" in \
        experiment and value > experiment["Docking Start Temp"]:
            text = "The " + attribute + " must be less than the temperature at "
            text += "the beginning of docking."
            raise FUNCTIONS.IPRO_IOError(text)
    # There aren't any other docking attributes
    else:
        text = attribute + " is not a docking-related IPRO Suite Experiment "
        text += "attribute."
        raise FUNCTIONS.IPRO_IOError(text)

def rotamer_info(attribute, value):
    """Check the values of attributes associated with the use of Rotamers"""
    # If this is the location of the rotamer library
    if attribute == "Rotamer Library":
        # Make sure it matches expected formats
        try:
            if has_whitespace(value) or "independ" not in os.listdir(value) or \
            "gly1.pdb" not in os.listdir(value + "independ/"):
                raise OSError
        except OSError:
            text = "The specified folder does not meet the expected rotamer "
            text += "library formatting requirements."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it must be a positive number
    elif attribute in ['Rotamer Window', 'Packing Cutoff']:
        if not isinstance(value, (int, float)) or value <= 0:
            text = "The " + attribute + " value must be a positive number."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it must be a positive integer
    elif attribute == "Max Rotamer Number":
        if not isinstance(value, int) or value < 50:
            text = "The " + attribute + " value must be an integer greater than"
            text += " or equal to 50."
            raise FUNCTIONS.IPRO_IOError(text)
    # Other values
    elif attribute == "Packing Method":
        if value not in ['Sequence', 'Distance']:
            text = "The Packing Method must be either 'Sequence' or 'Distance' "
            text += "based"
            raise FUNCTIONS.IPRO_IOError(text)
    elif attribute == "Packing Selection":
        if value not in ["RESTRAINED", "FREE"]:
            text = "The Packing Selection must be either 'RESTRAINED' Residues "
            text += "or 'FREE' ones, too."
            raise FUNCTIONS.IPRO_IOError(text)
    # There are no other attributes
    else:
        text = attribute + " is not a rotamer-related IPRO Suite Experiment "
        text += "attribute."
        raise FUNCTIONS.IPRO_IOError(text)

def solvation_info(attribute, value):
    """Check the values of solvation related attributes"""
    # If it is whether or not to use solvation in general or in a particular
    # procedure
    if attribute in ['Use Solvation', 'Relaxation Solvation', \
                     'Perturbation Solvation', 'Energy Solvation']:
        if not isinstance(value, bool):
            text = "The " + attribute + " must be either True or False."
            raise FUNCTIONS.IPRO_IOError(text)
    # What type of solvation to use during non-bonded energy calculations
    elif attribute == "Solvation Type":
        if value != None and value not in supportedSolvations:
            text = "The " + attribute + " value must be either None or an "
            text += "entry in the supportedSolvations list in the STANDARDS "
            text += "module."
            raise FUNCTIONS.IPRO_IOError(text)
    # The list of LK solvation files
    elif attribute == "LK Solvation Files":
        if not isinstance(value, list) or len(value) == 0:
            text = "The " + attribute + " must be a list of strings."
            raise FUNCTIONS.IPRO_IOError(text)
        # If we got past that, check each entry
        errors = ''
        for fileName in value:
            try:
                for_file(fileName)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
        if errors != '':
            raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Everything else should raise an error
    else:
        text = attribute + " is not an implicit solvation related attribute "
        text += "of an IPRO Suite Experiment."
        raise FUNCTIONS.IPRO_IOError(text)

def list_check(attribute, value):
    """Check attributes that should be gathered by the get_list function"""
    if attribute in ['CHARMM Topology Files', 'CHARMM Parameter Files', \
                     'CHARMM Energy Terms']:
        CHARMM_info(attribute, value)
    elif attribute == "LK Solvation Files":
        solvation_info(attribute, value)
    else:
        text = "The list check function does not support checking the "
        text += attribute + " attribute."
        raise FUNCTIONS.IPRO_IOError(text)

def IPRO_info(attribute, value):
    """Check the values of IPRO associated attributes"""
    # If it should be a positive integer
    if attribute == "IPRO Iterations":
        if not isinstance(value, int) or value <= 0:
            text = "The " + attribute + " value must be a positive integer."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it should be a positive number
    elif attribute == "IPRO Annealing Temperature":
        if not isinstance(value, (float, int)) or value <= 0:
            text = "The " + attribute + " value must be a positive number."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it should be a boolean value
    elif attribute == "Annealing Sharing":
        if value not in [True, False]:
            text = "The " + attribute + " value must be either True or False."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it has a specific string value
    elif attribute == "Energy Calculation":
        if value not in ['Binding', 'Interaction']:
            text = "The " + attribute + " must be either 'Binding' or "
            text += "'Interaction'."
            raise FUNCTIONS.IPRO_IOError(text)
    # Otherwise just raise an error
    else:
        text = attribute + " is not an attribute associated with how to run "
        text += "IPRO during IPRO Suite experiments."
        raise FUNCTIONS.IPRO_IOError(text)

def Molecule(molecule):
    """Check that a list of information about a Molecule is correct."""
    # Make sure the information is in a list with an acceptable number of
    # entries
    if not isinstance(molecule, list) or not 3 <= len(molecule) <= 4:
        text = "The information about a Molecule to use in an IPRO Suite "
        text += "Experiment must be stored in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check the entries in the list
    elif not isinstance(molecule[0], str):
        text = "The name of a Molecule's file must be a string."
        raise FUNCTIONS.IPRO_IOError(text)
    elif not isinstance(molecule[1], str):
        text = "The name of a Molecule must be a string."
        raise FUNCTIONS.IPRO_IOError(text)
    elif not isinstance(molecule[2], MOLECULES.Molecule):
        text = "The third entry in a list of information about a Molecule must "
        text += "be a Molecule class object."
        raise FUNCTIONS.IPRO_IOError(text)
    # If there is a fourth entry, make sure the label is valid
    elif len(molecule) == 4 and molecule[3] not in ['Design Molecule', \
    'Target Molecule', 'Enzyme', 'Protein', 'Antibody', 'Substrate', 'Ligand', \
    'Cofactor', 'Antigen']:
        text = molecule[3] + " is not an acceptable label for what a Molecule "
        text += "is."
        raise FUNCTIONS.IPRO_IOError(text)

def Molecules(molecules):
    """Check that the Molecules information is compatible."""
    # Make sure it is a non-empty list
    if not isinstance(molecules, list) or len(molecules) == 0:
        text = "The information about the Molecules used in an IPRO Suite "
        text += "Experiment must be stored in a list and at least one Molecule "
        text += "must be used."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check each entry
    errors = ''
    for molecule in molecules:
        try:
            Molecule(molecule)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # If there are any problems with the Molecules, say that
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Confirm that no Molecule is used twice and that each Molecule name is
    # unique
    for i in range(len(molecules) - 1):
        for j in range(i+1, len(molecules)):
            if molecules[i][0] == molecules[j][0] and molecules[i][1] == \
            molecules[j][1]:
                errors += "\nMolecule " + molecules[i][1] + " from file "
                errors += molecules[i][0] + " is used more than once."
            elif molecules[i][2].name == molecules[j][2].name:
                errors += "\n" + molecules[i][2].name + " is used as a Molecule"
                errors += " name more than once."
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def Dimer(info, experiment):
    """Check that information about a dimer is valid."""
    # make sure the it is a list
    if not isinstance(info, list) or len(info) != 2:
        text = "The information about a pair of Molecules that are a Dimer must"
        text += " be stored in a list of TWO items."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check that each Molecule is unique
    elif info[0] == info[1]:
        text = "The same Molecule name cannot be used for both Molecules in a "
        text += "dimer."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check that the each Molecule is part of the Experiment
    else:
        for mn in info:
            # Find the Molecule
            have = False
            for data in experiment["Molecules"]:
                if data[2].name == mn:
                    have = True
                    # If the Molecule isn't a Design Molecule, raise an error
                    if not data[2].design:
                        text = "Molecule " + mn + " is not a Design Molecule, "
                        text+= "so it cannot be listed in the Dimer information"
                        raise FUNCTIONS.IPRO_IOError(text)
                    break
            # If it wasn't found, raise an error
            if not have:
                text = "There is no " + mn + " Molecule, so it may not be "
                text += "listed as a Dimer."
                raise FUNCTIONS.IPRO_IOError(text)

def Dimers(dimers, experiment):
    """Check that the dimers information is valid."""
    # Make sure the information is stored in a list
    if not isinstance(dimers, list):
        text = "The information about dimers must be stored in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any generated errors here
    errors = ''
    # Check each dimer entry
    for dimer in dimers:
        try:
            Dimer(dimer, experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # If there is a problem raise an error
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Check that each dimer is unique
    for i in range(len(dimers) - 1):
        for j in range(i+1, len(dimers)):
            if dimers[i][0] in dimers[j] and dimers[i][1] in dimers[j]:
                errors += "\nDimer listings " + str(i+1) + " and " + str(j+1)
                errors += " are identical."
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def DesignPosition(mn, rn, experiment):
    """Determine if a Residue may be a Design Position or not."""
    # Find the molecule
    molecule = None
    for data in experiment["Molecules"]:
        if data[2].name == mn:
            molecule = data[2]
            break
    # If the Molecule couldn't be found
    if molecule == None:
        text = "There is no " + mn + " Molecule, so it may not contain a Design"
        text += " Position."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure it is a Design Molecule
    elif not molecule.design:
        text = "Molecule " + mn + " is not a Design Molecule, so it may not "
        text += "contain Design Positions."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the Residue isn't in the Molecule
    elif rn not in molecule:
        text = "Molecule " + mn + " does not contain Residue " + rn + ", so it "
        text += "may not be a Design Position."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the Residue is an amino acid
    elif molecule[rn].kind not in aminoAcids[molecule.fileFormat]:
        text = "Residue " + rn + " in Molecule " + mn + " is not an amino acid,"
        text += " so it may not be a Design Position."
        raise FUNCTIONS.IPRO_IOError(text)

def DesignPositions(positions, experiment):
    """Check the validity of all of the Design Positions."""
    # Make sure they are stored in a Dictionary
    if not isinstance(positions, dict):
        text = "The Design Positions must be stored in a dictionary."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any generated errors here
    errors = ''
    count = 0
    for mn in positions:
        # If it isn't a list, raise an error
        if not isinstance(positions[mn], list):
            errors += "\nThe Design Positions for a Molecule must be stored in "
            errors += "a list."
            continue
        # Check each Residue
        for rn in positions[mn]:
            try:
                DesignPosition(mn, rn, experiment)
                count += 1
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # If there are no Design Positions
    if count == 0 and errors == '':
        errors += "\nThere are no provided Design Positions"
    # If there are errors, raise an error
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def EpitopePosition(mn, rn, experiment):
    """Determine if a Residue may be a Epitope Position or not."""
    # Find the molecule
    molecule = None
    for data in experiment["Molecules"]:
        if data[2].name == mn:
            molecule = data[2]
            break
    # If the Molecule couldn't be found
    if molecule == None:
        text = "There is no " + mn + " Molecule, so it may not contain a Epitope"
        text += " Position."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure it is a Design Molecule
    elif molecule.design:
        text = "Molecule " + mn + " is not a Target Molecule, so it may not "
        text += "contain Epitope Positions."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the Residue isn't in the Molecule
    elif rn not in molecule:
        text = "Molecule " + mn + " does not contain Residue " + rn + ", so it "
        text += "may not be a Epitope Position."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the Residue is an amino acid
    elif molecule[rn].kind not in aminoAcids[molecule.fileFormat]:
        text = "Residue " + rn + " in Molecule " + mn + " is not an amino acid,"
        text += " so it may not be a Epitope Position."
        raise FUNCTIONS.IPRO_IOError(text)

def EpitopePositions(positions, experiment):
    """Check the validity of all of the Epitope Positions."""
    # Make sure they are stored in a Dictionary
    if not isinstance(positions, dict):
        text = "The Epitope Positions must be stored in a dictionary."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any generated errors here
    errors = ''
    count = 0
    for mn in positions:
        # If it isn't a list, raise an error
        if not isinstance(positions[mn], list):
            errors += "\nThe Epitope Positions for a Molecule must be stored in "
            errors += "a list."
            continue
        # Check each Residue
        for rn in positions[mn]:
            try:
                EpitopePosition(mn, rn, experiment)
                count += 1
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # If there are no Design Positions
    if count == 0 and errors == '':
        errors += "\nThere are no provided Epitope Positions"
    # If there are errors, raise an error
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def DesignGroup(group, experiment):
    """Check that the information describing a Design Group is valid."""
    # Make sure it is a list
    if not isinstance(group, list) or len(group) == 0:
        text = "The information about a Design Group must be stored in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the first entry is a binding description
    elif group[0] not in ['improve', 'maintain', 'reduce', 'eliminate']:
        text = str(group[0]) + " is not a recognized binding specification for "
        text += "a Design Group."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check each of the following entries
    for i in range(1, len(group)):
        # Find the Molecule
        have = False
        for data in experiment["Molecules"]:
            if data[2].name == group[i]:
                have = True
                # If it is a Design Molecule, that is bad (here)
                if data[2].design:
                   text="Design Molecules should not be listed in Design Groups"
                   raise FUNCTIONS.IPRO_IOError(text)
                break
        # If the Molecule doesn't exist, that's a problem
        if not have:
            text = "There is no " + mn + " Molecule, so it may not be in a "
            text += "Design Group."
            raise FUNCTIONS.IPRO_IOError(text)

def DesignGroups(groups, experiment):
    """Check all Design Groups."""
    # Make sure it is a list with at least one entry
    if not isinstance(groups, list) or len(groups) == 0:
        text = "The information about the Design Groups must be stored in a "
        text += "list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any generated errors here
    errors = ''
    # Check each Design Group
    for group in groups:
        try:
            DesignGroup(group, experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # If there are errors, be done
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Make sure the Design Groups are unique
    for i in range(len(groups) - 1):
        for j in range(i+1, len(groups)):
            # If they're different lengths, skip the comparison (unique by
            # definition)
            if len(groups[i]) != len(groups[j]):
                continue
            # Check each molecule
            same = True
            for mn in groups[i][1:]:
                if mn not in groups[j]:
                    same = False
                    break
            # If they are the same, say that error
            if same:
                errors += "\nDesign Groups " + str(i+1) + ' and ' + str(j+1)
                errors += " contain all the same Molecules."
    # Also make sure that each Target Molecule is used in at least one Design
    # Group
    targets = {}
    for data in experiment["Molecules"]:
        if not data[2].design:
            targets[data[2].name] = False
    # Go through the groups
    for group in groups:
        # Go through each Molecule
        for i in range(1, len(group)):
            # Indicate that the Molecule is used
            targets[group[i]] = True
    # Include any Molecule that is not used in the errors
    for mn in targets:
        if not targets[mn]:
            errors += "\nTarget Molecule " + mn + " is not used in any Design "
            errors += "Group."
    # if there are errors
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(text)

def OptzymeGroup(group, experiment):
    """Check that the information describing an Optzyme Group is valid."""
    # Make sure that it is a list of nonzero length
    if not isinstance(group, list) or len(group) == 0:
        print group
        text = "The information about the Optzyme Groups must be stored in a "
        text += "list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure that the first entry is a catalytic property
    elif group[0] not in ["km", "kcat/km", "kcat"]:
        text = str(group[0]) + " is not a recognized catalytic property for "
        text += "an Optzyme Group."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the second entry is a binding description      
    elif group[1] not in ['improve', 'maintain', 'reduce', 'eliminate']:
        text = str(group[1]) + " is not a recognized binding specification "
        text += "for an Optzyme Group."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the weight is appropriate, if applicable
    elif group[0] == "kcat" and float(group[2]) <= 0.0:
        text = str(group[2]) + " is not a valid value to be used for the "
        text += "Interaction Energy difference weight."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the proper number of Binding Assemblies are included
    elif (group[0] == "kcat" and len(group[3]) != 2) or \
         (group[0] in ["km", "kcat/km"] and len(group[3]) != 1):
        text = "There cannot be " + str(len(group[3])) + " Binding Assemblies " 
        text += "within an Optzyme group modifying " + str(group[0]) + "."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the target molecules are not the same within a kcat Group
    elif (group[0] == "kcat") and (group[3][0][1:] == group[3][1][1:]):
        text = "There cannot be an Optzyme Group with redundant Target "
        text += "Molecules between 2 Binding Assemblies."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check each of the following entries
    for i in range(len(group[3])):
        # Check the binding description for each Binding Assembly
        if group[3][i][0] not in ['improve', 'eliminate']:
            text = str(group[3][i][0]) + " is not a recognized binding "
            text += "specification for a Binding Assembly."
            raise FUNCTIONS.IPRO_Error(text)
        # Check each of the Target Molecules    
        for j in range(1, len(group[3][i])):
            # Find the Molecule
            have = False
            for data in experiment["Molecules"]:
                if data[2].name == group[3][i][j]:
                    have = True
                    # It should be a Target Molecule, not a Design Molecule
                    if data[2].design:
                        text = "Design Molecules should not be listed within "
                        text += "a Binding Assembly."
                        raise FUNCTIONS.IPRO_IOError(text)
                    break
            # If the Molecule does not exit, that is a problem
            if not have:
                text = "There is no " + group[3][i][j] + " Molecule, so it may "
                text += "not be used within a Binding Assembly."
                raise FUNCTIONS.IPRO_IOError(text)

def OptzymeGroups(groups, experiment):
    """Check all of the Optzyme Groups."""
    # Make sure it is a list with at least one entry
    if not isinstance(groups, list) or len(groups) == 0:
        text = "The information about the Optzyme Groups must be stored in a "
        text += "list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any generated errors here
    errors = ''
    # Check each Optzyme Group
    for group in groups:
        try:
            OptzymeGroup(group, experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # If there are errors, stop
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])
    # Make sure the Optzyme Groups are unique
    for i in range(len(groups) - 1):
        for j in range(i+1, len(groups)):
            # If the two objects are identical, raise an error
            if groups[i] == groups[j]:
                errors += "\nOptzyme Groups " + str(i+1) + " and " + str(j+1)
                errors += " are equivalent."
            # Make sure the same set of Target Molecules with different
            # states are not used
            for k in [0, 1]:
                for l in [0, 1]:
                    # Skip non-existent entries
                    if ((k+1) > len(groups[i])) or ((l+1) > len(groups[j])):
                        continue
                    if (groups[i][k][1:] == groups[j][l][1:]) and \
                       (groups[i][k][0] != groups[j][l][0]):
                        errors += "\nOptzyme Groups " + str(i+1) + " and " 
                        errors += str(j+1) + " have the same Target Molecules "
                        errors += "at different states."
    # Also make sure that at each Target Molecule is used in at least one
    # Optzyme Group
    targets = []
    for data in experiment["Molecules"]:
        if not data[2].design:
            # State that the Target Molecule was not used
            used = False
            # Go through the Optzyme Groups 
            for group in groups:
                for assembly in group[3]:
                    if data[2].name in assembly:
                        used = True
                        break
            # If the Target Molecule was unused, include it as a missing
            # Molecule
            if not used:
                targets.append(data[2].name)
    # If any Target Molecules are not used, include them in the errors
    for mn in targets:
        errors += "\nTarget Molecule " + mn + " is not used in any Optzyme "
        errors += "Group."
    # If there are errors
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors)

def residue_permittedKinds(mn, rn, permitted, experiment):
    """Check the permitted kinds list for a Residue."""
    # Find the relevant Molecule
    molecule = None
    for data in experiment["Molecules"]:
        if data[2].name == mn:
            molecule = data[2]
            break
    if molecule == None:
        text = "There is no Molecule " + mn + ", so its Residues may not have "
        text += "specific kinds of mutations allowed."
        raise FUNCTIONS.IPRO_IOError(text)
    # Confirm that it is a Design Molecule
    if not molecule.design:
        text = "Molecule " + mn + " is not a Design Molecule, so its Residues "
        text += "may not have specified kinds of mutations allowed."
        raise FUNCTIONS.IPRO_IOError(text)
    # Confirm the Residue is in the Molecule
    if rn not in molecule:
        text = "Molecule " + mn + " does not contain Residue " + rn
        raise FUNCTIONS.IPRO_IOError(text)
    # Confirm the Residue is a Design Position
    if mn not in experiment["Design Positions"] or rn not in \
    experiment["Design Positions"][mn]:
        text = "Residue " + rn + " in Molecule " + mn + " is not a Design "
        text += "Position, so it may not have specified kinds of mutations."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check the permissions
    if not isinstance(permitted, list):
        text = "The permitted kinds of amino acids for a Residue must be stored"
        text += " in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check each entry
    for kind in permitted:
        if kind not in aminoAcids[molecule.fileFormat]:
            text = "Residues are only permitted to mutate to known amino acids"
            raise FUNCTIONS.IPRO_IOError(text)

def PermittedKinds(permitted, experiment):
    """Check all permitted kinds of amino acids"""
    # Make sure they are stored in a dictionary
    if not isinstance(permitted, dict):
        text = "The permitted kinds of amino acids for mutations must be "
        text += "stored in a dictionary."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any further errors generated here
    errors = ''
    # Go through the Molecules
    for mn in permitted:
        # Make sure htis is a dictionary
        if not isinstance(permitted[mn], dict):
            errors += "\nThe permitted kinds of amino acids for Residues in a "
            errors += "Molecule must be stored in a dictionary."
            continue
        # Go through the Residues
        for rn in permitted[mn]:
            try:
                residue_permittedKinds(mn, rn, permitted[mn][rn], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # If there were any errors, share them all
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors)

def Mutation(mn, rn, aa, experiment):
    """Make sure a mutation is acceptable."""
    # Find the Molecule
    molecule = None
    for data in experiment["Molecules"]:
        if data[2].name == mn:
            molecule = data[2]
            break
    # If the Molecule doesn't exist
    if molecule == None:
        text = "There is no " + mn + " Molecule, so it may not be mutated."
        raise FUNCTIONS.IPRO_IOError(text)
    # If it isn't a Design Molecule
    if not molecule.design:
        text = "Molecule " + mn + " is not a Design Molecule, so it may not be "
        text += "mutated."
        raise FUNCTIONS.IPRO_IOError(text)
    # If there is no rn Residue
    if rn not in molecule:
        text = "Molecule " + mn + " does not contain a " + rn + " Residue, so "
        text += "it may not be mutated."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the Residue isn't an amino acid
    if molecule[rn].kind not in aminoAcids[molecule.fileFormat]:
        text = "Residue " + rn + " is not an amino acid, so it may not be "
        text += "mutated."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the mutation is to something that isn't an amino acid
    if aa not in aminoAcids[molecule.fileFormat]:
        text = aa + " is not an amino acid, so a Residue may not be mutated "
        text += "to it."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the amino acid isn't changing
    if aa == molecule[rn].kind:
        text = "Residue " + rn + " in Molecule " + mn + " is already a " + aa
        text += ", so it cannot be mutated to that."
        raise FUNCTIONS.IPRO_IOError(text)

def Mutant(mutant, experiment):
    """Check that a particular mutant is acceptable."""
    # Make sure it is a list
    if not isinstance(mutant, list) or len(mutant) == 0:
        text = "The information about a mutant must be stored in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any further errors here
    errors = ''
    # Go through the mutations
    for mutation in mutant:
        # Make sure it is a list of 3 items
        if not isinstance(mutation, list) or len(mutation) != 3:
            errors += "\nThe information about a mutation in a mutant must be "
            errors += "stored in a list of three items."
            continue
        # Check the mutation
        try:
            Mutation(mutation[0], mutation[1], mutation[2], experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # If there are no errors, check to make sure each mutation is to a different
    # Residue
    if errors == '':
        for i in range(len(mutant) - 1):
            for j in range(i+1, len(mutant)):
                # If it is the same Residue, that's a problem
                if mutant[i][:2] == mutant[j][:2]:
                    errors += "\nResidue " + mutant[i][1] + " in Molecule "
                    errors += mutant[i][0] + " is mutated more than once in the"
                    errors += " same mutant."
    # If there are errors, tell the user
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def Mutants(mutants, experiment):
    """Check all mutations of an experiment."""
    # Make sure it is a list of items
    if not isinstance(mutants, list) or len(mutants) == 0:
        text = "The information about mutants must be stored in a list."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any further errors here
    errors = ''
    for i, mutant in enumerate(mutants):
        # Check the mutant
        try:
            Mutant(mutant, experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += "\nMutant " + str(i+1) + " has the following problems:"
            errors += str(error)
    # If there are no problems, make sure each mutant is unique
    if errors == '':
        for i in range(len(mutants) - 1):
            for j in range(i+1, len(mutants)):
                # If they aren't the same length, they're different
                if len(mutants[i]) != len(mutants[j]):
                    continue
                # Find out if they're the same or not
                same = True
                for mutation in mutants[i]:
                    if mutation not in mutants[j]:
                        same = False
                        break
                # If they're the same, tell the user
                if same:
                    errors += "\nMutants " + str(i+1) + " and " + str(j+1)
                    errors += " are the same."
    # If there are errors, tell the user
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def refinement_info(attribute, value):
    """Make sure Refinement related attributes are correct."""
    # If it should be a boolean value
    if attribute == "Do Refinement":
        if not isinstance(value, bool):
            text = "The " + attribute + " value must be True or False."
            raise FUNCTIONS.IPRO_IOError(text)
    # If it must be a positive integer
    elif attribute in ['Refinement Iterations', 'Ensemble Number']:
        if not isinstance(value, int) or value <= 0:
            text = "The " + attribute + " must be a positive integer."
            raise FUNCTIONS.IPRO_IOError(text)
    else:
        text = attribute + " is not a refinement-related IPRO Suite Experiment "
        text += "attribute."
        raise FUNCTIONS.IPRO_IOError(text)

def Atom(gn, mn, rn, an, experiment, what):
    """Check whether an Atom may be used in a restraint."""
    # Check the group information
    if gn != 'all' and gn not in range(1, len(experiment["Design Groups"]) + 1):
        text = str(gn) + " is not a valid Design Group specification"
        raise FUNCTIONS.IPRO_IOError(text)
    # Find the Molecule
    molecule = None
    for data in experiment["Molecules"]:
        if data[2].name == mn:
            molecule = data[2]
            break
    if molecule == None:
        text = "There is no " + mn + " Molecule, so its Atoms may not be " +what
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the Molecule is in the relevant Design Group
    if not molecule.design and gn != 'all':
        if mn not in experiment["Design Groups"][gn - 1]:
            text = "Molecule " + mn + " is not a part of Design Group " +str(gn)
            text += ", so its Atoms may not be " + what
            raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the Residue is in the Molecule
    if rn not in molecule:
        text = "Molecule " + mn + " does not contain a " + rn + " Residue, so "
        text += "its Atoms may not be " + what
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure the Atom is in the Residue
    if an not in molecule[rn]:
        text = "Residue " + rn + " in Molecule " + mn + " does not contain Atom"
        text += " " + an + ", so it may not be " + what
        raise FUNCTIONS.IPRO_IOError(text)
    # If Residue rn is a Design Position and Atom an is in the side chain
    if mn in experiment["Design Positions"] and \
    rn in experiment["Design Positions"][mn] and aa not in \
    backboneAtoms[molecule.fileFormat]:
        text = "Residue " + rn + " in Molecule " + mn + " is a Design Position "
        text += "and Atom " + an + " is not a backbone Atom, so it may not be "
        text += what
        raise FUNCTIONS.IPRO_IOError(text)

def fixedAtoms(fixed, experiment):
    """Check the definitions of the Atoms to always fix in place."""
    # Make sure fixed is a dictionary
    if not isinstance(fixed, dict):
        text = "The information about what Atoms must always be fixed in place "
        text += "must be stored in a dictionary."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store any remaining errors here
    errors = ''
    group = False
    molecule = False
    residue = False
    # Go through the Design Groups
    for gn in fixed:
        # Make sure this is a dictionary
        if not isinstance(fixed[gn], dict):
            # if the error message hasn't been said before
            if not group:
                errors += "\nThe information about what Atoms to fix in place "
                errors += "in a Design Group must be stored in a dictionary."
            # Make sure that error message isn't repeated
            group = True
            continue
        # Go through the Molecules
        for mn in fixed[gn]:
            # Check that it is a dictionary
            if not isinstance(fixed[gn][mn], dict):
                # If the Molecule error message hasn't been created before
                if not molecule:
                    errors += "\nThe information abut what Atoms to fix in "
                    errors += "place in a Molecule must be stored in a "
                    errors += "dictionary"
                molecule = True
                continue
            # Residues
            for rn in fixed[gn][mn]:
                if not isinstance(fixed[gn][mn][rn], list):
                    if not residue:
                        errors += "\nThe information about what Atoms to fix in"
                        errors += " place in a Residue must be stored in a list"
                    residue = True
                    continue
                # Atoms
                for an in fixed[gn][mn][rn]:
                    try:
                        Atom(gn, mn, rn, an, experiment, "fixed in place.")
                    except FUNCTIONS.IPRO_IOError as error:
                        errors += str(error)
                        continue
    # If there are errors, raise an error
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def position_restraint(restraint, experiment):
    """Check that a single position restraint is valid."""
    # Determine how many entries the restraint should have
    if experiment["Force Field"] == 'CHARMM' and experiment['File Format'] == \
    "PDB":
        N = 5
    else:
        text = "The " + experiment["Force Field"] + " force field and the "
        text += experiment["File Format"] + " file format combination is not "
        text += "supported by the IO.CHECK.position_restraint function."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check that the restraint is a list
    if not isinstance(restraint, list):
        text = "The information about a position restraint must be stored in a "
        text += "list."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the list is the wrong length
    if len(restraint) not in [N, N+1]:
        text = "When using the " + experiment["Force Field"] + " force field "
        text += "and " + experiment["File Format"] + " file format combination,"
        text += " a position restraint should contain " + str(N) + " items, not"
        text += ' ' + str(len(restraint))
        raise FUNCTIONS.IPRO_IOError(text)
    # Check that the Atom is validily accessible. But first check that the
    # Molecule name, Residue name, and Atom name aren't 'all' (which is allowed
    # - such restraints aren't applied to the side chains of Residues)
    if restraint[1] == 'all' and restraint[2] != 'all':
        text="When the Molecule specification in a position restraint is 'all',"
        text += " the Residue specification must also be 'all'"
        raise FUNCTIONS.IPRO_IOError(text)
    elif restraint[1] == 'all' and restraint[3] != 'all':
        text = "When the Molecule specification in a position restraint is "
        text += "'all', the Atom specification must also be 'all'"
        raise FUNCTIONS.IPRO_IOError(text)
    elif restraint[2] == 'all' and restraint[3] != 'all':
        text = "When the Residue specification in a position restraint is 'all'"
        text += ", the Atom specification must also be 'all'"
        raise FUNCTIONS.IPRO_IOError(text)
    # Only check the Atom if the atom specification isn't 'all'
    if restraint[3] != 'all':
        Atom(restraint[0], restraint[1], restraint[2], restraint[3],experiment,\
             "used in a position restraint.")
    # Check that the remaining entries are numbers
    for i in range(4, N):
        if not isinstance(restraint[i], float):
            text = "The " + str(i+1) + "th entry in a position restraint must "
            text += "be a floating point number."
            raise FUNCTIONS.IPRO_IOError(text)
    # If there is a last entry
    if len(restraint) == N+1:
        if not isinstance(restraint[-1], int) or restraint[0] == 'all' or \
        restraint[0] < restraint[-1]:
            text = "When a position restraint contains an extra entry, that "
            text += "entry MUST contain a Design Group number that is lower "
            text += "than the number of the Design Group specified at the "
            text += "beginning of the restraint."
            raise FUNCTIONS.IPRO_IOError(text)

def position_restraints(restraints, experiment):
    """Check a container of position restraints"""
    # Store any generated errors here
    errors = ''
    # If the container is not a list
    if not isinstance(restraints, list):
        errors += "\nPosition Restraints must be stored in a list"
    # Otherwise, check each restraint
    else:
        for restraint in restraints:
            try:
                position_restraint(restraint, experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # If there are errors, tell the user
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def distance_restraint(restraint, experiment):
    """Check the validity of a Distance Restraint"""
    # Determine how many entries should be in this restraint
    if experiment["Force Field"] == "CHARMM" and experiment["File Format"] == \
    "PDB":
        N = 8
    else:
        text = "The " + experiment["Force Field"] + " force field and "
        text += experiment["File Format"] + " combination is not supported by "
        text += "the IO.CHECK.distance_restraint function."
        raise FUNCTIONS.IPRO_IOError(text)
    # If the container is not a list
    if not isinstance(restraint, list):
        text = "A distance restraint must be stored in a list"
        raise FUNCTIONS.IPRO_IOError(text)
    # If the list isn't the right length
    if len(restraint) != N:
        text="When using the " + experiment["Force Field"] + " force field and "
        text += experiment["File Format"] + " file format, a distance restraint"
        text += " must contain " + str(N) + " entries."
        raise FUNCTIONS.IPRO_IOError(text)
    # If either of the Atom specifications isn't useful
    if not isinstance(restraint[1], list) or len(restraint[1]) != 3 or \
       not isinstance(restraint[2], list) or len(restraint[2]) != 3:
        text = "The Atom specifications in this distance restraint are wrong:\n"
        text += str(restraint)
        raise FUNCTIONS.IPRO_IOError(text)
    # Check both Atoms
    Atom(restraint[0], restraint[1][0], restraint[1][1], restraint[1][2], \
         experiment, "used in a distance restraint.")
    Atom(restraint[0], restraint[2][0], restraint[2][1], restraint[2][2], \
         experiment, "used in a distance restraint.")
    if restraint[1] == restraint[2]:
        text = "The same Atom is specified twice in this Dihedral Restraint:\n"
        text += str(restraint)
        raise FUNCTIONS.IPRO_IOError(text)
    # Check that the remaining entries are numbers
    for i in range(3, N):
        if not isinstance(restraint[i], float):
            text = "The " + str(i+1) + "th entry in a distance restraint must "
            text += "be a floating point number."
            raise FUNCTIONS.IPRO_IOError(text)

def distance_restraints(restraints, experiment):
    """Check all Distance Restraints"""
    # Store errors here
    errors = ''
    # Make sure it is a list
    if not isinstance(restraints, list):
        errors += "\nDistance Restraints must be stored in a list."
    else:
        for restraint in restraints:
            try:
                distance_restraint(restraint, experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def dihedral_restraint(restraint, experiment):
    """Check a Dihedral Restraint for validity."""
    # Determine how many entries it should have
    if experiment["Force Field"] == "CHARMM" and experiment["File Format"] == \
    "PDB":
        N = 7
    else:
        text = "The " + experiment["Force Field"] + " force field and "
        text += experiment["File Format"] + " file format combination is not "
        text += "supported by the IO.CHECK.dihedral_restraint function."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make sure it is a list
    if not isinstance(restraints, list) or len(restraints) != N:
        text = "When using the " + experiment["Force Field"]+" force field and "
        text += experiment["File Format"] + " file format combination, a "
        text += "Dihedral Restraint must be a list of " + str(N) + " items."
        raise FUNCTIONS.IPRO_IOError(text)
    # Check the Atom definitions
    for i in range(1, 5):
        if not isinstance(restraint[i], list) or len(restraint[i]) != 3:
            text = "This is not a valid Atom specification for a Dihedral "
            text += "Restraint:\n" + str(restraint[i])
            raise FUNCTIONS.IPRO_IOError(text)
        Atom(restraint[0], restraint[i][0], restraint[i][1], restraint[i][2], \
             experiment, "used in a dihedral restraint.")
        # Make sure no Atoms match
        for j in range(1, i):
            if restraint[j] == restraint[i]:
                text = "The same Atom is used more than once in this Dihedral "
                text += "Restraint:\n" + str(restraint)
                raise FUNCTIONS.IPRO_IOError(text)
    # Check that the remaining entries are numbers
    for i in range(5, N):
        if not isinstance(restraint[i], float):
            text = "The " + str(i+1) + "th entry in a Dihedral Restraint must "
            text += "be a floating point number."
            raise FUNCTIONS.IPRO_IOError(text)

def dihedral_restraints(restraints, experiment):
    """Check the validity of a set of Dihedral Angle Restraints"""
    errors = ''
    if not isinstance(restraints, list):
        errors += "\nDihedral Restraints must be stored in a list."
    else:
        for restraint in restraints:
            try:
                dihedral_restraint(restraint, experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def Restraints(restraints, experiment):
    """Check that a Restraints entry is acceptable."""
    # Make sure it is a dictionary
    if not isinstance(restraints, dict):
        error = "Structure restraints must be stored in a dictionary."
        raise FUNCTIONS.IPRO_IOError(error)
    # Store other errors here
    errors = ''
    # go through each restraint type in the dictionary
    for type in restraints:
        if type == "Fixed Atoms":
            try:
                fixedAtoms(restraints['Fixed Atoms'], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
        elif type == "Position":
            try:
                position_restraints(restraints["Position"], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
        elif type == "Distance":
            try:
                distance_restraints(restraints["Distance"], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
        elif type == "Dihedral":
            try:
                dihedral_restraints(restraints["Dihedral"], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
        else:
            errors += "\n" + type + " is not a valid structure restraint type"
    # If there are errors, tell teh user
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def structure(molecule):
    """Check a Molecule for gaps in its sequence"""
    # Store Residues that are missing Atoms here
    missing = []
    # Store Gaps that are too long here
    gaps = []
    # Search for different things based on file format
    if molecule.fileFormat == "PDB":
        n1 = 'C'
        n2 = 'N'
    else:
        text = "The CHECK structure function does not support the "
        text += molecule.fileFormat + " file format"
        raise FUNCTIONS.IPRO_IOError(text)
    # Loop through the Molecule's Residues
    for i in range(len(molecule) - 1):
        # If either this Residue or the next is not an amino acid, don't check
        # it
        if molecule[i].kind not in aminoAcids[molecule.fileFormat] or \
           molecule[i+1].kind not in aminoAcids[molecule.fileFormat]:
            continue
        # If the current Residue is missing the first Atom for the distance
        # calculation, store the Residue as having a structural problem
        if n1 not in molecule[i]:
            # If the Residue isn't already listed as being incomplete, do so
            if molecule[i].name not in missing:
                missing.append(molecule[i].name)
            continue
        # If the next Residues is missing the second Atom
        if n2 not in molecule[i+1]:
            missing.append(molecule[i+1].name)
            continue
        # Calculate the distance between the Atoms
        distance = MOLECULES.calculate_distance(molecule[i][n1], \
                                                molecule[i+1][n2])
        # If the Atoms are too far apart, store that
        if distance > 2.5:
            text = "Residues " + molecule[i].name + " and " + molecule[i+1].name
            text += " are " + format(distance, '.3f') + " angstroms apart."
            gaps.append(text)
    # Also identify an Residues that are possibly not supposed to be included in
    # IPRO
    hetatms = []
    for residue in molecule:
        # Check differently based on file format
        if residue.fileFormat == "PDB":
            if residue[0].kind == "HETATM":
                hetatms.append(residue)
    # Return the problem information
    return missing, gaps, hetatms
