#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Docking Functions"
# Documentation
__docking__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This module contains functions for carrying out the local, rigid-body docking
that is used in the IPRO Suite of Programs."""

# Include standard PYTHON modules
import os
import sys
import math
import random
# Also include IPRO Suite Modules
from STANDARDS import *
import ROTAMERS

class DockingError(IPRO_Error):
    """An error class for problems in the Docking Functions Module."""
    def __init__(self, error = ''):
        """The initialization function of the Docking Error class."""
        IPRO_Error.__init__(error)

def find_docking_group(dockingGroups, gn, mn):
    """Find the docking group that a particular Molecule is in."""
    # Loop through them
    for i in range(len(dockingGroups)):
        # If the Molecule is in this group
        if [gn, mn] in dockingGroups[i]:
            return i
    # Raise an error if it couldn't be found
    print dockingGroups
    print gn
    print mn
    text = "Somehow there is no docking group for Molecule " + mn + " in "
    text += "Design Group " + str(gn)
    raise DockingError(text)

def combine_docking_groups(dockingGroups, i, j):
    """Combine docking groups i and j."""
    # Store them in the lower numbered group
    if i < j:
        # Combine them
        dockingGroups[i].extend(dockingGroups[j])
        # And delete the extra
        del dockingGroups[j]
    elif j < i:
        dockingGroups[j].extend(dockingGroups[i])
        del dockingGroups[i]
    # If they're equal, nothing is done

def calculate_docking_groups(experiment):
    """Calculate all groups of Molecules that move together during docking"""
    # Initialize the groups as each Molecule being in an independent group
    dockingGroups = []
    # Go through the Design Groups and Molecules
    for group in experiment:
        for molecule in group:
            dockingGroups.append([[group.number, molecule.name]])
    print dockingGroups
    for group in experiment:
        print group.number
        print group
    # Go through the various types of restraints, and use them to combine
    # docking groups
    if "Restraints" in experiment:
        # The fixed Atoms specifications will be checked later. 
        # Check restraints on positions of individual atoms
        if "Position" in experiment["Restraints"]:
            for info in experiment["Restraints"]["Position"]:
                # This ONLY matters when the last value is an integer (instead
                # of a float), indicating a position restraint relative to
                # another Design Group
                if isinstance(info[-1], int):
                    i = find_docking_group(dockingGroups, info[0], info[1])
                    j = find_docking_group(dockingGroups, info[-1], info[1])
                    combine_docking_groups(dockingGroups, i, j)
        # Check restraints on distances between 2 Atoms
        if "Distance" in experiment["Restraints"]:
            # Go through the restraints
            for info in experiment["Restraints"]["Distance"]:
                # Get the appropriate Design Groups (it is possibly 'all')
                if info[0] == 'all':
                    gns = []
                    for group in experiment:
                        gns.append(group.number)
                else:
                    gns = [info[0]]
                # Loop through the Design Groups
                for gn in gns:
                    # Get the two Molecules involved
                    i = find_docking_group(dockingGroups, gn, info[1][0])
                    j = find_docking_group(dockingGroups, gn, info[2][0])
                    # And combine their groups
                    combine_docking_groups(dockingGroups, i, j)
        # Do the same thing for Dihedral restraints
        if "Dihedral" in experiment["Restraints"]:
            for info in experiment["Restraints"]["Dihedral"]:
                if info[0] == 'all':
                    gns = []
                    for group in experiment:
                        gns.append(group.number)
                else:
                    gns = [info[0]]
                for gn in gns:
                    # There are four Atoms, so make sure the first is in the
                    # same docking group as each of the others. That will ensure
                    # ALL are in the same group
                    i = find_docking_group(dockingGroups, gn, info[1][0])
                    for J in range(2, 5):
                        j = find_docking_group(dockingGroups, gn, info[J][0])
                        combine_docking_groups(dockingGroups, i, j)
    # Create a final list of docking groups that have Molecules that can
    # actually move during docking
    final = []
    for dockingGroup in dockingGroups:
        # Keep track of whether or not it can move
        use = True
        # Go through each Molecule in the docking group
        for mol in dockingGroup:
            # If the Molecule is a Design Molecule, the group may not move
            if experiment[mol[0]][mol[1]].design:
                use = False
                break
            # If the Molecule has any Residues with Atoms that are fixed
            # permanently in place, it may not move
            for residue in experiment[mol[0]][mol[1]]:
                if len(residue.fixedAtoms) > 0:
                    use = False
                    break
            if not use:
                break
        # If the docking group should move during docking, store it
        if use:
            final.append(dockingGroup)
    # Store the docking groups in the experiment
    experiment["Docking Groups"] = final

def collect_docking_groups(experiment, gn = None):
    """Collect the Molecules that will actually move during docking."""
    # Store the Molecules in this list
    dockingGroups = []
    # Loop through the Experiment's docking groups
    for group in experiment["Docking Groups"]:
        # Only keep groups that have all Molecules in allowed Design Groups
        use = True
        for molID in group:
            if gn not in [None, molID[0]]:
                use = False
                break
        # Store the group if it should be used
        if use:
            dockingGroups.append(group)
    return dockingGroups

def create_moving_structures(experiment, dockingGroups, fileNames):
    """Create the structures that move during docking."""
    # Start creating the docking_information.txt file's data with this string
    information = str(len(dockingGroups)) + "\n"
    # Store each molecule that moves during docking in this list
    movingMolecules = []
    # Loop through the docking groups
    for dockingGroup in dockingGroups:
        # Inlcude the number of Molecules in this docking group
        information += str(len(dockingGroup)) + "\n"
        # Go through each molecule
        for molID in dockingGroup:
            # Get that Molecule
            molecule = experiment[molID[0]][molID[1]]
            # Modify its Residues' permissions so that all Atoms will be output
            for residue in molecule:
                residue.permission = "FIXED"
            # Parameterize the Molecule
            ROTAMERS.parameterize(molecule, experiment)
            # Get the energy calculation - ready text
            contents = format(molecule, "energy")
            # Put those contents in a file
            fileName = "docking_" + str(molID[0]) + str(molID[1]) + ".txt"
            f = open(fileName, "w")
            f.write(contents)
            f.close()
            # Store that file name so it can be deleted later
            fileNames.append(fileName)
            # Store this Molecule in the moving Molecules list
            movingMolecules.append(molID)
            # Update information with data about this Molecule
            information += fileName + " " + str(molID[0]) + "\n"
    # Return what needs to be returned
    return information, movingMolecules

def create_static_structures(experiment, movingMolecules, fileNames):
    """Create the structures that don't move during docking."""
    # First, determine what Design Groups have structures that move during
    # docking
    designNumbers = []
    for molID in movingMolecules:
        if molID[0] not in designNumbers:
            designNumbers.append(molID[0])
    # Assemble the constant portions of each of those complexes
    for gn in designNumbers:
        # Store the contents of all of the Molecules in this list
        contents = ''
        # Loop through the Molecules in this Design Group
        for molecule in experiment[gn]:
            # If this Molecule moves during docking, don't consider it
            if [gn, molecule.name] in movingMolecules:
                continue
            # Set each Residue's permission to FIXED so all atoms are output
            for residue in molecule: 
                residue.permission = "FIXED"
            # Parameterize the Molecule
            ROTAMERS.parameterize(molecule, experiment)
            # Add its contents to contents
            contents += format(molecule, "energy")
        # Write these contents to a file, even if they're empty 
        fileName = "constant" + str(gn) + ".txt"
        f = open(fileName, "w")
        f.write(contents)
        f.close()
        # Store the file name so it can be deleted later
        fileNames.append(fileName)
    # Nothing needs to be returned from this function

def make_docking_information(experiment, dockingGroups, information, fileNames):
    """Make the file that describes how to run docking."""
    # Determine how long to run docking for
    if "Docking Iterations" in experiment:
        N = experiment["Docking Iterations"]
    else:
        N = defaultDockingIterations
    # Calculate the denominator for calculating the cooling schedule of
    # simulated annealing
    if N > 1:
        D = N - 1
    else:
        D = 1
    # The cooling schedule
    dT = (experiment["Docking Start Temp"] - experiment["Docking End Temp"])/D
    # Include the number of iterations of docking in information
    information += str(N) + "\n"
    # Include the data for each docking group for each iteration
    for i in range(N):
        # Calculate the Temperature for this docking iteration
        T = format(GasConstant*(experiment["Docking Start Temp"] - \
                                          i*dT), '.3f')
        # Loop through the docking groups
        for dockingGroup in dockingGroups:
            # Generate the random numbers for rotations
            for j in range(3): 
                information += format(math.radians(random.gauss(0, \
                               experiment["Docking SD Rotations"])), '.4f') +" "
            # Do the same for the perturbations
            for j in range(3): 
                information += format(random.gauss(0, \
                               experiment["Docking SD Movements"]), '.4f') + " "
            # Include the Temperature and a random number for simulated
            # annealing 
            information += T + " " + format(random.random(), '.4f') + "\n"
    # Write all of this information to a file
    fileName = "docking_information.txt"
    f = open(fileName, "w")
    f.write(information)
    f.close()
    # Store the file name
    fileNames.append(fileName)

def load_docking_results(experiment, dockingGroups, fileNames):
    """Load the results of docking."""
    # Try to open the docking results file
    try:
        fileName = "docking_results.txt"
        f = open(fileName, "r")
        # Since the file could be opened, store its name for deletion
        fileNames.append(fileName)
    # If that didn't work, return a positive integer so that an error can be
    # raised
    except IOError:
        text = "The docking_results.txt file was never created."
        raise DockingError(text)
    # Read in the results
    results = []
    for line in f:
        # There should be 8 numbers for each docking group: X rotation, Y
        # rotation, Z rotation, X perturbation, Y perturbation, Z perturbation,
        # final energy, initial energy
        items = line.split()
        # skip lines that don't have the right number of results
        if len(items) != 8:
            continue
        # Store the float of each result
        data = []
        for item in items:
            try:
                data.append(float(item))
            # If that didn't work, return 2 to indicate there was a problem with
            # a particular results
            except ValueError:
                f.close()
                text = "The contents of the docking_results.txt file are not "
                text += "properly formated."
                raise DockingError(text)
        # Store the data in results
        results.append(data)
    # Close the file
    f.close()
    # Make sure there are the proper number of results
    if len(results) != len(dockingGroups):
        text += "The docking_results.txt file does not have the right number "
        text += "of results."
        raise DockingError(text)
    # Create a message to store in the experiment's summary
    message = ''
    # Summarize the results for each docking group
    for I, dockingGroup in enumerate(dockingGroups):
        # Say what group this is
        message += "Docking Group " + str(I+1) + "\n"
        # Include each Molecule
        for molID in dockingGroup:
            message += "Molecule " + str(molID[1]) + " from Design Group "
            message += str(molID[0]) + "\n"
        # Include the results of how it was moved
        message += "Initial Energy: " + format(results[I][7], '.3f') + "\n"
        message += "Final Energy:   " + format(results[I][6], '.3f') + "\n"
        message += "X Rotation:     " + format(math.degrees(results[I][0]), \
                                               '.3f') + "\n"
        message += "Y Rotation:     " + format(math.degrees(results[I][1]), \
                                               '.3f') + "\n"
        message += "Z Rotation:     " + format(math.degrees(results[I][2]), \
                                               '.3f') + "\n"
        message += "X Perturbation: " + format(results[I][3], '.3f') + "\n"
        message += "Y Perturbation: " + format(results[I][4], '.3f') + "\n"
        message += "Z Perturbation: " + format(results[I][5], '.3f') + "\n"
    return message

def load_docking_structures(experiment, dockingGroups):
    """Load the structures of Molecules after docking."""
    # Loop through the docking groups
    for dockingGroup in dockingGroups:
        # Loop through the Molecules
        for molID in dockingGroup:
            # Generate the expected name for the Molecule's structure
            fileName = "docking_" + str(molID[0]) + str(molID[1]) + ".txt"
            # Try to open the file
            try:
                f = open(fileName, "r")
            except IOError:
                text = "The structure of a Molecule that moves during docking "
                text += "is missing."
                raise DockingError(text)
            # Determine the number of entries that are expected for describing
            # an Atom
            if experiment["File Format"] == "PDB":
                NUMBER = 8
            else:
                text = "The load_docking_structures function does not support "
                text += "the " + str(experiment["File Format"]) + " file format"
                raise ExperimentError(text)
            # Store the coordinates of the atoms in this dictionary
            data = {}
            # Keep track of the number of atoms found
            count = 0
            # Loop through the file's data
            for line in f:
                # Split the contents by white space
                items = line.split()
                if len(items) != NUMBER:
                    continue
                # Get the residue number, atom name, and coordinates
                try:
                    rn = int(items[1])
                    an = items[4]
                    x = float(items[5])
                    y = float(items[6])
                    z = float(items[7])
                # If any of the information isn't formatted properly
                except ValueError:
                    f.close()
                    text = "The contents of the file of a Molecule that moves "
                    text += "during docking are not formatted correctly."
                    raise DockingError(text)
                # Store the information about this Atom
                count += 1
                if rn not in data:
                    data[rn] = {}
                data[rn][an] = [x, y, z]       
            # Close the file
            f.close()                          
            # If there aren't the right number of entries, indicate a problem
            molecule = experiment[molID[0]][molID[1]] 
            if count != molecule.atomLength:
                text = "The file of a Molecule that moves during docking does "
                text += "not have the proper number of Atom entries."
                raise DockingError(text)
            # Store the coordinates of each Atom
            for residue in molecule:
                for atom in residue:
                    for i in range(3):
                        try:
                            atom[i] = data[residue.number][atom.name][i]
                        except (KeyError, IndexError, IPRO_Error):
                            text = "The file of a Molecule that moves during "
                            text += "docking contains wrong Atom entries."
                            raise DockingError(text)
