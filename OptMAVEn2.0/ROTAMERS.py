#!/usr/bin/env python

# The name of this module
__name__ = "IPRO Suite Rotamers Module"
# The file's documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for collecting, positioning, and choosing rotamers
for replacing the side chains of amino acid Residues."""

# Import necessary PYTHON modules
import os
import sys
import math
import copy
import time
# Import functions and classes from other IPRO suite modules that are needed.
# Start with the STANDARDS module, which contains all "default" variables, the
# "supported" lists, the aminoAcids and backboneAtoms dictionaries, and the
# IPRO_Error class
from STANDARDS import *
# From MOLECULES include the Residue, Molecule, and DesignGroup classes, as well
# as the MoleculeError class. Also provide access to the functions from
# MOLECULES
import MOLECULES
# Include the Experiment, GAMS, and Input / Output modules
import EXPERIMENT
import GAMS
import CPLEX
# Create an error class for this module
class RotamerError(IPRO_Error):
    """An Error for problems in the ROTAMERS module."""
    def __init__(self, error = ''):
        """The initialization function of the RotamerError class."""
        IPRO_Error.__init__(self, error)

# There are 3 main functions in this module - Closest_Rotamer (finds the most
# similar rotamer from the rotamer library for each Residue), Glycines (mutates
# Residues to glycine before Backbone Perturbations) and Optimal_Rotamers (which
# selects the best combination of rotamers after a Backbone Relaxation). This
# function validates the Residues given to those functions to make sure they are
# permissible.
def validate_residues(residues):
    """Make sure the input 'residues' to functions are allowed."""
    # Keep track of whether or not the Residues are acceptable. They need to be
    # in a nested iterable object that has Residues at the third level
    # (DesignGroup, Molecule, Residues)
    error = False
    # A Design Group is the expected input, and nothing needs to be done if that
    # is what is provided
    if isinstance(residues, MOLECULES.DesignGroup):
        pass
    # Put a Molecule in a list
    elif isinstance(residues, MOLECULES.Molecule):
        residues = [residues]
    # Put a single Residue in a list of a list
    elif isinstance(residues, MOLECULES.Residue):
        residues = [[residues]]
    # If it is a non-empty list, things are more complicated
    elif isinstance(residues, list) and len(residues) > 0:
        # Make a new list for the output of this function
        new = []
        # Loop through the list
        for obj in residues:
            # if it is a Molecule, store it
            if isinstance(obj, MOLECULES.Molecule):
                new.append(obj)
            # If it is a Residue, store it in a list
            elif isinstance(obj, MOLECULES.Residue):
                new.append([obj])
            # A list of lists of Residues is also permitted, so check to see if
            # this object is a list of Residues
            elif isinstance(obj, list) and len(obj) > 0:
                for obj2 in obj:
                    if not isinstance(obj2, MOLECULES.Residue):
                        error = True
                        break
                # If the object is a list of Residues, store it
                if not error:
                    new.append(obj)
            # If the input isn't one of those things, break the for loop
            else:
                error = True
                break
            # If there's an error, break the for loop
            if error:
                break
        # If the input has all been sorted properly, store residues as the new
        # list
        if not error:
            residues = new
    # If the input is anything else, say there's an error
    else:
        error = True
    # If the input is not permitted, raise an error
    if error:
        text = "The 'residues' input to a function in the ROTAMERS module must "
        text += "be either a DesignGroup, a Molecule, a Residue, or a list "
        text += "containing Molecules or Residues."
        raise RotamerError(text)
    # Make sure that all Residues use the proper file format. First get that
    # file format
    fileFormat = residues[0][0].fileFormat
    # Also make sure they use the same force field
    forceField = residues[0][0].forceField
    # Loop through the Residues, making sure they are numbered sequentially as
    # the loops proceed
    rn = 1
    an = 1
    for molecule in residues:
        for residue in molecule:
            # If the file format doesn't match up, raise an error
            if residue.fileFormat != fileFormat:
                text = "All Residues used by a function in the ROTAMERS module "
                text += "must use the same file format."
                raise RotamerError(text)
            elif residue.forceField != forceField:
                text = "All Residues used by a function in the ROTAMERS module "
                text += "must use the same force field."
                raise RotamerError(text)
            # Renumber the residues sequentially
            residue.renumber(an, rn)
            # Update the Atom and Residue numbers
            an = residue.lastAtom + 1
            rn = residue.number + 1
    # Return the residues
    return residues

def generate_library_path(residue, experiment = None):
    """Identify the proper folder in the rotamer library."""
    # This function is used internally in the ROTAMERS module, so it can be
    # 'safely' assumed that the input is a Residue class object
    # If either dihedral angle is missing, use the independ/ folder
    if residue.phi == None or residue.psi == None:
        end = "independ/"
    # Otherwise round the phi and psi dihedral angles to the nearest 10 and put
    # a sign (+, -) at the front of the string.
    else:
        phi = format(int(round(residue.phi, -1)), "+")
        psi = format(int(round(residue.psi, -1)), "+")
        end = phi + psi + "/"
    # Try to use the Experiment input
    try:
        path = experiment["Rotamer Library"] + end
    # If that can't be done, use the default
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        path = RotamerLibraryPath + end
    return path

def orient_residue(residue, how = False):
    """Move a Residue to a standard position.

    This function moves the Residue so that its CA is at the origin, its N is on
    the Z axis, and its CB (or HA1) in the XZ plane."""
    # Make sure the residue is a Residue
    if not isinstance(residue, MOLECULES.Residue):
        text = "The orient_residue function only works for Residue objects."
        raise RotamerError(text)
    # List the names of the expected Atoms for positioning
    if residue.fileFormat == "PDB":
        N1 = "CA"
        N2 = "N"
        if residue.kind == "GLY":
            N3 = "HA1"
        else:
            N3 = "CB"
    # If the file format isn't supported
    else:
        text = "The orient_residue function does not support the "
        text += str(residue.fileFormat) + " file format."
        raise RotamerError(text)
    # Confirm that the Residue has all of the needed Atoms
    if N1 not in residue or N2 not in residue or N3 not in residue:
        text = "This Residue is missing Atoms needed to orient it:\n"
        text += str(residue)
        raise RotamerError(text)
    # Make a new copy of the alpha carbon's coordinates. List does this
    coors = list(residue[N1])
    # Move the Residue so that Atom is at the origin
    MOLECULES.move(residue, coors, '-')
    # Rotate N into the XZ plane. Project it into the XY plane and calculate
    # an angle of rotation around the Z axis. Work out the trigonometry
    # yourself for why this works and explore how the acos function works to
    # determine why the if statement is necessary. Or just trust me that it
    # works.
    N = residue[N2]
    angle1 = math.fabs(math.acos(N['X'] /\
             math.sqrt(math.pow(N['X'], 2) + math.pow(N['Y'], 2))))
    if N['Y'] > 0:
        angle1 = -angle1
    rmatrix = MOLECULES.calculate_rmatrix(angle1, [0.0, 0.0, 1.0])
    MOLECULES.rotate(residue, rmatrix)
    # Now rotate the N onto the Z axis by rotating around the Y axis
    angle2 = math.fabs(math.acos(N['Z'] /\
             math.sqrt(math.pow(N['X'], 2) + math.pow(N['Z'], 2))))
    if N['X'] > 0:
        angle2 = -angle2
    rmatrix = MOLECULES.calculate_rmatrix(angle2, [0.0, 1.0, 0.0])
    MOLECULES.rotate(residue, rmatrix)
    # Rotate the first side chain atom into the XZ plane
    F = residue[N3]
    angle3 = math.fabs(math.acos(F['X'] /\
             math.sqrt(math.pow(F['X'], 2) + math.pow(F['Y'], 2))))
    if F['Y'] > 0:
        angle3 = -angle3
    rmatrix = MOLECULES.calculate_rmatrix(angle3, [0.0, 0.0, 1.0])
    MOLECULES.rotate(residue, rmatrix)
    # If this function has been run with a request to access how to position
    # a Rotamer into a Molecule's backbone, generate the list of reverse
    # movements and put the Residue back where it came from
    if how:
        reverse = []
        reverse.append(MOLECULES.calculate_rmatrix(-angle3, [0.0, 0.0, 1.0]))
        reverse.append(MOLECULES.calculate_rmatrix(-angle2, [0.0, 1.0, 0.0]))
        reverse.append(MOLECULES.calculate_rmatrix(-angle1, [0.0, 0.0, 1.0]))
        reverse.append(coors)
        unorient_residue(residue, reverse)
        return reverse

def unorient_residue(residue, movements):
    """Position a Residue in a Molecule's backbone."""
    # Assume no one would call this function on something other than a Residue
    try:
        # Rotate the first side chain atom out of the XZ plane
        MOLECULES.rotate(residue, movements[0])
        # Rotate the N off the Z axis
        MOLECULES.rotate(residue, movements[1])
        # Rotate the N out of the XZ plane
        MOLECULES.rotate(residue, movements[2])
        # Move the Residue
        MOLECULES.move(residue, movements[3], "+")
    except (IndexError, MOLECULES.MoleculeError):
        text = "The unorient_residue function has failed because of a "
        text += "problem with these inputs:\n" + str(movements)
        raise RotamerError(text)

def special_rotamers(residue, kind, require = False):
    """Retrieve special rotamers of the indicated kind for the Residue."""
    # Currently this function retrieves only a Proline Residue, but it could
    # theoretically be expanded in the future
    if residue.fileFormat == "PDB" and kind == "PRO":
        # Please review the IPRO suite paper for the specific details of how
        # this rotamer was generated. The Ramachandran values came from data in:
        # S.C. Lovell, et al., "Structure validation by Ca geometry: phi, psi
        # and Cb deviation," Proteins, 50:437-450 (2003)
        # The Proline rotamer itself is a repositioning of residue 14 in chain B
        # in 1CFQ.pdb, with hydrogens added by CHARMM c34b1
        # Keep track of whether or not to use the Proline rotamer
        use = False
        # If the user is requiring that the Proline residue be used, regardless
        # of Ramachandran clashes
        if require:
            use = True
        # If either phi or psi dihedral angle is missing, permit the rotamer
        elif residue.phi == None or residue.psi == None:
            use = True
        # If use hasn't been specified yet
        if not use:
            # Get the dihedral angles
            phi = int(round(residue.phi, -1))
            psi = int(round(residue.psi, -1))
            # Go through the permissible dihedral angles sets
            if phi == -40 and psi in [-60, -50, -40, -30, 110, 120]:
                use = True
            elif phi == -50 and psi in [-60, -50, -40, -30, 110, 120, 130, 140,\
            150]:
                use = True
            elif phi == -60 and psi in [-60, -50, -40, -30, -20, -10, 0, 100, \
            110, 120, 130, 140, 150, 160, 170]:
                use = True
            elif phi == -70 and psi in [-40, -30, -20, -10, 0, 10, 20, 30, 40, \
            50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170]:
                use = True
            elif phi == -80 and psi in [-180, 180, -170, -30, -20, -10, 0, 10, \
            20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, \
            170]:
                use = True
            elif phi == -90 and psi in [-180, 180, -20, -10, 0, 10, 20, 30, \
            120, 130, 140, 150, 160, 170]:
                use = True
            elif phi == -100 and psi in [-10, 0, 10, 20, 120, 130, 140, 150, \
            160, 170]:
                use = True
        # If the rotamer should not be gathered, return None
        if not use:
            return None
        # Otherwise, make a rotamer from this (it has already been oriented)
        text = """
ATOM      1  N   PRO P   1       0.000   0.000   1.451  1.00 52.87      MLP
ATOM      2  CD  PRO P   1       1.380  -0.040   1.994  1.00 52.80      MLP
ATOM      3  HD1 PRO P   1       0.768   0.827   2.331  1.00  0.00      MLP
ATOM      4  HD2 PRO P   1       2.021  -0.413   2.827  1.00  0.00      MLP
ATOM      5  CA  PRO P   1       0.000   0.000   0.000  1.00 51.90      MLP
ATOM      6  HA  PRO P   1      -0.455   0.913  -0.367  1.00  0.00      MLP
ATOM      7  CB  PRO P   1       1.468   0.000  -0.396  1.00 52.39      MLP
ATOM      8  HB1 PRO P   1       1.843   1.048  -0.400  1.00  0.00      MLP
ATOM      9  HB2 PRO P   1       1.643  -0.449  -1.394  1.00  0.00      MLP
ATOM     10  CG  PRO P   1       2.236   0.357   0.825  1.00 52.57      MLP
ATOM     11  HG1 PRO P   1       3.233   0.784   0.599  1.00  0.00      MLP
ATOM     12  HG2 PRO P   1       2.357  -0.550   1.460  1.00  0.00      MLP
ATOM     13  C   PRO P   1      -0.731  -1.192  -0.590  1.00 50.77      MLP
ATOM     14  O   PRO P   1      -0.717  -2.296  -0.043  1.00 50.66      MLP"""
        rot = MOLECULES.Residue(text[1:], 1, residue.number, residue.name, \
              residue.moleculeName, residue.forceField, residue.fileFormat)
        return rot
    # If the requested rotamer doesn't exist, just return None
    return None

def gather_rotamers(residue, experiment = None, procedure = "optimization", \
                    probability = False):
    """Collect all of the rotamers for a Residue."""
    # Currently, this function only supports the PDB file format, because that
    # is the format of the rotamer library. 
    if residue.fileFormat != "PDB":
        text = "The gather_rotamers function does not support the "
        text += str(residue.fileFormat) + " file format."
        raise RotamerError(text)
    # Make sure the procedure is recognized
    if procedure not in ['closest', 'glycine', 'optimization']:
        text = "The gather_rotamers function does not recognize the "
        text += str(procedure) + " procedure."
        raise RotamerError(text)
    # Determine if the Residue should actually receive a rotamer
    receive = True
    # Only amino acids should be mutated, Residues that are Fixed in place
    # shouldn't be modified, and Prolines that aren't being mutated
    if residue.kind not in aminoAcids[residue.fileFormat] or residue.permission\
    == "FIXED" or (residue.kind == "PRO" and residue.permission != "ROTAMER"):
        receive = False
    # Check to see if it has any non-backbone atoms that should be fixed in
    # place
    else:
        for name in residue.fixedAtoms:
            if name not in backboneAtoms[residue.fileFormat]:
                receive = False
                break
    # If rotamers should not be collected, return None
    if not receive:
        if probability:
            return None, None
        else:
            return None
    # Determine what types of amino acids the Residue may get rotamers of
    # If the procedure is just to find the most similar Rotamer to the current
    # side chain, use the Residue's current kind
    if procedure == "closest":
        kinds = [residue.kind]
    # If the Residue should be mutated to glycine before a backbone perturbation
    elif procedure == "glycine":
        kinds = ["GLY"]
    elif procedure == "optimization":
        # If the Residue is eligible for mutation
        if residue.permission == "ROTAMER" and residue.design:
            # Use the list of permitted mutations
            kinds = residue.permittedKinds
            # The default setting for that list for design Residues is NOT an
            # empty list. If it is empty, something odd is going on and simply
            # prevent the Residue from mutating
            if kinds == []:
                kinds = [residue.currentKind]
        else:
            # If it isn't being mutated, only use amino acids that match the
            # Residue's current kind attribute
            kinds = [residue.currentKind]
    # Store the rotamers (and probabilities) in a list
    rotamers = []
    if probability:
        probabilities = []
    # Get the movements needed to position a rotamer
    movements = orient_residue(residue, True)
    # Get the path to the proper portion of the rotamer library
    path = generate_library_path(residue, experiment)
    # Gather rotamers
    for kind in kinds:
        # Adjust for histidine
        if kind in ['HIS']:
            if "HSD" not in kinds:
                kind = "HSD"
            else:
                continue
        # There are an unknown number of rotamers
        for i in range(1, 1001):
            # The rotamers are numbered sequentially
            try:
                f = open(path + kind.lower() + str(i) + ".pdb", "r")
            except IOError:
                # If the file couldn't be opened then there are no more rotamers
                # of that type
                break
            lines = f.readlines()
            f.close()
            # Make and position a rotamer
            rot = MOLECULES.Residue(lines, 1, residue.number, residue.name, \
                  residue.moleculeName, residue.forceField, residue.fileFormat)
            orient_residue(rot)
            unorient_residue(rot, movements)
            # Store the rotamer
            rotamers.append(rot)
            # If appropriate, get probabilities
            if probability:
                try:
                    f = open(path + kind.lower() + str(i) + ".prob", "r")
                    line = f.readline()
                    f.close()
                    p1 = -math.log(float(line.split()[0]))
                    p2 = -math.log(float(line.split()[1]))
                # not all rotamers (e.g. GLY, ALA) have probabilities. Use 0 and
                # 0 for those values
                except IOError:
                    p1 = 0.0
                    p2 = 0.0
                probabilities.append([p1, p2])
    # Try to retrieve special rotamers
    specialKinds = ["PRO"]
    for kind in specialKinds:
        # If the Residue should not receive rotamers of this kind
        if kind not in kinds:
            continue
        # If the Residue already has at least one rotamer of this kind
        flag = False
        for rot in rotamers:
            if rot.kind == kind:
                flag = True
                break
        if flag:
            continue
        # Determine if the rotamer MUST be used
        if len(kinds) == 1:
            require = True
        else:
            require = False
        # Get the rotamer
        rot = special_rotamers(residue, kind, require)
        # If it is None, don't do anything
        if rot == None:
            continue
        # position it
        unorient_residue(rot, movements)
        # Store it
        rotamers.append(rot)
        # If appropriate, store probability values
        if probability:
            probabilities.append([0.0, 0.0])
    # If no rotamers were collected, raise an error
    if len(rotamers) == 0:
        text = "No rotamers were gathered for this Residue:\n" + str(residue)
        text += "File Format: " + str(residue.fileFormat)
        text += "\nForce Field: " + str(residue.forceField)
        text += "\nDesign Status: " + str(residue.design)
        text += "\nPhi Dihedral: " + format(residue.phi, '.3f')
        text += "\nPsi Dihedral: " + format(residue.psi, '.3f')
        raise RotamerError(text)
    if probability:
        return rotamers, probabilities
    else:
        return rotamers

def rotamer_patcher(residue, rotamer):
    """Put the side chain of a rotamer into a Residue."""
    # Assume no one would be foolish enough to mis-use this function. 
    # have a different method for each supported file format
    if residue.fileFormat == "PDB":
        # Create an appropriate end attribute for the Atoms
        end = "  1.00  0.00      ML" + residue.moleculeName.upper()
        # Make sure the Residue ends up properly numbered
        atomNum = residue[0].number
        resNum = residue.number
        # Modify the kind of Residue to match the Rotamer. This automatically
        # updates the Residue's Atoms
        residue.kind = rotamer.kind
        # Gather the backbone atoms from the Residue
        atoms = []
        for atom in residue:
            # Doing things this way instead of residue["backbone"] allows
            # non-canonical amino acids to be mutated, although getting rotamers
            # for them will be difficult
            if atom.name in backboneAtoms[atom.fileFormat]:
                # Make sure alpha carbon hydrogens are properly named on
                # glycines
                if residue.kind == "GLY" and atom.name == "HA":
                    atom.name = "HA2"
                elif atom.name == "HA2" and residue.kind != "GLY":
                    atom.name = "HA"
                # Skip HN atoms when mutating to Proline
                elif residue.kind == "PRO" and atom.name == "HN":
                    continue
                # Do something similar for N-terminal Proline mutations
                elif residue.kind == "PRO" and atom.name in ['HT1','HT2','HT3']:
                    # Skip the third one (bond taken up by the CD of the PRO)
                    if atom.name == "HT3":
                        continue
                    else:
                        atom.name = "HN" + atom.name[2]
                # Do the opposite when mutating away from N-terminal prolines
                elif atom.name in ['HN1', 'HN2'] and residue.kind != "PRO":
                    atom.name = "HT" + atom.name[2]
                # Update the Atom's end attribute
                atom.end = end
                # Store the Atom in the list of Residues
                atoms.append(atom)
        # Loop through the rotamer and get its side chain Atoms
        for atom in rotamer:
            if atom.name not in backboneAtoms[rotamer.fileFormat]:
                # Side chain atoms won't ever need to be renamed, so just update
                # the Atom's end attribute
                atom.end = end
                atoms.append(atom)
        # Have the Residue load the new Atoms
        residue.load(atoms, atomNum, resNum)
    # If the file format is not permitted
    else:
        text = "The rotamer_patcher function does not support the "
        text += str(residue.fileFormat) + " file format."
        raise RotamerError(text)

# Have functions that gather the non-bonded energy parameters for non-bonded
# energy calculations
def validate_files(files, kind):
    """Make sure a list of files is a list of strings without spaces."""
    # Keep track of whether the list is acceptable
    error = False
    # If it isn't a list or is empty, that's bad
    if not isinstance(files, list) or len(files) == 0:
        error = True
    # If it is a non-empty list, check each entry
    else:
        for f in files:
            if not isinstance(f, str) or ' ' in f or '\t' in f or '\n' in f:
                error = True
                break
    # If the list doesn't match expected patterns
    if error:
        text = "The following is not an acceptable list of " + kind + "files:\n"
        text += str(files)
        raise RotamerError(text)

def open_nb_file(path, fileName, fileType):
    """Attempt to open the indicated file."""
    # First, try and open it in an input_files folder
    try:
        f = open(path + "input_files/" + fileName, "r")
    except IOError:
        # Next try to open it in the current folder
        try:
            f = open(path + fileName, "r")
        except IOError:
            # Raise an error
            text = "The " + fileName + " " + fileType + " file could not be "
            text += "opened in this folder:\n" + path 
            raise RotamerError(text)
    return f

def get_CHARMM_topology(path, experiment = None):
    """Read in CHARMM topology file information."""
    # Get the names of the files
    try:
        files = experiment["CHARMM Topology Files"]
    # use defaults if necessary
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        files = defaultCHARMMTopologies
    # Validate the files
    validate_files(files, "CHARMM Topology")
    # Store the gathered information in this dictionary
    information = {}
    for fileName in files:
        # Open the file
        f = open_nb_file(path, fileName, "CHARMM Topology")
        res = False
        for line in f:
            # If it starts the declaration of a Residue
            if line.startswith(("RESI", "PRES")):
                try:
                    # Get the name of the Residue and create an entry in the
                    # information dictionary for it
                    res = line.split()[1]
                    if res not in information:
                        information[res] = {}
                # If the name of the Residue is missing, don't store any
                # following Atoms
                except IndexError:
                    res = False
            elif res != False and line.startswith("ATOM"):
                # Get the information about the Atom
                try:
                    items = line.split()
                    information[res][items[1]] = {"type":items[2], \
                                                  "Charge":float(items[3])}
                # If there is a problem getting the information about the Atom,
                # don't
                except (IndexError, ValueError):
                    pass
        # Close the file
        f.close()
    return information

def get_CHARMM_parameter(path, experiment = None):
    """Read in CHARMM parameter file information."""
    # Get the names of the files
    try:
        files = experiment["CHARMM Parameter Files"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        files = defaultCHARMMParameters
    # Validate them
    validate_files(files, "CHARMM Parameter")
    # Store the VDW information
    vdw = {}
    for fileName in files:
        f = open_nb_file(path, fileName, "CHARMM Parameter")
        # Have a flag for when the non-bonded portion of the file is reached
        flag = False
        for line in f:
            line = line.strip()
            # if the line is empty, skip it
            if len(line) == 0:
                continue
            # If it is a non-bonded parameter line
            if flag and not line.upper().startswith(("!", "CUTNB", "END", \
            "HBOND")):
                items = line.split()
                try:
                    # Get the information about the Atom
                    vdw[items[0]]={"eps":float(items[2]),"rad":float(items[3])}
                    # Try to get the 1-4 parameters, as well
                    try:
                        if items[4] == "!":
                            raise IndexError
                        vdw[items[0]]['14E'] = float(items[5])
                        vdw[items[0]]['14R'] = float(items[6])
                    except (IndexError, ValueError):
                        # If they can't be gathered (not all Atoms have them),
                        # then store them as the standard eps and rad
                        vdw[items[0]]['14E'] = vdw[items[0]]['eps']
                        vdw[items[0]]['14R'] = vdw[items[0]]['rad']
                except (IndexError, ValueError):
                    # Or skip it if there is a problem
                    pass
            # If the line starts the non-bonded parameters section
            elif line.upper().startswith(("NONBONDED", "NBONDED")):
                flag = True
        f.close()
    return vdw

def get_LK_solvation(path, experiment = None):
    """Get Lazaridis-Karplus non-bonded solvation information."""
    # Get the file names
    try:
        files = experiment["LK Solvation Files"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        files = defaultLKSolvationFiles
    # Validate them
    validate_files(files, "LK Solvation")
    # Store the information
    lk = {}
    for fileName in files:
        f = open_nb_file(path, fileName, "LK Solvation")
        for line in f:
            # Split the line on whitespace
            items = line.split()
            # Skip lines that shouldn't be used
            if len(items) != 6 or items[0] == "ATOM":
                continue
            # Get the Residue and store it in lk
            res = items[1]
            if res not in lk:
                lk[res] = {}
            # Store the information
            try:
                lk[res][items[0]] = {"LK Radius":float(items[2]), "LK Lambda":\
                    float(items[3]), "LK Volume":float(items[4]), \
                    "LK Gibbs Free Energy":float(items[5])}
            except (IndexError, ValueError):
                pass
        f.close()
    return lk

def get_nb_parameters(experiment = None):
    """Retrieve the non-bonded parameters from files."""
    # Determine the force field that is being used
    try:
        forceField = experiment["Force Field"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        forceField = defaultField
    # have a unique method for each supported force field
    if forceField == "CHARMM":
        # Identify where the parameter files should be located
        try:
            path = experiment["Folder"]
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            path = "./"
        # Get the topology file information
        information = get_CHARMM_topology(path, experiment)
        # And VDW from parameter files
        vdw = get_CHARMM_parameter(path, experiment)
        # Store the VDW information in information
        for res in information:
            for atom in information[res]:
                type = information[res][atom]["type"]
                # If there is VDW information for this atom type
                if type in vdw:
                    information[res][atom]["VDW Radius"] = vdw[type]["rad"]
                    information[res][atom]["VDW Epsilon"] = vdw[type]["eps"]
                    information[res][atom]["VDW 1-4 Rad"] = vdw[type]["14R"]
                    information[res][atom]["VDW 1-4 Eps"] = vdw[type]["14E"]
        # Determine what type of solvation, if any, is being used
        try:
            if experiment["Use Solvation"]:
                solvation = experiment["Solvation Type"]
            else:
                solvation = False
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            if defaultUseSolvation:
                solvation = defaultSolvationType
            else:
                solvation = False
        # LK solvation
        if solvation == "Lazaridis-Karplus":
            data = get_LK_solvation(path, experiment)
            # Store the data in information
            for res in data:
                for atom in data[res]:
                    if res in information and atom in information[res]:
                        for par in data[res][atom]:
                            information[res][atom][par] = data[res][atom][par]
        # If this is an Experiment, store the information in it
        if isinstance(experiment, (EXPERIMENT.Experiment, dict)):
            experiment["Nonbonded Parameters"] = information
        return information
    # If the force field is not supported
    else:
        text = "The get_nb_parameters function does not support the "
        text += str(forceField) + " force field."
        raise RotamerError(text)

def parameterize_Residue(residue, experiment = None, N = False, C = False):
    """Assign non-bonded energy parameters to the Residue's Atoms."""
    # Confirm that residue is a Residue
    if not isinstance(residue, MOLECULES.Residue):
        text = "The parameterize_Residue function requires a Residue class "
        text += "object as its 'residue' input."
        raise RotamerError(text)
    # Have a different method for each combination of force field and file
    # format
    if residue.forceField == "CHARMM" and residue.fileFormat == "PDB":
        # Determine what non-bonded energies to use
        try:
            energies = experiment["CHARMM Energy Terms"]
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            energies = defaultCHARMMEnergies
        # Determine what type of solvation is being used
        try:
            if experiment["Use Solvation"]:
                solvation = experiment["Solvation Type"]
            else:
                solvation = False
        # If settings haven't been explicitly stated, use the defaults
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            if defaultUseSolvation:
                solvation = defaultSolvationType
            else:
                solvation = False
        # Get the non-bonded energy parameters
        try:
            information = experiment["Nonbonded Parameters"]
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            information = get_nb_parameters(experiment)
        # Keep a record of whether or not any errors occurred while
        # parameterizing the Residue
        errorText = ''
        # Make sure the N and C terminus flags are valid by checking if atoms
        # that ONLY appear in those Residues exist in this Residue
        if not N and residue.kind in aminoAcids[residue.fileFormat]:
            for name in ['HT1', 'HT2', 'HT3', 'HN1', 'HN2']:
                if name in residue:
                    N = True
                    break
        if not C and residue.kind in aminoAcids[residue.fileFormat]:
            if "OT1" in residue or "OT2" in residue:
                C = True
        # Loop through the Residue's Atoms
        for atom in residue:
            # Generate a tag for accessing the proper non-bonded parameters
            if N and residue.kind == "GLY" and atom.name in ['N', 'HT1', 'HT2',\
            'HT3', 'CA', 'HA1', 'HA2']:
                res = "GLYP"
            elif N and residue.kind == "PRO" and atom.name in ['N', 'HN1', \
            'HN2', 'CA', 'HA', 'CD', 'HD1', 'HD2']:
                res = "PROP"
            elif N and residue.kind in aminoAcids[residue.fileFormat] and \
            atom.name in ['N', 'HT1', 'HT2', 'HT3', 'CA', 'HA']:
                res = "NTER"
            elif C and residue.kind in aminoAcids[residue.fileFormat] and \
            atom.name in ['C', 'OT1', 'OT2']:
                res = "CTER"
            elif residue.kind == "HIS":
                res = "HSD"
            else:
                res = residue.kind
            # Assign the appropriate non-bonded energy parameters
            if "vdw" in energies:
                # Try to assign the VDW attributes of the Atom
                try:
                    atom.vdwRadius = information[res][atom.name]["VDW Radius"]
                    atom.vdwEpsilon = information[res][atom.name]["VDW Epsilon"]
                    atom.vdw14R = information[res][atom.name]["VDW 1-4 Rad"]
                    atom.vdw14E = information[res][atom.name]["VDW 1-4 Eps"]
                    if experiment is not None:
                        atom.ljPhi = experiment["Lj Phi"] 
                # If that fails, store that error information
                except (KeyError, MOLECULES.MoleculeError):
                    errorText += "Residue: " + res.ljust(5) + "Atom: " 
                    errorText += atom.name.ljust(5) + "Parameters: VDW\n"
            # Similarly, if electrostatics is being used, try to store the value
            if "elec" in energies:
                try:
                    atom.charge = information[res][atom.name]["Charge"]
                except (KeyError, MOLECULES.MoleculeError):
                    #errorText += "Residue: " + res.ljust(5) + "Atom: "
                    errorText += "Residue: " + atom.residueName + "Atom: "
                    errorText += atom.name.ljust(5) + "Parameters: Charge\n"
            # Currently, Lazaridis-Karplus implicit solvation is the only
            # supported solvation type
            if solvation == "Lazaridis-Karplus":
                try:
                    atom.lkRadius = information[res][atom.name]["LK Radius"]
                    atom.lkLambda = information[res][atom.name]["LK Lambda"]
                    atom.lkVolume = information[res][atom.name]["LK Volume"]
                    atom.lkGibbs = information[res][atom.name]["LK Gibbs Free Energy"]
                except (KeyError, MOLECULES.MoleculeError):
                    errorText += "Residue: " + res.ljust(5) + "Atom: "
                    errorText += atom.name.ljust(5)+"Parameters: LK Solvation\n"
        return errorText
    # if the force field / file format combination is not supported
    else:
        text = "The parameterize_Residue function does not support the "
        text += str(residue.fileFormat) + " file format and "
        text += str(residue.forceField) + " force field combination."
        raise RotamerError(text)

def parameterize(residues, experiment = None):
    """Assign non-bonded energy parameters to a group of Residues."""
    # Validate the residues
    residues = validate_residues(residues)
    # Keep track of whether or not any errors occurred
    errorText = ''
    # Loop through at the molecule level
    for molecule in residues:
        # Loop through the Residues
        for residue in molecule:
            # Determine if this Residue is the first or last in a Molecule
            N = False
            C = False
            if isinstance(molecule, MOLECULES.Molecule):
                if residue.name == molecule[0].name:
                    N = True
                if residue.name == molecule[-1].name:
                    C = True
            # Parameterize the Residue
            errorText += parameterize_Residue(residue, experiment, N, C)
    # If there has been an error, raise it
    if errorText != '':
        text = "Assignment of the following non-bonded energy parameters has "
        text += "failed:\n" + errorText
        raise RotamerError(text)

def calculate_rc_energies(residues, rotamers, probability = False, \
                          probabilities = []):
    """Calculate the energies between rotamers and the rest of a system."""
    # Output the rotamers to a file. Store the output in this string
    output = ''
    # Keep track of the rotamer number of each rotamer
    resNum = None
    rotNum = 0
    # Keep track of which Residues have rotamers
    spots = []
    for rotamer in rotamers:
        # If this is the first rotamer for a new Residue, update the numbering
        # information
        if rotamer.number != resNum:
            resNum = rotamer.number
            rotNum = 0
            spots.append([rotamer.moleculeName, rotamer.name])
        # Increment the rotamer number
        rotNum += 1
        # Include the rotamer's side chain atoms in the output
        output += format(rotamer, "sidechain - " + str(rotNum))
    # Write the rotamer output to the appropriate file
    rotamerFile = "ROTAMERS.txt"
    f = open(rotamerFile, "w")
    f.write(output)
    f.close()
    # Now output the appropriate constant portions of the system
    output = ''
    # Loop through the Residues
    for molecule in residues:
        for residue in molecule:
            # If this Residue has at least one rotamer, only output its backbone
            # Atoms
            if [residue.moleculeName, residue.name] in spots:
                output += format(residue, "backbone - 1")
            # Otherwise output all of its Atoms
            else:
                output += format(residue, "all - 1")
    # Write all of this to a file
    constantFile = "CONSTANT.txt"
    f = open(constantFile, "w")
    f.write(output)
    f.close()
    # Run the CPP program that calculates the energies
    i = os.system("./rotamer_constant.out")
    # If everything worked fine, read in the results
    if i == 0:
        resultsFile = "RC_ENERGIES.txt"
        energies = []
        try:
            # Open the file
            f = open(resultsFile, "r")
            # Read in the energies
            for line in f:
                items = line.split()
                # Skip lines that don't have the proper number of entries
                if len(items) != 3:
                    continue
                # Get the information
                energies.append(float(items[2]))
            # Close the file
            f.close()
            # If the number of energies is not correct
            if len(energies) != len(rotamers):
                i = 1
        # If there was an error
        except (IOError, ValueError):
            i = 1
    # If there was a problem
    if i != 0:
        text = "There was a problem using C++ to calculate the energies of "
        text += "rotamers with the constant portions of the system of "
        text += "Molecules. Good luck solving this error."
        raise RotamerError(text)
    # If appropriate, modify the values of the energies with probability values
    if probability and len(probabilities) == len(energies):
        for i in range(len(probabilities)):
            energies[i] += probabilities[i][0] + probabilities[i][1]
            # Also modify the energies by statistical terms for the type of
            # amino acid, but have a flag to turn this off if wanted
            flag = True
            if flag and rotamers[i].fileFormat == "PDB":
                if rotamers[i].kind == "ALA":
                    energies[i] -= 0.05
                elif rotamers[i].kind == "CYS":
                    pass # The energy for CYS is 0.0
                elif rotamers[i].kind == "ASP":
                    energies[i] -= 1.43
                elif rotamers[i].kind == "GLU":
                    energies[i] -= 1.44
                elif rotamers[i].kind == "PHE":
                    energies[i] += 1.20
                elif rotamers[i].kind == "GLY":
                    energies[i] -= 0.37
                elif rotamers[i].kind in ["HIS", "HSD"]:
                    energies[i] += 1.92
                elif rotamers[i].kind == "ILE":
                    energies[i] += 0.45
                elif rotamers[i].kind == "LYS":
                    energies[i] -= 1.30
                elif rotamers[i].kind == "LEU":
                    energies[i] += 1.62
                elif rotamers[i].kind == "MET":
                    energies[i] -= 0.25
                elif rotamers[i].kind == "ASN":
                    energies[i] -= 0.17
                elif rotamers[i].kind == "PRO":
                    pass # The rotamer library doesn't have PRO, so use 0
                elif rotamers[i].kind == "GLN":
                    energies[i] -= 0.14
                elif rotamers[i].kind == "ARG":
                    energies[i] -= 0.60
                elif rotamers[i].kind == "SER":
                    energies[i] -= 0.14
                elif rotamers[i].kind == "THR":
                    energies[i] -= 0.31
                elif rotamers[i].kind == "VAL":
                    energies[i] += 0.25
                elif rotamers[i].kind == "TRP":
                    energies[i] += 1.85
                elif rotamers[i].kind == "TYR":
                    energies[i] += 1.08
    # Remove the files made by this function
    for fileName in [rotamerFile, constantFile, resultsFile]:
        try:
            os.remove(fileName)
        except OSError:
            pass
    # Return the calculated energies
    return energies

def sort_rotamers(rotamers, energies, experiment = None):
    """Retain only the most promising rotamers.

    Due to memory limitations, there is a maximum number of rotamers that can be
    considered at once. Use this function to collect those rotamers once the
    rotamer-constant energy calculations have been completed."""
    # First, sort and store the rotamers by position and energy
    sorted = {}
    for i in range(len(rotamers)):
        if rotamers[i].number not in sorted:
            sorted[rotamers[i].number] = []
        sorted[rotamers[i].number].append([rotamers[i], energies[i]])
    # Actually do the sorting now
    order = sorted.keys()
    order.sort()
    for spot in order:
        flag = False
        while not flag:
            flag = True
            for i in range(len(sorted[spot]) - 1):
                if sorted[spot][i][1] > sorted[spot][i+1][1]:
                    flag = False
                    temp = sorted[spot][i]
                    sorted[spot][i] = sorted[spot][i+1]
                    sorted[spot][i+1] = temp
    # Keep a minimum of N of each kind of rotamer, as well as all rotamers that
    # are within 'window' (in kcal/mol) of the lowest energy rotamer at each
    # position. Get that window now
    try:
        WINDOW = experiment["Rotamer Window"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        WINDOW = defaultRotamerWindow
    try:
        MAX = experiment["Max Rotamer Number"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        MAX = defaultMaxRotamers
    # In systems where the same protein chain is present more than once (e.g.
    # dimers), it is possible that MANY positions may be assigned rotamers. In
    # those cases, the maximum number of permissible rotamers must be decreased
    # to account for the presence of MANY positions being considered.
    if isinstance(experiment, (EXPERIMENT.Experiment, dict)) and "Dimers" in \
    experiment:
        # Get the names of all Molecules that have rotamers
        moleculeNames = []
        for rotamer in rotamers:
            if rotamer.moleculeName not in moleculeNames:
                moleculeNames.append(rotamer.moleculeName)
        # If any pair of dimers BOTH have rotamers, shrink the maximum number of
        # rotamers
        for pair in experiment["Dimers"]:
            if pair[0] in moleculeNames and pair[1] in moleculeNames:
                MAX -= 300
                break
    # N is the minimum number of each kind of amino acid to keep for each
    # position
    N = 3
    total = 2*MAX
    window = copy.deepcopy(WINDOW)
    # Identify when the while loop is starting
    start = time.time()
    # Keep only the best rotamers
    while total > MAX:
        # Reset the total to 0
        total = 0
        # Make lists of the chosen rotamers and their energies
        chosen_rotamers = {}
        chosen_energies = {}
        # Loop through the various positions
        for n in order:
            # Initialize the dictionaries
            chosen_rotamers[n] = []
            chosen_energies[n] = []
            # Determine what kinds of amino acids there are rotamers of
            kinds = []
            for group in sorted[n]:
                if group[0].kind not in kinds:
                    kinds.append(group[0].kind)
            # Keep track of which rotamers have been retained
            kept = []
            # Keep N of each kind of rotamer
            for kind in kinds:
                # Count how many of this type have been kept
                count = 0
                for i in range(len(sorted[n])):
                    if sorted[n][i][0].kind == kind:
                        chosen_rotamers[n].append(sorted[n][i][0])
                        chosen_energies[n].append(sorted[n][i][1])
                        kept.append(i)
                        count += 1
                        if count == N:
                            break
            # Now keep all rotamers within 'window' of the best
            for i in range(len(sorted[n])):
                # If that Rotamer has already been kept, skip it
                if i in kept:
                    continue
                # If it has an acceptable energy, store it
                elif sorted[n][i][1] <= sorted[n][0][1] + window:
                    chosen_rotamers[n].append(sorted[n][i][0])
                    chosen_energies[n].append(sorted[n][i][1])
                # Since the rotamers are sorted by energy, once one is not
                # acceptable no more that are will be found
                else:
                    break
            # Add the number of kept rotamers to the total
            total += len(chosen_rotamers[n])
        # If the while loop will not exit, modify the search parameters
        if total > MAX:
            # If window has been shrunk a lot, decrease N
            if window < 0.15*WINDOW and N > 1:
                N -= 1
                window = copy.deepcopy(WINDOW)
            # Otherwise just keep shrinking the window
            else:
                window = 0.95*window
        # If this while loop has been running for a long time, raise an error
        # because there is some problem with this system
        if time.time() - start > 300:
            text = "The sort_rotamers function is stuck in a while loop. It is "
            text += "impossible for it to shrink the number of rotamers "
            text += "sufficiently to make the problem computationally "
            text += "feasible. Good luck adjusting your system to fix this "
            text += "problem."
            raise RotamerError(text)
    # Return the chosen rotamers and energies
    return chosen_rotamers, chosen_energies

def calculate_rr_energies(rotamers):
    """Calculate the energies between rotamers."""
    # Create an output for the rotamers
    output = ''
    order = rotamers.keys()
    order.sort()
    for spot in order:
        for i, rotamer in enumerate(rotamers[spot]):
            output += format(rotamer, "sidechain - " + str(i+1))
    rotamerFile = "ROTAMERS.txt"
    f = open(rotamerFile, "w")
    f.write(output)
    f.close()
    # Run the CPP program
    i = os.system("./rotamer_rotamer.out")
    # Try to load the energies
    if i == 0:
        energies = []
        resultsFile = "RR_ENERGIES.txt"
        try:
            f = open(resultsFile, "r")
            for line in f:
                items = line.split()
                if len(items) != 5:
                    continue
                energies.append([int(items[0]), int(items[1]), int(items[2]), \
                                 int(items[3]), float(items[4])])
            f.close()
        # If there is a problem, just change the value of 
        except (IOError, ValueError):
            i = 1
    # If there was a problem raise an error
    if i != 0:
        text = "There was a problem using C++ to calculate the energies between"
        text += " rotamers. Good luck figuring this problem out."
        raise RotamerError(text)
    # Remove the files made by this function
    for fileName in [rotamerFile, resultsFile]:
        try:
            os.remove(fileName)
        except OSError:
            pass
    # Return the calculated energies
    return energies

def is_hydrogen(atom):
    """Determine if an Atom is a hydrogen or not."""
    # Have a different method for each supported file format, as the naming
    # scheme for Atoms may change
    if atom.fileFormat == "PDB":
        # If the first letter in the Atom's name is an H, it is a hydrogen
        for char in atom.name:
            # Check to see if it is a letter
            if char.isalpha():
                if char == "H":
                    return True
                else:
                    return False
        # If for some reason there were no letters, something is odd with the
        # Atom's name. Still, just return False and be fine with it
        return False
    # If the file format isn't supported, raise an error
    else:
        text = "The is_hydrogen function doesn't support the "
        text += str(atom.fileFormat) + " file format."
        raise RotamerError(text)

# Create the 3 functions that should actually be called from outside of this
# module
def Closest_Rotamer(residues, experiment = None):
    """Identify the most similar rotamer to each Residue."""
    # Validate the input Residues
    residues = validate_residues(residues)
    # Loop through the Molecules
    for molecule in residues:
        # If it is a Molecule, calculate its dihedral angles.
        if isinstance(molecule, MOLECULES.Molecule):
            # If this isn't a design molecule, it can be safely skipped
            if not molecule.design:
                continue
            # Otherwise calculate the dihedral angles
            molecule.calculate_dihedrals()
        # Loop through the Residues
        for residue in molecule:
            # Make sure the permissions attribute of the residues are properly
            # set
            residue.permission = "ISOMER"
            # Get the rotamers for this Residue
            rotamers = gather_rotamers(residue, experiment, "closest")
            # If there are no rotamers, skip this Residue
            if rotamers == None:
                continue
            # Identify the rotamer with the smallest side chain RMSD with the
            # Residue's current side chain
            best = None
            bestRmsd = 50000
            for rotamer in rotamers:
                # Calculate its RMSD
                rmsd = 0
                count = 0
                # Loop through the side chain Atoms
                for atom in rotamer["side chain"]:
                    if atom.name in residue:
                        rmsd += MOLECULES.calculate_distance(atom, \
                                          residue[atom.name], True)
                        count += 1.0
                # If any Atoms were used to calculate the RMSD
                if count > 0:
                    rmsd = math.sqrt(rmsd/count)
                    # If this is the best rotamer identified so far
                    if rmsd < bestRmsd:
                        bestRmsd = rmsd
                        best = rotamer
            # If no best rotamer could be found, raise an error
            if best == None:
                text = "The Closest_Rotamer function failed to identify a "
                text += "rotamer for this Residue:\n" + str(residue)
                raise RotamerError(text)
            # Patch the best rotamer into the Residue
            rotamer_patcher(residue, best)
    # Nothing needs to be returned, because everything is handled by pointers

def Glycine(residues, experiment = None):
    """Mutate appropriate Residues to Glycine before a Backbone Perturbation."""
    # Validate the residues
    residues = validate_residues(residues)
    # Loop through the Molecules
    for molecule in residues:
        # If it is a Molecule, calculate its dihedral angles
        if isinstance(molecule, MOLECULES.Molecule):
            molecule.calculate_dihedrals()
        # Loop through the Residues
        for residue in molecule:
            # Get rotamers for the Residue
            rotamers = gather_rotamers(residue, experiment, "glycine")
            # If there are no rotamers for this Residue, skip it
            if rotamers == None:
                continue
            # There is only ever a single glycine Residue, so just use the first
            rotamer_patcher(residue, rotamers[0])
    # Thanks to pointers, nothing needs to be returned

def prep_optimal_rotamers(residues, experiment = None, probability = True):
    """Do all of the preparations for an Optimal Rotamer Selection."""
    # Validate the residues
    residues = validate_residues(residues)
    # Collect the rotamers for the system. As this is done, parameterize the
    # Residues and the rotamers for the non-bonded energy calculations
    rotamers = []
    probabilities = []
    # Loop through the Molecules
    for molecule in residues:
        # If it is a Molecule, calculate its dihedral angles
        if isinstance(molecule, MOLECULES.Molecule):
            molecule.calculate_dihedrals()
        # Loop through the Residues
        for i in range(len(molecule)):
            # Get the residue
            residue = molecule[i]
            # Determine if this is the N or C terminus of a Molecule
            N = False
            C = False
            if isinstance(molecule, MOLECULES.Molecule):
                if i == 0:
                    N = True
                if i == len(molecule) - 1:
                    C = True
            # Parameterize the Residue
            parameterize_Residue(residue, experiment, N, C)
            # Gather rotamers for the residue
            if probability:
                rots, probs = gather_rotamers(residue, experiment, \
                                              "optimization", probability)
            else:
                rots = gather_rotamers(residue, experiment, "optimization", \
                                       probability)
            # If there are no rotamers for the residue, continue the for loop
            if rots == None:
                continue
            # Parameterize each rotamer
            for rot in rots:
                parameterize_Residue(rot, experiment, N, C)
            # Store the rotamers and probabilities
            rotamers.extend(rots)
            if probability:
                probabilities.extend(probs)
    # Calculate the rotamer-constant energies
    RCE = calculate_rc_energies(residues, rotamers, probability, probabilities)
    # Sort the rotamers so that the problem is feasible in a computationally
    # reasonable period of time
    rotamers, RCE = sort_rotamers(rotamers, RCE, experiment)
    # Calculate the energies of these rotamers with each other
    RRE = calculate_rr_energies(rotamers)
    # Return the rotamers and energies so that an Optimal Rotamer Selection can
    # be run
    return rotamers, RCE, RRE

def Optimal_Rotamers(residues, experiment = None, probability = True):
    """Select an optimal combination of rotamers for a particular system."""
    # Do the preparations for this calculation
    rotamers, RCE, RRE =prep_optimal_rotamers(residues, experiment, probability)
    if experiment != None and experiment["Type"] in ['OptMAVEn']:
        # Use CPLEX to select an optimal combination of rotamers
        # Build the humanization integer cuts
        objective, solution = CPLEX.optimal_rotamer_selector(residues, rotamers, \
                                                        RCE, RRE, experiment)
    else:    
        # Use GAMS to select an optimal combination of rotamers
        objective, solution = GAMS.optimal_rotamer_selector(residues, rotamers, \
                                                        RCE, RRE, experiment)
    # Loop through the Residues and patch in the selected rotamers
    for molecule in residues:
        for residue in molecule:
            # If this Residue received a rotamer, identify it
            if residue.number in solution:
                # Extract the appropriate rotamer
                rotamer = rotamers[residue.number][solution[residue.number]]
                # Patch in the rotamer
                rotamer_patcher(residue, rotamer)
    # Return the objective value and the solution
    return objective, solution
