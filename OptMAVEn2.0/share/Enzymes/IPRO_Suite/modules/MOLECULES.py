#!/usr/bin/env python

# The name of this module
__name__ = "IPRO Suite Molecules Module"
# The documentation string of the module
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains classes and functions for the storage, validation, movement,
modification and output of molecules and in the IPRO suite of programs."""

# Include needed PYTHON modules
import os
import sys
import math
import numpy
import string
import datetime
# Include contents from the IPRO suite standards module. This includes all
# 'default' variables, the 'supported' lists, the amino acids and backbone atoms
# dictionarites, the IPRO_Error class, and the screen_formatting function
from STANDARDS import *
# Include the I/O CHECK module, which will allow the MoleculeFile class to
# search for files
import IO_CHECK

# Define an error for problems in this module. 
class MoleculeError(IPRO_Error):
    """An error for problems in the MOLECULES module of the IPRO suite."""
    def __init__(self, error = ''):
        """The initialization of the MoleculeError class."""
        IPRO_Error.__init__(self, error)

# Create functions for calculating the properties of molecules
def calculate_distance(atom1, atom2, squared = False):
    """Calculate the distance between two Atoms."""
    # Try to do the calculation
    try:
        # Store the calculated distance in this variable
        distance = 0
        # Loop through the first Atom's coordinates
        for i in range(len(atom1)):
            distance += math.pow(atom1[i] - atom2[i], 2)
        # If desired, take the square root of this value
        if not squared:
            distance = math.sqrt(distance)
        return distance
    # If there was an error
    except (TypeError, ValueError, IndexError, AttributeError, IPRO_Error):
        # Create a Molecule Error
        text = "There was an error calculating the distance between these Atoms"
        text += ":\n" + str(atom1).strip() + "\n" + str(atom2).strip()
        raise MoleculeError(text)

def make_unit_vector(atom1, atom2 = None):
    """Create a unit vector between two Atoms."""
    # Use a try statement to catch errors
    try:
        # Do this differently depending on whether or not 2 Atoms have been
        # provided
        if not atom2:
            vector = list(atom1)
        else:
            vector = []
            for i in range(len(atom1)):
                vector.append(atom1[i] - atom2[i])
        # Calculate the magnitude of the vector
        magnitude = 0
        for value in vector:
            magnitude += math.pow(value, 2)
        magnitude = math.sqrt(magnitude)
        # Make a unit vector
        for i in range(len(vector)):
            vector[i] /= magnitude
        return vector
    # If there was an error
    except (TypeError, IndexError, ValueError, AttributeError, IPRO_Error):
        text = "There was an error creating a unit vector from these inputs:\n"
        text += str(atom1).strip() + "\n" + str(atom2).strip()
        raise MoleculeError(text)

def calculate_cross_product(atom1, atom2):
    """Calculate the cross product between 2 containers of 3 coordinates."""
    # Use a try statement to catch errors
    try:
        # Make sure both inputs have three coordinates
        if len(atom1) != 3 or len(atom2) != 3:
            raise TypeError
        # Calculate the value
        cross = [atom1[1]*atom2[2] - atom1[2]*atom2[1], \
                 atom1[2]*atom2[0] - atom1[0]*atom2[2], \
                 atom1[0]*atom2[1] - atom1[1]*atom2[0]]
        return cross
    # If there is an error
    except (TypeError, IndexError, ValueError, AttributeError, IPRO_Error):
        text = "There was an error calculating the cross product from these "
        text += "inputs:\n" + str(atom1).strip() + "\n" + str(atom2)
        raise MoleculeError(text)

def calculate_dot_product(atom1, atom2):
    """Calculate the dot product between two containers of coordinates."""
    # use a try statement
    try:
        # Make sure the containers have the same length
        if len(atom1) != len(atom2):
            raise TypeError
        # calculate the dot product
        value = 0
        for i in range(len(atom1)):
            value += atom1[i]*atom2[i]
        return value
    # If there was a problem
    except (TypeError, IndexError, ValueError, AttributeError, IPRO_Error):
        text = "There was an error calculating the dot product from these "
        text += "inputs:\n" + str(atom1).strip() + "\n" + str(atom2)
        raise MoleculeError(text)

def calculate_dihedral_angle(atom1, atom2, atom3, atom4):
    """Calculate the dihedral angle between 4 Atoms."""
    # All of the functions used by this function have error checks built into
    # them, so there is no need to do so here.
    # Make unit vectors between the Atoms (the h vector is SUPPOSED to be
    # between Atom4 and Atom3, not the other way around)
    f = make_unit_vector(atom1, atom2)
    g = make_unit_vector(atom2, atom3)
    h = make_unit_vector(atom4, atom3)
    # Calculate the cross products between the vectors, and ensure they are unit
    # vectors
    a = make_unit_vector(calculate_cross_product(f, g))
    b = make_unit_vector(calculate_cross_product(h, g))
    # Calculate the dot product of these two cross products
    angle = math.degrees(math.acos(calculate_dot_product(a, b)))
    # The sign of the angle MAY need to be adjusted
    if calculate_dot_product(g, calculate_cross_product(a, b)) > 0:
        angle = - angle
    return angle

def calculate_rmatrix(in1, in2):
    """Calculate a rotation matrix.
        
    This function has two methods of working, depending on the provided inputs:
    1) Rodriguez's Rotation Formula calculate a rotation matrix that rotates
    counter-clockwise by "in1" RADIANS around the "in2" vector.
    2) Singular Value Decomposition minimizes the RMSD between two sets of
    coordinates (or Atoms) when applied to the second set."""
    # Use a try statement to catch any errors
    try:
        # If the inputs match the requirements for method 1
        if isinstance(in1, float) and len(in2) == 3 and \
        isinstance(in2, (list, Atom)):
            # Calculate some numbers that are needed in the calculations
            c = math.cos(in1)
            s = math.sin(in1)
            v = 1 - c
            # Calculate the entries in the rotation matrix
            t1 = c + v*in2[0]*in2[0]
            t2 = -s*in2[2] + v*in2[0]*in2[1]
            t3 = s*in2[1] + v*in2[0]*in2[2]
            t4 = s*in2[2] + v*in2[1]*in2[0]
            t5 = c + v*in2[1]*in2[1]
            t6 = -s*in2[0] + v*in2[1]*in2[2]
            t7 = -s*in2[1] + v*in2[2]*in2[0]
            t8 = s*in2[0] + v*in2[2]*in2[1]
            t9 = c + v*in2[2]*in2[2]
            # And assemble them into a 3 x 3 list
            rmatrix = [[t1, t2, t3], [t4, t5, t6], [t7, t8, t9]]
            return rmatrix
        # Or if they match the inputs for the second method
        elif len(in1) == len(in2) and len(in1) >= 3 and len(in1[0]) == 3:
            # Extract the coordinates in a way that is useful in numpy
            list1 = []
            list2 = []
            for i in range(len(in1)):
                list1.append([in1[i][0], in1[i][1], in1[i][2]])
                list2.append([in2[i][0], in2[i][1], in2[i][2]])
            # Create a 3 x 3 matrix for SVD
            matrix = numpy.dot(numpy.transpose(list1), list2)
            # Do the SVD
            v, s, w = numpy.linalg.svd(matrix)
            # If appropriate, modify the v matrix
            if numpy.linalg.det(v)*numpy.linalg.det(w) < 0:
                for i in range(len(v)):
                    v[i][-1] = -v[i][-1]
            # Calculate the rotation matrix
            matrix = numpy.matrix(numpy.dot(v, w))
            # Convert this into a 3 x 3 list
            rmatrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            A = matrix.A1
            i = 0
            j = 0 
            for k in range(9):
                if j == 3:
                    i += 1
                    j = 0
                rmatrix[i][j] = A[k]
                j += 1
            return rmatrix
        # Otherwise raise an error
        else:
            raise TypeError
    # If there is a problem, raise a MoleculeError
    except (TypeError, IndexError, ValueError, AttributeError, IPRO_Error):
        text = "There was a problem calculating a rotation matrix using these "
        text += "inputs:\n" + str(in1).strip() + "\n" + str(in2).strip()
        raise MoleculeError(text)

def move(obj, coordinates, sign = '-'):
    """Move a class from the MOLECULES module.

    This function calls the private _move function of the classes. This is where
    all error checking is done to make sure the move will happen as expected."""
    # Do error checking, which is the entire purpose of this function. First,
    # make sure the obj has a _move function
    if "_move" not in dir(obj) or "dimensions" not in dir(obj):
        text = "This object does not have the necessary methods and attributes "
        text += "to support the move function:\n" + str(obj)
        raise MoleculeError(text)
    # Make sure the coordinates are in an acceptable format
    elif not isinstance(coordinates, (list, Atom)):
        text = "These coordinates are not acceptable for moving an object:\n"
        text += str(coordinates)
        raise MoleculeError(text)
    # Make sure there are an appropriate number of coordinates
    elif len(coordinates) != obj.dimensions:
        text = "The following list of coordinates is not compatible with "
        text += "moving the object below it:\n" + str(coordinates) + "\n"
        text += str(obj)
        raise MoleculeError(text)
    # check the input sign
    elif sign not in ['-', '+']:
        text = str(sign) + " is not a supported sign for the move function."
        raise MoleculeError(text)
    # Check each coordinate in the coordinates list
    for coor in coordinates:
        if not isinstance(coor, (float, int)):
            text = "The move function requires all coordinates used to move an "
            text += "object to be numbers:\n" + str(coordinates)
            raise MoleculeError(text)
    # Call the object's move function
    obj._move(coordinates, sign)

def rotate(obj, rmatrix):
    """Rotate a class from the MOLECULES module.

    This function calls the private _rotate function of the classes. This is
    where all error checking for those functions should be carried out."""
    # Make sure the obj can use the _rotate function and pass the upcoming
    # validations
    if "_rotate" not in dir(obj) or "dimensions" not in dir(obj):
        text = "This object does not have the necessary methods and attributes "
        text += "for the rotate function to work properly:\n" + str(obj)
        raise MoleculeError(text)
    # Evaluate the rmatrix for appropriateness
    error = False
    if not isinstance(rmatrix, list) or len(rmatrix) != obj.dimensions:
        error = True
    else:
        for row in rmatrix:
            if not isinstance(row, list) or len(row) != obj.dimensions:
                error = True
                break
            else:
                for coor in row:
                    if not isinstance(coor, (float, int)):
                        error = True
                        break
                if error:
                    break
    # If the rotation matrix is not acceptable
    if error:
        text = "An attempt was made to rotate this object:\n" + str(obj)
        text += "\nThis failed because the rotation matrix was not a "
        text += str(obj.dimensions) + " x " + str(obj.dimensions) + " array "
        text += "of numbers:\n" + str(rmatrix)
        raise MoleculeError(text)
    # Otherwise just rotate the object
    obj._rotate(rmatrix)

def summarize(obj, exclude = []):
    """Summarize the contents of a MOLECULES class object."""
    # Make sure the object has a _summarize function
    if "_summarize" not in dir(obj):
        text = "This object does not have a built-in summarize function:\n"
        text += str(obj)
        raise MoleculeError(text)
    # Or if exclude isn't a list of strings
    error = False
    if not isinstance(exclude, list):
        error = True
    else:
        for value in exclude:
            if not isinstance(value, str):
                error = True
                break
    # If there is a problem with the exclude list
    if error:
        text = "The exclude value provided to the summarize function must be a "
        text += "list of strings, and this is not:\n" + str(exclude)
        raise MoleculeError(text)
    # Get the summary
    return obj._summarize(exclude)

# Have a function that checks the values assigned to attributes in the classes
# to make sure they are acceptable. Start with checking attributes associated
# with the PDB file format
def check_PDB_attribute(type, name, value):
    """Check the values of PDB file format attributes."""
    # Type is the class, name is the attribute name, and value is its assigned
    # value.
    # Keep track of whether or not an error has occured
    error = False
    # An Atom's kind should be ATOM or HETATM
    if type == "Atom" and name == "kind":
        if value not in ['ATOM', 'HETATM']:
            error = True
    # An Atom's number and a Residue's number must be no more than 5 characters
    # (i.e. between -9999 and 99999)
    elif (type == "Atom" and name in ['number', 'residueNumber']) or (type == \
    "Residue" and name in ['number', 'lastAtom']) or (type == "Molecule" and \
    name in ['lastAtom', 'lastResidue']):
        if not isinstance(value, int) or not -9999 <= value <= 99999:
            error = True
    # An Atom's name should be a capitalized string of 1 to 4 characters with no
    # spaces (isupper and isalnum confirm that there is at least 1 character).
    # Additionally, at least one of the characters must be a letter (guaranteed
    # by isupper).
    elif type == "Atom" and name == "name":
        if not isinstance(value, str) or not value.isupper() or not \
        value.isalnum() or len(value) > 4:
            error = True
    # A Residue's kind must be a capitalized string of 3 to 4 characters without
    # spaces
    elif (type == "Atom" and name == "residueKind") or (type == "Residue" and \
    name == "kind"):
        if not isinstance(value, str) or not value.isupper() or not \
        value.isalnum() or not 3 <= len(value) <= 4:
            error = True
    # A Residue's name must be a capitalized string of 5 or fewer characters,
    # and only the last character is permitted to be a letter, although it does
    # not have to be
    elif (type == "Atom" and name == "residueName") or (type == "Residue" and \
    name == "name"):
        # Check that it is a capitalized string of 1 to 5 characters
        if not isinstance(value, str) or value.upper() != value or not 1 <=\
        len(value) <= 5:
            error = True
        # If it passes those checks, try and make an integer from the value
        else:
            try:
                v = int(value)
            except ValueError:
                # If the entire thing isn't an integer, make sure everything but
                # the last character is an integer
                try:
                    v = int(value[:-1])
                    # And make sure the last character is a letter
                    if not value[-1].isalpha():
                        raise ValueError
                # If none of those things worked there's a problem
                except ValueError:
                    error = True
    # A Molecule's name must be a single, capitalized character but it MAY be a
    # space or an integer
    elif (type in ['Atom', 'Residue'] and name == "moleculeName") or (type == \
    "Molecule" and name == "name"):
        if not isinstance(value, str) or len(value) != 1 or value not in \
        string.ascii_uppercase + string.digits + " ":
            error = True
    # An Atom's coordinate should be a floating point number between -999.999
    # and 9999.000
    elif type == "Atom" and name == "coordinate":
        # If the value is an integer, change it to a float
        if isinstance(value, int):
            value = float(value)
        if not isinstance(value, float) or not -999.999 <= value <= 9999.000:
            error = True
    # The Atom's end attribute should be a string
    elif type == "Atom" and name == "end":
        if not isinstance(value, str):
            error = True
    # That is all of the attributes associated with the PDB file format, so if
    # the name is something else an error should be raised
    else:
        text = "The PDB file format does not support the " + str(name)
        text += " attribute of the " + str(type) + " class."
        raise MoleculeError(text)
    # If there's a problem, raise an error
    if error:
        text = "The PDB file format does not permit the " + str(name) 
        text += " attribute of the " + str(type) + " class to have a value of "
        text += str(value) + "."
        raise MoleculeError(text)

def check_CHARMM_attribute(type, name, value):
    """Check attributes associated with the CHARMM force field."""
    # Keep track of whether or not an error has occured
    error = False
    # These attributes should all be floating point numbers
    if type == "Atom" and name in ['vdwRadius', 'vdwEpsilon', 'vdw14R', \
    'vdw14E', 'charge', 'lkRadius', 'lkLambda', 'lkVolume', 'lkGibbs', \
    'ljPhi']:
        if value != None and not isinstance(value, float):
            error = True
    # No other attributes are supported
    else:
        text = "The CHARMM force field does not support the " + str(name)
        text += " attribute of the " + str(type) + " class."
        raise MoleculeError(text)
    # If there was a problem
    if error:
        text = "The CHARMM force field does not permit the " + str(name)
        text += " attribute of the " + str(type) + " class to have a value of "
        text += str(value) + "."
        raise MoleculeError(text)

def check_IPRO_attribute(type, name, value):
    """Check the values of IPRO decision making attributes."""
    # Keep track of whether or not an error has been detected
    error = False
    # A Residue's permission attribute specifies how it may be mutated when
    # interacting with the rotamer library
    if type == "Residue" and name == "permission":
        if value not in ['ROTAMER', 'ISOMER', 'FIXED']:
            error = True
    # A Residue's freedom attribute specifies how a Residue may move during
    # energy minimizations
    elif type == "Residue" and name == "freedom":
        if value not in ['FREE', 'RESTRAINED', 'FIXED']:
            error = True
    # A Residue's phi and psi dihedral angles must be floating point numbers
    # between -180 and 180
    elif type == "Residue" and name in ['phi', 'psi']:
        if value != None and (not isinstance(value, float) or not \
        -180 <= value <= 180):
            error = True
    # A Residue or Molecule's design attribute must be a boolean True or False
    # value
    elif type in ['Residue', 'Molecule'] and name == "design":
        if value not in [True, False]:
            error = True
    # A DesignGroup's number must be an integer >= 0
    elif type in ["DesignGroup", "OptzymeGroup"] and name == "number":
        if not isinstance(value, int) or value < 0:
            error = True
    # How a DesignGroup's interactions with its Target Molecules are being
    # modified
    elif type in ["DesignGroup", "OptzymeGroup"] and name == "objective":
        if value not in ['improve', 'maintain', 'reduce', 'eliminate']:
            error = True
    # The catalytic property of an OptZyme Group
    elif type == "OptzymeGroup" and name == "property":
        if value not in ["km", "kcat", "kcat/km"]:
            error = True 
    # The value of (RT)TSA/(RT)S used for kcat optimization in OptZyme        
    elif type == "OptzymeGroup" and name == "weight":
        if not isinstance(value, float) or value <= 0.0:
            error = True
    # Those are all of the IPRO decision making attributes
    else:
        text = "The IPRO suite of programs does not support the " + str(name)
        text += " attribute of the " + str(type) + " class."
        raise MoleculeError(text)
    # If there was an error
    if error:
        text = "The IPRO suite of programs does not permit the " + str(name)
        text += " attribute of the " + str(type) + " class to have a value of "
        text += str(value) + "."
        raise MoleculeError(text)

def check_fileFormat(value):
    """Make sure the file format is supported by the IPRO suite."""
    # supportedFormats, aminoAcids, and backboneAtoms are all in the STANDARDS
    # module
    if value not in supportedFormats:
        text = str(value) + " is not a file format supported by the IPRO suite "
        text += "of programs."
        raise MoleculeError(text)
    elif value not in aminoAcids:
        text = "The " + str(value) + " file format does not have an entry in "
        text += "the aminoAcids dictionary in the STANDARDS module."
        raise MoleculeError(text)
    elif value not in backboneAtoms:
        text = "The " + str(value) + " file format does not have an entry in "
        text += "the backboneAtoms dictionary in the STANDARDS module."
        raise MoleculeError(text)

def check_forceField(value):
    """Make sure the force field is supported by the IPRO suite of programs."""
    if value not in supportedFields:
        text = "The " + str(value) + " force field is not supported by the "
        text += "IPRO suite of programs."
        raise MoleculeError(text)

def is_hydrogen(atom):
    """Check to see if an Atom is a hydrogen or not."""
    # Have a different method based on the file format
    if atom.fileFormat == "PDB":
        # Go through the characters in the Atoms name. If the first letter is an
        # H, it is a hydrogen. Otherwise it is not
        for char in atom.name:
            # If the first character is a number, skip it
            if char.isdigit():
                continue
            # If it is H, return True
            elif char == "H":
                return True
            # Otherwise it is not a hydrogen
            else:
                return False
        # I don't know how the function could have reached this point, but if it
        # did say False
        return False
    # If the file format isn't supported, raise an error
    else:
        text = "The " + atom.fileFormat + " file format is not supported by the"
        text += " is_hydrogen function."
        raise MoleculeError(text)

class MolIter(object):
    """A basic class for iterating through the contents of MOLECULES classes"""
    def __init__(self, items):
        """The initialization of the MolIter class."""
        # Store the items
        self.items = items
        # Create an index for iteration
        self.index = -1
        # Store the number of items
        self.length = len(items)

    def next(self):
        """The next function used in iterating through a MolIter object."""
        # Increment the iteration index
        self.index += 1
        # Stop iterating when all items have already been returned
        if self.index == self.length:
            raise StopIteration
        # Otherwise return the proper item
        return self.items[self.index]

# Start creating the classes that store information about molecules. An Atom is
# a container of coordinates and each following class is a progressively more
# complicated container of Atoms
class Atom(object):
    """A class that stores information about a single atom in a molecule."""
    def __init__(self, line, number = None, resNum = None, resName = None, \
                 molName = None, forceField = defaultField, fileFormat = \
                 defaultFormat):
        """The initialization routine of the Atom class."""
        # Store the Atom's file format and force field
        self.fileFormat = fileFormat
        self.forceField = forceField
        # Have a different method of initialization for each supported file
        # format.
        if self.fileFormat == "PDB":
            self._PDB_init(line, number, resNum, resName, molName)
        else:
            text = "The Atom class does not support initialization using the "
            text += str(self.fileFormat) + " file format."
            raise MoleculeError(text)
        # Create the Atom's dimension attribute. This value is always controlled
        # by __setattr__
        self.dimensions = 0
        # Have additional initialization steps depending on the force field
        if self.forceField == "CHARMM":
            self._CHARMM_init()
        else:
            text = "The Atom class does not support initialization using the "
            text += str(self.forceField) + " force field."
            raise MoleculeError(text)

    def _PDB_init(self, line, number, resNum, resName, molName):
        """Store the contents of a PDB formatted line in the Atom."""
        # This function is private (it starts with a _), so it should never be
        # accessed from outside of the class and thus we can be confident that
        # the PDB file format is in fact being used.
        # Use a try statement to catch errors
        try:
            # The line must be a string of the right length that is not a
            # duplicated Atom
            if not isinstance(line, str) or not 60 <= len(line) <= 100 or \
            line[16] not in [' ', 'A']:
                raise ValueError
            # Get the Atom's kind
            self.kind = line[:6].strip()
            # Store the Atom's number
            if not number:
                self.number = int(line[6:11])
            else:
                self.number = number
            # The Atom's name
            self.name = line[12:16].strip()
            # The Atom's Residue's kind
            self.residueKind = line[17:21].strip()
            # The Atom's Molecule's name
            if not molName:
                self.moleculeName = line[21]
            else:
                self.moleculeName = molName
            # The Atom's Residue's name
            if not resName:
                self.residueName = line[22:28].strip()
            else:
                self.residueName = resName
            # The Atom's Residue's number
            if not resNum:
                # If the Residue's name is a string of all digits, make an int
                # of it
                if self.residueName.isdigit() or (self.residueName[0] == "-" \
                and self.residueName[1:].isdigit()):
                    self.residueNumber = int(self.residueName)
                # Otherwise only the last character shouldn't be an integer
                else:
                    self.residueNumber = int(self.residueName[:-1])
            else:
                self.residueNumber = resNum
            # The Atom's _coordinates list cannot directly have a value assigned
            # to it, so bypass the __setattr__ function of the Atom class to
            # create it.
            object.__setattr__(self, "_coordinates", [float(line[30:38]), \
                               float(line[38:46]), float(line[46:54])])
            # Store everything else in the end attribute
            self.end = line[54:].rstrip()
        # If there was an error
        except (TypeError, IndexError, ValueError, AttributeError):
            text = "The PDB file format does not accept the following as an "
            text += "input:\n" + str(line).strip()
            raise MoleculeError(text)

    def _CHARMM_init(self):
        """Initialize CHARMM associated attributes."""
        # This is a private function, so the Atom's file format doesn't need to
        # be checked
        for name in ['vdwRadius', 'vdwEpsilon', 'charge', 'lkRadius', \
        'lkLambda', 'lkVolume', 'lkGibbs', 'vdw14R', 'vdw14E', 'ljPhi']:
            self.__setattr__(name, None)

    def __setattr__(self, name, value):
        """Assign values to the Atom's attributes.

        This function exists to make sure that only acceptable values are
        assigned to the attributes."""
        # If the name is a method of the class, don't allow value assignment
        if name in ['__init__', '_PDB_init', '_CHARMM_init', '__setattr__', \
        '__getattribute__', '__getattr__', '__len__', '__iter__', \
        '__contains__', '__getitem__', '__setitem__', '__eq__', '__ne__', \
        '__format__', '_PDB_format', '_CHARMM_format', '_energy_format', \
        '__str__', '__repr__', 'duplicate', '_move', '_rotate', '_summarize']:
            text = name + " is a method of the Atom class and cannot be "
            text += "assigned a value."
            raise MoleculeError(text)
        # If the name is a private attribute that should not be directly
        # accessed
        elif name in ['_coordinates']:
            text = name + " is a private variable of the Atom class and cannot "
            text += "be assigned a value."
            raise MoleculeError(text)
        # If the Atom is being assigned to a Residue
        elif name == "residue":
            # Confirm that the value is a Residue
            if not isinstance(value, Residue):
                text = "The residue attribute of an Atom MUST be a Residue "
                text += "class object."
                raise MoleculeError(text)
            # Delete all of the attributes of the Atom that are references to
            # information the Residue should contain.
            for term in ['fileFormat', 'forceField', 'residueName', \
                         'residueKind', 'residueNumber', 'moleculeName']:
                if term in self.__dict__:
                    del self.__dict__[term]
            # Store the attribute
            object.__setattr__(self, name, value)
        # If the Atom belongs to a Residue and the attribute is an attribute of
        # a higher MOLECULES class
        elif "residue" in self.__dict__ and name in ['fileFormat','forceField',\
        'residueKind', 'residueNumber', 'residueName', 'moleculeName', \
        'groupNumber']:
            # Modify the attribute name, if necessary
            if name.startswith("residue"):
                name = name[7:].lower()
            # Store the attribute in the Residue
            self.residue.__setattr__(name, value)
        # If this is the file format
        elif name == 'fileFormat':
            # Validate the new value
            check_fileFormat(value)
            # If this is the first declaration of the file format, store it
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            # If the file format isn't changing, don't do anything
            elif value == self.fileFormat:
                pass
            # Otherwise have methods to convert between formats
            else:
                text = "The Atom class does not support converting the "
                text += str(self.fileFormat) + " file format to the "
                text += str(value) + " file format."
                raise MoleculeError(text)
        # Do the same for the force field
        elif name == "forceField":
            check_forceField(value)
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            elif value == self.forceField:
                pass
            else:
                text = "The Atom class does not support converting the "
                text += str(self.forceField) + " force field to the "
                text += str(value) + " force field."
                raise MoleculeError(text)
        # The dimensions attribute should always be the length of the
        # _coordinates list
        elif name == "dimensions":
            # Get the _coordinates list
            coors = object.__getattribute__(self, "_coordinates")
            object.__setattr__(self, name, len(coors))
        # PDB associated attributes
        elif self.fileFormat == "PDB" and name in ['kind', 'number', 'name', \
        'residueKind', 'moleculeName', 'residueName', 'residueNumber', 'end']:
            check_PDB_attribute("Atom", name, value)
            object.__setattr__(self, name, value)
        # CHARMM associated attributes
        elif self.forceField == "CHARMM" and name in ['vdwRadius','vdwEpsilon',\
        'vdw14R', 'vdw14E', 'charge', 'lkRadius', 'lkLambda', 'lkVolume', \
        'lkGibbs', 'ljPhi']:
            check_CHARMM_attribute("Atom", name, value)
            object.__setattr__(self, name, value)
        # Otherwise raise an error
        else:
            text = "The " + str(self.fileFormat) + " file format and the "
            text += str(self.forceField) + " force field do not support the "
            text += str(name) + " attribute in the Atom class."
            raise MoleculeError(text)

    def __getattribute__(self, name):
        """Access attribute values in the Atom class."""
        # Raise an error if there is an attempt to access the _coordinates
        if name == "_coordinates":
            error = "The _coordinates attribute of the Atom class is protected "
            error += "from direct access. Please use __getitem__."
            raise MoleculeError(error)
        # Try to access it using the normal method
        try:
            return object.__getattribute__(self, name)
        # And use __getattr__ if that fails
        except AttributeError:
            return self.__getattr__(name)

    def __getattr__(self, name):
        """Retrieve attributes from higher MOLECULE classes."""
        # Only do anything if the Atom belongs to a Residue
        if "residue" in self.__dict__ and name in ['fileFormat', 'forceField', \
        'residueName', 'residueKind', 'residueNumber', 'moleculeName', \
        'groupNumber']:
            # If the attribute name starts with "residue", adjust it
            if name.startswith("residue"):
                name = name[7:].lower()
            # Get the object from the Atom's Residue
            return self.residue.__getattribute__(name)
        # Otherwise raise an error
        else:
            error = "The Atom class does not have a " + name + " attribute."
            raise MoleculeError(error)
    
    def __len__(self):
        """The number of coordinates in the Atom."""
        # Update dimensions
        self.dimensions = 0 # Controlled by __setattr__
        return self.dimensions

    def __iter__(self):
        """The iterator of the Atom class."""
        # Get the coordinates of the Atom
        coors = object.__getattribute__(self, "_coordinates")
        # Make a new list of the coordinates
        items = []
        for i in range(len(coors)):
            items.append(coors[i])
        # Make a MolIter object
        obj = MolIter(items)
        return obj

    def __contains__(self, coor):
        """Determine if the specified coordinate is in the Atom."""
        # If it is a floating point number, compare with _coordinates
        if isinstance(coor, float):
            return coor in object.__getattribute__(self, "_coordinates")
        # If it is an integer
        elif isinstance(coor, int):
            return -self.dimensions <= coor < self.dimensions
        # Or a string
        elif isinstance(coor, str):
            return coor.upper() in ['X', 'Y', 'Z']
        # If it is a list, search it recursively
        elif isinstance(coor, list):
            if len(coor) == 0:
                return False
            for obj in coor:
                if not self.__contains__(obj):
                    return False
            return True
        # otherwise return false
        else:
            return False

    def __getitem__(self, coor):
        """Retrieve the specified coordinate from the Atom."""
        # Access the Atom's coordinates
        coordinates = object.__getattribute__(self, "_coordinates")
        # If the input is an integer in the proper range
        if isinstance(coor, int) and -self.dimensions <= coor < self.dimensions:
            return coordinates[coor]
        # If it is a string
        elif isinstance(coor, str) and coor.upper() in ['X', 'Y', 'Z']:
            return coordinates[['X', 'Y', 'Z'].index(coor.upper())]
        # Otherwise raise an error
        else:
            text = str(coor) + " is not a valid label to retrieve a coordinate "
            text += "from this Atom:\n" + self.__str__()
            raise MoleculeError(text)

    def __setitem__(self, coor, value):
        """Store the coordinate in the Atom."""
        # Check the value of the coordinate based on the file format
        if self.fileFormat == "PDB":
            check_PDB_attribute("Atom", "coordinate", value)
        # If the file format isn't supported
        else:
            text = "The __setitem__ function of the Atom class does not support"
            text += " the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)
        # If coor is an integer
        if isinstance(coor, int) and -self.dimensions <= coor < self.dimensions:
            object.__getattribute__(self, "_coordinates")[coor] = value
        # Or a string
        elif isinstance(coor, str) and coor.upper() in ['X', 'Y', 'Z']:
            object.__getattribute__(self, "_coordinates")[\
                                    ['X', 'Y', 'Z'].index(coor.upper())] = value
        # Raise an error
        else:
            text = str(coor) + " is not a valid label to store a coordinate in "
            text += "this Atom:\n" + self.__str__()
            raise MoleculeError(text)

    def __eq__(self, other):
        """The equality check of the Atom class.

        This function compares the Atoms' labelling information, NOT their
        coordinates."""
        # Confirm that the other is an Atom with the same file format and
        # dimensions and force field
        if not isinstance(other, Atom) or other.fileFormat != self.fileFormat \
        or other.dimensions != self.dimensions or other.forceField != \
        self.forceField:
            return False
        # Have a different check for each file format
        elif self.fileFormat == "PDB":
            return self.name == other.name and self.residueName == \
            other.residueName and self.residueKind == other.residueKind and \
            self.moleculeName == other.moleculeName
        # Raise an error if the file format is not supported
        else:
            text = "The __eq__ function of the Atom class does not support the "
            text += str(self.fileFormat) + " file format."
            raise MoleculeError(text)

    def __ne__(self, other):
        """The inequality check of the Atom class."""
        return not self.__eq__(other)

    def __format__(self, how = ''):
        """Create a formatted string of text describing the Atom."""
        # If how to format the Atom hasn't been specified, use the Atom's file
        # format as default
        if not isinstance(how, str) or how == '':
            how = self.fileFormat
        # Have different functions for each supported format
        if how == "PDB":
            return self._PDB_format()
        elif how == "CHARMM":
            return self._CHARMM_format()
        elif isinstance(how, str) and how.startswith("energy"):
            return self._energy_format(how)
        else:
            text = "The __format__ function of the Atom class does not support "
            text += "the " + str(how) + " file format."
            raise MoleculeError(text)

    def _PDB_format(self):
        """Create a string of PDB formatted text."""
        # First, confirm that the Atom uses the PDB file format
        if self.fileFormat != "PDB":
            text = "An Atom using the " + str(self.fileFormat) + " file format "
            text += "cannot be used to create a PDB formatted output."
            raise MoleculeError(text)
        # Put this in a try statement in case there is an error. __setattr__
        # should make that impossible, but be cautious anyway
        try:
            # Include the Atom's kind and number
            text = self.kind.ljust(6) + str(self.number).rjust(5)
            # The Atom's name
            if len(self.name) > 3:
                text += " " + self.name.ljust(5)
            else:
                text += "  " + self.name.ljust(4)
            # Use HIS for all histidines
            if self.residueKind in ['HSD']:
                text += "HIS "
            else:
                text += self.residueKind.ljust(4)
            # The Molecule's name
            text += self.moleculeName
            # The Atom's Residue's name, with proper spacing
            if self.residueName[-1].isdigit():
                n = self.residueName
                s = ''
            else:
                n = self.residueName[:-1]
                s = self.residueName[-1]
            if len(n) == 5:
                text += n + s.ljust(3)
            else:
                text += n.rjust(4) + s.ljust(4)
            # The Atom's coordinates
            for coor in self:
                text += format(coor, '.3f').rjust(8)
            # The end attribute of the atom
            text += self.end + "\n"
            return text
        # If there was an error
        except (TypeError, IndexError, ValueError, AttributeError):
            text = "There was an error formatting an Atom, which should be "
            text += "impossible. Good luck identifying the error."
            raise MoleculeError(text)

    def _CHARMM_format(self):
        """Create a string of text formatted for use in CHARMM."""
        # Have a different method for each supported FILE FORMAT
        if self.fileFormat == "PDB":
            # Use a try statement
            try:
                # This will be almost the same as the _PDB_format function, so
                # I'm only going to comment on differences. Starting with only
                # using ATOM for Atoms' kind attributes
                text = "ATOM  " + str(self.number).rjust(5)
                if len(self.name) > 3:
                    text += " " + self.name.ljust(5)
                else:
                    text += "  " + self.name.ljust(4)
                # Use HSD for all histidines
                if self.residueKind in ['HIS']:
                    text += "HSD "
                else:
                    text += self.residueKind.ljust(4)
                text += self.moleculeName
                # Use the Residue's number, not its name
                n = str(self.residueNumber)
                if len(n) >= 4:
                    text += n.rjust(5) + ''.ljust(3)
                else:
                    text += n.rjust(4) + ''.ljust(4)
                for coor in self:
                    text += format(coor, '.3f').rjust(8)
                text += self.end + "\n"
                return text
            # If there was an error
            except (TypeError, IndexError, ValueError, AttributeError):
                text = "There was an error formatting an Atom, which should "
                text += "be impossible. Good luck finding the problem."
                raise MoleculeError(text)
        # If the file format is not supported
        else:
            text = "The Atom class does not support creating a formatted line "
            text += "of text for CHARMM when the " + str(self.fileFormat) 
            text += " file format is being used."
            raise MoleculeError(text)

    def _energy_format(self, how):
        """Create a string of text for non-bonded energy calculations."""
        # Store the text in this string
        text = ''
        # Have a different method for each file format
        if self.fileFormat == "PDB":
            # Extract the rotamer number from the how input
            try:
                rotNum = str(int(how.split('-')[1].strip()))
            except (IndexError, ValueError):
                rotNum = '1'
            # Now assemble the PDB labelling information - molecule name,
            # residue number, residue kind, rotamer number, atom name, and
            # coordinates
            text += "PDB "
            if self.moleculeName == ' ':
                text += "_ "
            else:
                text += self.moleculeName + " "
            text += str(self.residueNumber) + " " + self.residueKind + " "
            text += rotNum + " " + self.name
            for coor in self:
                text += " " + format(coor, '.3f')
        # If the file format is not supported
        else:
            text = "Formatting for non-bonded energy calculations does not "
            text += "support the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)
        # Have a different method for each force field
        if self.forceField == "CHARMM":
            text += " CHARMM"
            # Check each set of parameters
            if isinstance(self.vdwRadius, float) and \
            isinstance(self.vdwEpsilon, float) and \
            isinstance(self.vdw14R, float) and isinstance(self.vdw14E, float) \
            and isinstance(self.ljPhi, float):
                text += " VDW " + format(self.vdwRadius, '.4f')
                text += " " + format(self.vdwEpsilon, '.4f')
                text += " " + format(self.vdw14R, '.4f')
                text += " " + format(self.vdw14E, '.4f')
                text += " " + format(self.ljPhi, '.4f')
            if isinstance(self.charge, float):
                text += " ELEC " + format(self.charge, '.4f')
            if isinstance(self.lkRadius, float) and \
            isinstance(self.lkLambda, float) and isinstance(self.lkGibbs,float)\
            and isinstance(self.lkVolume, float):
                text += " LK " + format(self.lkRadius, '.4f')
                text += " " + format(self.lkLambda, '.4f')
                text += " " + format(self.lkVolume, '.4f')
                text += " " + format(self.lkGibbs, '.4f')
        # If the force field is not supported
        else:
            text = "Formatting for non-bonded energy calculations does not "
            text += "support the " + str(self.forceField) + " force field."
            raise MoleculeError(text)
        # Add an end line character to the text and return it
        text += "\n"
        return text

    def __str__(self):
        """Create a formatted string of text for the Atom."""
        # Just call the Atom's format function, using the Atom's file format
        return self.__format__(self.fileFormat)

    def __repr__(self):
        """Create a formatted string of text for the Atom."""
        # Just call the __str__ function
        return self.__str__()

    def duplicate(self):
        """Create a new Atom class object with identical attributes."""
        # Get the formatted text for the Atom
        text = self.__str__()
        # Create a new atom
        return Atom(text, self.number, self.residueNumber, self.residueName, \
                    self.moleculeName, self.forceField, self.fileFormat)

    def _move(self, coordinates, sign):
        """A function to move the Atom's coordinates.

        This function should only be called from the general move function at
        the beginning of the MOLECULES module for error-checking reasons."""
        # Wrap everything in a try statement in the event of an error
        try:
            if sign == '-':
                for i in range(self.dimensions):
                    object.__getattribute__(self, "_coordinates")[i] -= \
                                                            coordinates[i]
            elif sign == "+":
                for i in range(self.dimensions):
                    object.__getattribute__(self, "_coordinates")[i] += \
                                                            coordinates[i]
            else:
                raise TypeError
            # Validate that the new coordinates are all acceptable
            if self.fileFormat == "PDB":
                for coor in self:
                    check_PDB_attribute("Atom", "coordinate", coor)
            else:
                text = "The _move function of the Atom class does not support "
                text += "the " + str(self.fileFormat) + " file format."
                raise MoleculeError(text)
        # If there was an error
        except (TypeError, IndexError, ValueError, AttributeError):
            text = "There was an error moving an Atom. Use the general move "
            text += "function at the start of the MOLECULES module to find out "
            text += "why instead of the private function of the Atom class."
            raise MoleculeError(text)

    def _rotate(self, rmatrix):
        """Rotate the Atom."""
        # Put everything in a try statement
        try:
            # Store the new coordinates of the Atom in this list
            new = []
            # Loop through the dimensions of the Atom
            for i in range(self.dimensions):
                coor = 0
                for j in range(self.dimensions):
                    coor += self[j]*rmatrix[i][j]
                new.append(coor)
            # Check the new coordinates
            if self.fileFormat == "PDB":
                for coor in new:
                    check_PDB_attribute("Atom", "coordinate", coor)
            else:
                text = "The _rotate function of the Atom class does not support"
                text += " the " + str(self.fileFormat) + " file format."
                raise MoleculeError(text)
            # Store the new coordinates in the Atom
            object.__setattr__(self, "_coordinates", new)
        # If there was an error
        except (TypeError, IndexError, ValueError, AttributeError):
            text = "There was an error moving an Atom. Use the general rotate "
            text += "function at the start of the MOLECULES module to find out "
            text += "why instead of the private function of the Atom class."
            raise MoleculeError(text)

    def _summarize(self, exclude = []):
        """'Summarize' an Atom"""
        return self.__str__()

def gather_Atoms(atomLines, forceField, fileFormat, otherLines = False):
    """Collect all Atoms that make up a Residue or Molecule."""
    # if atomLines is a string, split it into a list
    if isinstance(atomLines, str):
        atomLines = atomLines.split("\n")
    # If atomLines is not a list now, raise an error
    if not isinstance(atomLines, list):
        raise TypeError
    # Store the collected Atoms in this list
    atoms = []
    # Store the lines that aren't Atoms in this list
    others = []
    # Loop through the Atoms
    for obj in atomLines:
        # If it is an Atom, keep it
        if isinstance(obj, Atom):
            atoms.append(obj)
        # Otherwise try and make an atom out of it
        else:
            try:
                atoms.append(Atom(obj, None, None, None, None, forceField, \
                                  fileFormat))
            except MoleculeError:
                others.append(obj)
    # If there are no Atoms, raise an error
    if len(atoms) == 0:
        raise TypeError
    # Return the appropriate lists
    if otherLines:
        return atoms, others
    else:
        return atoms

class Residue(object):
    """A subunit of a Molecule (e.g. an amino acid) that contains Atoms."""
    def __init__(self, atomLines, atomNum = None, number = None, name = None, \
        molName = None, forceField = defaultField, fileFormat = defaultFormat):
        """The initialization function of the Residue class."""
        # Store the Residue's file format and force field
        self.fileFormat = fileFormat
        self.forceField = forceField
        # Try to gather Atoms
        try:
            atomList = gather_Atoms(atomLines, forceField, fileFormat)
        except TypeError:
            text = "The following input is not valid to initialize a Residue:\n"
            text += str(atomLines)
            raise MoleculeError(text)
        # If no name has been specified for the Residue, use the first Atom's
        if not name:
            self.name = atomList[0].residueName
        else:
            self.name = name
        # Same thing for the Residue's Molecule's name
        if not molName:
            self.moleculeName = atomList[0].moleculeName
        else:
            self.moleculeName = molName
        # Store the number of spatial coordinates the Residue's Atoms have
        self.dimensions = atomList[0].dimensions
        # Load the Residue's Atoms
        self.load(atomList, atomNum, number)
        # Initialize the IPRO decision making variables
        # currentKind is controlled by __setattr__
        self.currentKind = ''
        # Whether this Residue can ever be mutated
        self.design = False
        # To what amino acids the Residue may be mutated
        self.permittedKinds = []
        # Whether the Residue may CURRENTLY be mutated. Don't allow any
        # Residue's side chain to be changed initially.
        self.permission = "FIXED"
        # Whether the Residue may move during energy minimizations. By default
        # allow Residues to move.
        self.freedom = "FREE"
        # The phi and psi dihedral angles
        self.phi = None
        self.psi = None
        # The list of Atoms in the Residue that are NEVER allowed to move
        self.fixedAtoms = []

    def load(self, atomLines, atomNum = None, number = None):
        """Store the Atoms in the Residue after validating them."""
        # Gather the atoms
        try:
            atomList = gather_Atoms(atomLines, self.forceField, self.fileFormat)
        except TypeError:
            text = "The following input is not valid to load a Residue from:\n"
            text += str(atomLines)
            raise MoleculeError(text)
        # Make sure every Atom uses the Residue's file format, force field, and
        # dimensions attributes
        for term in ['fileFormat', 'forceField', 'dimensions']:
            for atom in atomList:
                # Because these attributes may very well not ACTUALLY be stored
                # in the Atom / Residue any more, get their values and compare
                # them.
                v1 = atom.__getattribute__(term)
                v2 = self.__getattribute__(term)
                if v1 != v2:
                    text = "Each Atom used in the load function of the Residue "
                    text += "class must use the Residue's " + term+" attribute."
                    raise MoleculeError(text)
        # Have a different method for each supported file format
        if self.fileFormat == "PDB":
            self._PDB_load(atomLines, atomNum, number)
        else:
            text = "The load function of the Residue class does not support "
            text += "the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)
        # Store the Residue's length (number of Atoms) and renumber it
        self.length = 0 # controlled by __setattr__
        self.renumber(atomNum, number)
        # Make sure each Atom in the Residue knows it belongs to this Residue
        for atom in self:
            atom.residue = self

    def _PDB_load(self, atomList, atomNum, number):
        """Load a Residue's Atoms when the PDB file format is being used."""
        # This is a private function, so don't check to see if PDB is correct
        # Start by validating that each Atom has the same molecule name, residue
        # name, and residue kind, and has a unique name
        names = []
        for atom in atomList:
            # Check the molecule names, residue kinds and residue names
            for term in ['moleculeName', 'residueKind', 'residueName']:
                # Get the values of the appropriate Atoms
                v1 = atom.__getattribute__(term)
                v2 = atomList[0].__getattribute__(term)
                # Compare them
                if v1 != v2:
                    text = "Each Atom used in the load function of the Residue "
                    text += "class must have the same " + term + " attribute."
                    raise MoleculeError(text)
            # Check that each Atom has a unique name
            if atom.name not in names:
                names.append(atom.name)
            else:
                text = "Each Atom used in the load function of the Residue "
                text += "class must have a unique name."
                raise MoleculeError(text)
        # Store the Residue's kind based on the information in the Atom list
        self.kind = atomList[0].residueKind
        # __setattr__ makes accessing the _atomOrder list and _atoms dictionary
        # difficult, as a Residue is a container of Atoms
        object.__setattr__(self, "_atoms", {})
        object.__setattr__(self, "_atomOrder", [])
        # Loop through the Atoms and store them
        for atom in atomList:
            # Store the Atom in the Residue
            object.__getattribute__(self, "_atomOrder").append(atom.name)
            object.__getattribute__(self, "_atoms")[atom.name] = atom

    def renumber(self, atomNum = None, number = None):
        """Number the Residue's Atoms sequentially."""
        # If no atom number has been specified, use the number of the first Atom
        # in the Residue
        if not isinstance(atomNum, int):
            atomNum = self[0].number
        # Refer to __setattr__ for why these next two steps renumber the Atoms
        atomNum += (self.length - 1)
        self.lastAtom = atomNum
        # If no residue number has been specified, get it from the first Atom
        if not isinstance(number, int):
            number = self[0].residueNumber
        self.number = number

    def __setattr__(self, name, value):
        """Control attribute value assignment in the Residue class."""
        # If the name is a method of the class
        if name in ['__init__', 'load', '_PDB_load', 'renumber', '__setattr__',\
        '__getattribute__', '__getattr__', '__len__', '__iter__', \
        '__contains__', '__getitem__', '__setitem__', '__eq__', '__ne__', \
        '__format__', '__str__', 'duplicate', '__repr__', '_move', '_rotate', \
        '_summarize']:
            text = name + " is a method of the Residue class and does not "
            text += "support value assignment."
            raise MoleculeError(text)
        # If the name is a private attribute of the class.
        elif name in ['_atoms', '_atomOrder']:
            text = name + " is a private attribute of the Residue class and "
            text += "does not support value assignment."
            raise MoleculeError(text)
        # If the Residue is being assigned to a Molecule
        elif name == "molecule":
            # Make sure the value is a Molecule class object
            if not isinstance(value, Molecule):
                text = "The molecule attribute of a Residue must be a Molecule "
                text += "class object."
                raise MoleculeError(text)
            # Delete the values of higher level MOLECULES class attributes
            for term in ['fileFormat', 'forceField', 'moleculeName', \
            'groupNumber']:
                if term in self.__dict__:
                    del self.__dict__[term]
            # Store the Residue's Molecule
            object.__setattr__(self, name, value)
        # If the Residue belongs to a Molecule and the attribute is a higher
        # level MOLECULES class attribute, access it through the Residue's
        # Molecule
        elif "molecule" in self.__dict__ and name in ['fileFormat', \
        'forceField', 'moleculeName', 'groupNumber']:
            # If it is the name of the Residue's Molecule, modify the name
            if name == "moleculeName":
                name = "name"
            self.molecule.__setattr__(name, value)
        # The file format of the Residue
        elif name == "fileFormat":
            # Validate the value
            check_fileFormat(value)
            # If this is the first time a file format is stored
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            # If the file format isn't changing
            elif value == self.fileFormat:
                pass
            # Have a different method for converting between each pair of file
            # formats
            else:
                text = "The Residue class does not support converting the "
                text += str(self.fileFormat) + " file format to the "
                text += str(value) + " file format."
                raise MoleculeError(text)
        # Do the same for the force field
        elif name == "forceField":
            check_forceField(value)
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            elif self.forceField == value:
                pass
            else:
                text = "The Residue class does not support converting the "
                text += str(self.forceField) + " force field to the "
                text += str(value) + " force field."
                raise MoleculeError(text)
        # Attributes that are calculated properties of the Residue
        elif name == "length":
            # Number of Atoms in the Residue
            N = len(object.__getattribute__(self, "_atomOrder"))
            object.__setattr__(self, name, N)
        elif name == "dimensions":
            # Spatial coordinates of the Residue's Atoms. Guaranteed to all be
            # the same
            if "_atoms" not in self.__dict__:
                object.__setattr__(self, name, value)
            else:
                object.__setattr__(self, name, self[0].dimensions)
        # If the attribute is a decision making variable for IPRO
        elif name in ['permission', 'freedom', 'design', 'phi', 'psi']:
            check_IPRO_attribute("Residue", name, value)
            object.__setattr__(self, name, value)
            # If this is the design attribute and there are fixed in place
            # Atoms, make sure none of them are backbone Atoms
            if name == "design" and value == True and "fixedAtoms" in \
            self.__dict__:
                for label in self.fixedAtoms:
                    # If there are non-backbone Atoms that may never move, raise
                    # an error
                    if label not in backboneAtoms[self.fileFormat]:
                        text = "Design Positions are not permitted to contain "
                        text += "side chain Atoms that may never move."
                        raise MoleculeError(text)
        # Several other decision making variables must be handled independently
        elif name == "currentKind":
            # Always set this to the current kind attribute of the Residue
            object.__setattr__(self, name, self.kind)
        elif name == "permittedKinds":
            # Confirm the value is a list
            if not isinstance(value, list):
                text = "A Residue's permittedKinds attribute must be a list."
                raise MoleculeError(text)
            # The design attribute is initialized before the permittedKinds
            # list, so this check should always be valid
            if self.design:
                # Make sure each entry is an amino acid
                for label in value:
                    if label not in aminoAcids[self.fileFormat]:
                        text = "A Residue's permittedKinds attribute must be a "
                        text += "list containing only the names of amino acids "
                        text += "that are supported by the "
                        text += str(self.fileFormat) + " file format."
                        raise MoleculeError(text)
            else:
                # Must be an empty list
                if len(value) > 0:
                    text = "Only design Residues are permitted to have entries "
                    text += "in their permittedKinds attribute list."
                    raise MoleculeError(text)
            # If no error was raised, store the list
            object.__setattr__(self, name, value)
        elif name == "fixedAtoms":
            # Should be a list of Atoms that are in this Residue. This is only
            # declared AFTER the initial Atoms have been loaded.
            if not isinstance(value, list):
                text = "A Residue's fixedAtoms attribute must be a list."
                raise MoleculeError(text)
            for label in value:
                if label not in self:
                    text = "A Residue's fixedAtoms list may ONLY contain the "
                    text += "names of Atoms that are IN the Residue."
                    raise MoleculeError(text)
            # If this is a Design Position, make sure that all of the Atoms are
            # backbone atoms
            if "design" in self.__dict__ and self.design == True:
                for label in value:
                    if label not in backboneAtoms[self.fileFormat]:
                        text = "A Design Position may not contain side chain "
                        text += "Atoms that are never permitted to move."
                        raise MoleculeError(text)
            # Store the list
            object.__setattr__(self, name, value)
        # Have unique methods for storing values related to each file format.
        elif self.fileFormat == "PDB" and name in ['name', 'kind', 'number', \
        'moleculeName', 'lastAtom']:
            # Validate the value
            check_PDB_attribute("Residue", name, value)
            # Store the attribute
            object.__setattr__(self, name, value)
            # If this is the number of the last Atom in the Residue, adjust the
            # numbers of all Atoms accordingly
            if name == "lastAtom":
                # Adjust the number to that of the first Atom in the Residue
                value -= (self.length - 1)
                # Check its value
                check_PDB_attribute("Atom", "number", value)
                # If the first number is valid and the last number is valid, all
                # numbers in between must be. Assign the proper number to each
                # Atom
                for atom in self:
                    atom.number = value
                    value += 1
        # There are no other known attributes
        else:
            text = "The " + str(self.fileFormat) + " file format and the "
            text += str(self.forceField) + " force field do not support the "
            text += str(name) + " attribute in the Residue class."
            raise MoleculeError(text)

    def __getattribute__(self, name):
        """Get an attribute in the Residue."""
        # Control access to the _atoms and _atomOrder attributes
        if name in ["_atoms", '_atomOrder']:
            error = "The " + name + " attribute of the Residue class is "
            error += "protected from direct access. Please use __getitem__."
            raise MoleculeError(error)
        # Try to access an attribute from the __dict__ of the Residue
        try:
            return object.__getattribute__(self, name)
        # If that fails, try to use the __getattr__ function, which will
        # reference the Residue's Molecule
        except AttributeError:
            return self.__getattr__(name)

    def __getattr__(self, name):
        """Control access to higher level MOLECULES classes."""
        # Make sure the Residue has a Molecule
        if "molecule" in self.__dict__ and name in ['fileFormat', 'forceField',\
        'moleculeName', 'groupNumber']:
            # Modify moleculeName
            if name == "moleculeName":
                name = "name"
            # Get that attribute from the Molecule
            return self.molecule.__getattribute__(name)
        # Otherwise raise an error
        else:
            error = "The Residue class does not have a " + name + " attribute."
            raise MoleculeError(error)

    def __len__(self):
        """The number of Atoms in the Residue, (i.e. length attribute)."""
        # __setattr__ gets this correct, guaranteed
        self.length = 0
        return self.length

    def __iter__(self):
        """The iterator of the Residue class."""
        # Make a new list of atoms for the iteration
        new = []
        # Access the Atoms and their names directly
        atoms = object.__getattribute__(self, "_atoms")
        atomOrder = object.__getattribute__(self, "_atomOrder")
        # Loop through the names
        for atomName in atomOrder:
            new.append(atoms[atomName])
        # Make a MolIter object
        return MolIter(new)

    def __contains__(self, atom):
        """Determine if the Residue contains the specified Atom."""
        # If the atom input is a string
        if isinstance(atom, str):
            return atom in object.__getattribute__(self, "_atomOrder")
        # Or an integer
        elif isinstance(atom, int):
            return -self.length <= atom < self.length
        # Or an Atom class object
        elif isinstance(atom, Atom):
            # If any Atom in the Residue equals the provided Atom, return True
            for obj in self:
                if obj.__eq__(atom):
                    return True
            return False
        # If the input is a list
        elif isinstance(atom, list):
            # If the list is empy return False
            if len(atom) == 0:
                return False
            # Check each entry in the list
            for obj in atom:
                # If any of them aren't in the Residue, return False
                if not self.__contains__(obj):
                    return False
            # Otherwise return True
            return True
        # Otherwise just return False
        else:   
            return False

    def __getitem__(self, atom):
        """Retrieve the specified Atom from the Residue."""
        # Retrieve the Atoms and their order of appearance
        atoms = object.__getattribute__(self, "_atoms")
        atomOrder = object.__getattribute__(self, "_atomOrder")
        # If it is an integer in the proper range
        if isinstance(atom, int) and -self.length <= atom < self.length:
            return atoms[atomOrder[atom]]
        # Or the name of an atom
        elif atom in atomOrder:
            return atoms[atom]
        # If the input is a request for a list of Atoms
        elif isinstance(atom, str) and atom.lower() in ['all', 'backbone', \
        'side chain', 'sidechain']:
            # Store the atoms here
            ats = []
            # Have a unique option for each request
            if atom.lower() == "all":
                # If all Atoms are desired, store all of them
                for atom in self:
                    ats.append(atom)
            elif atom.lower() == "backbone":
                # Store everything that is not the side chain of an amino acid
                for atom in self:
                    if atom.name in backboneAtoms[self.fileFormat] or self.kind\
                    not in aminoAcids[self.fileFormat]:
                        ats.append(atom)
            else:
                # The third option is to store everything that is the side chain
                # of an amino acid
                for atom in self:
                    if atom.name not in backboneAtoms[self.fileFormat] and \
                    self.kind in aminoAcids[self.fileFormat]:
                        ats.append(atom)
            return ats
        # If it is a list that will be valid for any of these options
        elif isinstance(atom, list):
            ats = []
            for obj in atom:
                ats.append(self.__getitem__(obj))
            return ats
        # Otherwise raise an error
        else:
            text = "This is not an acceptable label:\n" + str(atom).strip()
            text += "\nTo retrieve an Atom from this Residue:\n" + str(self)
            raise MoleculeError(text)

    def __setitem__(self, atomName, atom):
        """Store the specified Atom in the Residue."""
        # Make sure it is an Atom
        if not isinstance(atom, Atom):
            text = "Only Atom class objects may be stored in a Residue."
            raise MoleculeError(text)
        # Confirm that the Atom's name is being used to store the Atom
        if atomName != atom.name:
            text = "An Atom's name must be used to store it in a Residue."
            raise MoleculeError(text)
        # Make sure the Atom matches the Residue's file format, force field, and
        # dimensions attributes
        for label in ['fileFormat', 'forceField', 'dimensions']:
            # Get the values of the attributes in the Atom and the Residue
            v1 = atom.__getattribute__(label)
            v2 = self.__getattribute__(label)
            # compare the values
            if v1 != v2:
                text = "All Atoms in a Residue must use the Residue's " + label
                text += " attribute."
                raise MoleculeError(text)
        # check the residueKind
        if atom.residueKind != self.kind:
            text = "All Atoms in a Residue must use the Residue's kind "
            text += "attribute as their residueKind attribute."
            raise MoleculeError(text)
        # Make sure the Atom knows it belongs to this Residue
        atom.residue = self
        # If this Atom is replacing a current Atom
        if atom.name in self:
            atom.number = self[atom.name].number
            object.__getattribute__(self, "_atoms")[atom.name] = atom
        # Otherwise things will have to be renumbered
        else:
            # Put this Atom last, with the correct number
            atom.number = self.lastAtom + 1
            object.__getattribute__(self, "_atomOrder").append(atom.name)
            object.__getattribute__(self, "_atoms")[atom.name] = atom
            # The validity of the number has been checked by the Atom, so
            # just store it
            object.__setattr__(self, "lastAtom", atom.number)

    def __format__(self, how = ''):
        """Create a formatted string of text describing the Residue's Atoms."""
        # Store the text in this string
        text = ''
        # Determine which Atoms to use in the formatting
        if isinstance(how, str) and '-' in how:
            # Get the atoms specification
            atoms = how.split("-")[0].strip()
            how = "energy - " + how.split("-")[1]
        # If there is no specification, use all Atoms
        else:
            atoms = "all"
        # Retrieve the specified Atoms
        for atom in self[atoms]:
            text += format(atom, how)
        return text

    def __str__(self):
        """Create a formatted string of text summarizing the Residue's Atoms."""
        return self.__format__(self.fileFormat)

    def __repr__(self):
        """Call the __str__ function of the Residue."""
        return self.__str__()

    def duplicate(self):
        """Create a new Residue with this one's attributes."""
        # Create the formatted text
        text = self.__str__()
        # Return a new Residue
        return Residue(text, self[0].number, self.number, self.name, \
                       self.moleculeName, self.forceField, self.fileFormat)

    def __eq__(self, other):
        """The equality comparison of the Residue class."""
        # Confirm that the other is a Residue with the same number of Atoms and
        # the same design status
        if not isinstance(other, Residue) or other.length != self.length or \
        other.design != self.design:
            return False
        # Because the Atoms of a Residue must use the Residue's attributes,
        # checking that the Atoms are equal is sufficient.
        for atom in self:
            # If any Atom doesn't match up, return False
            if atom.name not in other or atom.__ne__(other[atom.name]):
                return False
        # If all Atoms match up return True
        return True

    def __ne__(self, other):
        """The not equals check of the Residue class."""
        return not self.__eq__(other)

    def _move(self, coordinates, sign):
        """Move the Residue.

        This function should only be called by the general move function at the
        start of the MOLECULES module."""
        for atom in self:
            atom._move(coordinates, sign)

    def _rotate(self, rmatrix):
        """Rotate the Residue.

        This function should only be called by the general rotate function at 
        the start of the MOLECULES module."""
        for atom in self:
            atom._rotate(rmatrix)

    def _summarize(self, exclude = []):
        """List the Atoms in the Residue"""
        # Store the summary here
        summary = ''
        # If there are Atoms to leave out
        if len(exclude) > 0:
            summary += "Eligible "
        summary += "Atoms in Residue " + self.name + " in Molecule " + \
                   self.moleculeName + ":"
        # Get the names of the relevant Atoms
        use = []
        for atom in self:
            if atom.name not in exclude:
                use.append(atom.name)
        # List items is a function in the STANDARDS module
        summary += list_items(use)
        return summary

class Molecule(object):
    """A group of covalently bonded Atoms (e.g. a protein)."""
    def __init__(self, atomLines, atomNum = None, resNum = None, name = None, \
        design = False, forceField = defaultField, fileFormat = defaultFormat):
        """The initialization function of the Molecule class."""
        # Store the file format and force field
        self.fileFormat = fileFormat
        self.forceField = forceField
        # Try to gather Atoms
        try:
            atomList = gather_Atoms(atomLines, forceField, fileFormat)
        except TypeError:
            text = "The following is not a valid input to initialize a "
            text += "Molecule class object:\n" + str(atomLines)
            raise MoleculeError(text)
        # Store the Molecule's name
        if not name:
            self.name = atomList[0].moleculeName
        else:
            self.name = name
        # Store the Molecule's design status
        self.design = design
        # Store the Molecule's dimensions attribute
        self.dimensions = atomList[0].dimensions
        # Load the Molecule's structure
        self.load(atomList, atomNum, resNum)

    def load(self, atomLines, atomNum = None, resNum = None):
        """Load a Molecule's structure from a group of Atoms."""
        # Gather the needed Atoms
        try:
            atomList = gather_Atoms(atomLines, self.forceField, self.fileFormat)
        except TypeError:
            text = "The following is not a valid input for loading a "
            text += "Molecule's structure:\n" + str(atomLines)
            raise MoleculeError(text)
        # There are unique methods for loading depending on file format
        if self.fileFormat == "PDB":
            self._PDB_load(atomList)
        else:
            text = "The load function of the Molecule class does not support "
            text += "the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)
        # Number the Molecule sequentially
        self.renumber(atomNum, resNum)

    def _PDB_load(self, atomList):
        """Load PDB formatted Atoms into a Molecule."""
        # This is a private function, so error checking is not used because it
        # should only ever be used in controlled circumstances within the class.
        # First, make sure that each Atom uses the Molecule's file format, force
        # field, and dimensions and that they are all in the same Molecule
        for atom in atomList:
            for label in ['fileFormat', 'forceField', 'dimensions']:
                # Get the values of the Attributes, then compare them
                v1 = atom.__getattribute__(label)
                v2 = self.__getattribute__(label)
                if v1 != v2:
                    text = "All Atoms in a Molecule must use the Molecule's "
                    text += label + " attribute."
                    raise MoleculeError(text)
            if atom.moleculeName != atomList[0].moleculeName:
                text = "All Atoms in a Molecule must have the same moleculeName"
                text += " attribute."
                raise MoleculeError(text)
        # Separate the Atoms by Residue name
        residues = []
        residue = []
        oldName = atomList[0].residueName
        for atom in atomList:
            if atom.residueName != oldName:
                residues.append(residue)
                residue = []
                oldName = atom.residueName
            residue.append(atom)
        residues.append(residue)
        # Make sure that each Residue has a unique name
        names = []
        for res in residues:
            if res[0].residueName not in names:
                names.append(res[0].residueName)
            else:
                text = "Each Residue in a Molecule must have a unique name."
                raise MoleculeError(text)
        # If the Molecule is being initialized
        if "_residues" not in self.__dict__:
            # _residues and _residueOrder have special controls in __setattr__
            object.__setattr__(self, "_residues", {})
            object.__setattr__(self, "_residueOrder", [])
            # Store the Residues
            for res in residues:
                # Make a Residue
                residue = Residue(res)
                # Tell the Residue it belongs to this Molecule
                residue.molecule = self
                # Store the Residue
                object.__getattribute__(self, \
                                        "_residueOrder").append(residue.name)
                object.__getattribute__(self, "_residues")[residue.name] = \
                                                           residue
            # Calculate the number of Residues in the Molecule. This property is
            # controlled by __setattr__ to always be correct
            self.length = 0
        # If the structure is being reloaded
        else:
            # Make sure there are the proper number of Residues
            if len(residues) != self.length:
                text = "The number of Residues in a Molecule must not be "
                text += "changed when re-loading its structure."
                raise MoleculeError(text)
            # Loop through the Residues
            for i in range(self.length):
                self[i].load(residues[i])

    def renumber(self, atomNum = None, resNum = None):
        """Number the Molecule's Atoms and Residues sequentially."""
        # Make sure the Molecule knows how many Atoms it contains. This value is
        # controlled by __setattr__ to be always correct
        self.atomLength = 0
        # If no atom number has been specified
        if not isinstance(atomNum, int):
            atomNum = self[0][0].number
        # Use __setattr__ to do the numbering
        self.lastAtom = atomNum + self.atomLength - 1
        # Do the same for the Residues
        if not isinstance(resNum, int):
            resNum = self[0].number
        self.lastResidue = resNum + self.length - 1

    def __setattr__(self, name, value):
        """Control attribute assignment in the Molecule class."""
        # If the name is a method of the class
        if name in ['__init__', 'load', '_PDB_load', 'renumber', '__setattr__',\
        '__getattribute__', '__getattr__', '__len__', '__iter__', \
        '__contains__', '__getitem__', '__setitem__', '__eq__', '__ne__', \
        '__format__', '__str__', '__repr__', 'duplicate', '_move', '_rotate', \
        '_summarize', 'center', 'generate_name', 'output', \
        'calculate_dihedrals', 'identify_design_positions', \
        'identify_fixedAtoms', 'identify_permittedKinds']:
            text = name + " is a method of the Molecule class and does not "
            text += "support value assignment."
            raise MoleculeError(text)
        # If it is a private attribute of the class
        elif name in ['_residues', '_residueOrder']:
            text = name + " is a private attribute of the Molecule class and "
            text += "does not support value assignment."
            raise MoleculeError(text)
        # If it is assignment of the Molecule to a Design Group
        elif name == "group":
            # Make sure the value is a DesignGroup
            if not isinstance(value, DesignGroup):
                text = "The group attribute of a Molecule MUST be a DesignGroup"
                text += " class object."
                raise MoleculeError(text)
            # Delete the existing attributes that depend on the DesignGroup
            for term in ['fileFormat', 'forceField', 'groupNumber']:
                if term in self.__dict__:
                    del self.__dict__[term]
            # Store the Design Group the Molecule belongs to
            object.__setattr__(self, name, value)
        # If the Molecule belongs to a DesignGroup and one of the relevant
        # attributes is being modified
        elif "group" in self.__dict__ and name in ['fileFormat', 'forceField',\
        'groupNumber']:
            # modify groupNumber
            if name == "groupNumber":
                name = "number"
            # Set the value
            self.group.__setattr__(name, value)
        # The file format
        elif name == "fileFormat":
            # Validate the value
            check_fileFormat(value)
            # If this is the first time the file format is being stored
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            # If the value isn't being changed
            elif value == self.fileFormat:
                pass
            # Otherwise have unique methods for changing the Molecule's file
            # format
            else:
                text = "The Molecule class does not support converting the "
                text += str(self.fileFormat) + " file format to the "
                text += str(value) + " file format."
                raise MoleculeError(text)
        # Do the same for the force field
        elif name == "forceField":
            check_forceField(value)
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            elif value == self.forceField:
                pass
            else:
                text = "The Molecule class does not support converting the "
                text += str(self.forceField) + " force field to the "
                text += str(value) + " force field."
                raise MoleculeError(text)
        # If the attribute is a calculated property of the Molecule
        elif name == "length":
            # The number of Residues in the Molecule
            N = len(object.__getattribute__(self, "_residueOrder"))
            object.__setattr__(self, name, N)
        elif name == "atomLength":
            # Calculate and store the number of Atoms in the Molecule
            number = 0
            for residue in self:
                number += residue.length
            object.__setattr__(self, name, number)
        elif name == "dimensions":
            # All Residues are guaranteed to use the same value. If the Molecule
            # doesn't have Residues yet, just store the value
            if "_residues" not in self.__dict__:
                object.__setattr__(self, name, value)
            # If it does have Residues, use the value of the first one of them
            else:
                object.__setattr__(self, name, self[0].dimensions)
        # If the attribute is an IPRO decision making variable
        elif name == "design":
            check_IPRO_attribute("Molecule", name, value)
            object.__setattr__(self, name, value)
        # If the attribute is associated with the PDB file format
        elif self.fileFormat == "PDB" and name == "name":
            # Validate and store the value
            check_PDB_attribute("Molecule", name, value)
            object.__setattr__(self, name, value)
        # The numbering of the Molecule's Atoms
        elif self.fileFormat == "PDB" and name == "lastAtom":
            # Validate and store the number
            check_PDB_attribute("Molecule", name, value)
            object.__setattr__(self, name, value)
            # Adjust the number so it matches the first Atom's number
            value -= (self.atomLength - 1)
            # Check that new value
            check_PDB_attribute("Atom", "number", value)
            # Assign each Atom and Residue the proper numbers
            for residue in self:
                for atom in residue:
                    object.__setattr__(atom, "number", value)
                    value += 1
                object.__setattr__(residue, "lastAtom", value - 1)
        # The numbering of the Molecule's Residues
        elif self.fileFormat == "PDB" and name == "lastResidue":
            # Validate and store the number
            check_PDB_attribute("Molecule", name, value)
            object.__setattr__(self, name, value)
            # Modify it to match the first Residue's number
            value -= (self.length - 1)
            # Validate the new number
            check_PDB_attribute("Residue", "number", value)
            # Assign the proper value to the Residues (which the Atoms access)
            for residue in self:
                object.__setattr__(residue, "number", value)
                value += 1
        # That's all of the attributes, so raise an error if a problem is
        # encountered
        else:
            text = "The " + str(self.fileFormat) + " file format and the "
            text += str(self.forceField) + " force field do not support the "
            text += str(name) + " attribute of the Molecule class."
            raise MoleculeError(text)

    def __getattribute__(self, name):
        """Get an attribute of the Molecule."""
        # Control access to the _residues and _residueOrder attributes
        if name in ["_residues", "_residueOrder"]:
            error = "The " + name + " attribute of the Molecule class is "
            error += "protected from direct access. Please use __getitem__."
            raise MoleculeError(error)
        # Try the basic function from the object class
        try:
            return object.__getattribute__(self, name)
        # If that fails, use the __getattr__ function
        except AttributeError:
            return self.__getattr__(name)

    def __getattr__(self, name):
        """Control access to higher level MOLECULES attributes."""
        # Make sure the Molecule belongs to a DesignGroup
        if "group" in self.__dict__ and name in ['fileFormat', 'forceField', \
        'groupNumber']:
            # Modify groupNumber
            if name == "groupNumber":
                name = "number"
            # Get the value from the DesignGroup
            return object.__getattribute__(self.group, name)
        # Otherwise raise an error
        else:
            error = "The Molecule class does not have a " + name + " attribute."
            raise MoleculeError(error)

    def __len__(self):
        """The number of Residues in the Molecule."""
        self.length = 0 # __setattr__ gets this correct
        return self.length

    def __iter__(self):
        """The iterator of the Molecule class."""
        # Make a new list of Residues for iteration
        new = []
        # Access the Residues and their order directly
        residues = object.__getattribute__(self, "_residues")
        residueOrder = object.__getattribute__(self, "_residueOrder")
        # Store each Residue in new
        for residueName in residueOrder:
            new.append(residues[residueName])
        return MolIter(new)

    def __contains__(self, res):
        """Determine if the Molecule contains the specified Residue(s)."""
        # If the res input is an integer
        if isinstance(res, int):
            return -self.length <= res < self.length
        # Or a string of a Residue's name
        elif isinstance(res, str):
            return res in object.__getattribute__(self, "_residueOrder")
        # Or a Residue class object
        elif isinstance(res, Residue):
            for residue in self:
                if residue.__eq__(res):
                    return True
            return False
        # Or an Atom class object
        elif isinstance(res, Atom):
            for residue in self:
                for atom in residue:
                    if atom.__eq__(res):
                        return True
            return False
        # If it is a list of any of these things
        elif isinstance(res, list):
            # If the list is empty, return False
            if len(res) == 0:
                return False
            # Check each entry
            for obj in res:
                if not self.__contains__(obj):
                    return False
            return True
        # Otherwise say False
        else:
            return False

    def __getitem__(self, res):
        """Retrieve the specified Residue(s)."""
        # Access the Residues and the order they appear
        residues = object.__getattribute__(self, "_residues")
        residueOrder = object.__getattribute__(self, "_residueOrder")
        # If the residue input is an integer
        if isinstance(res, int) and -self.length <= res < self.length:
            return residues[residueOrder[res]]
        # If it is a string of a Residue's name
        elif isinstance(res, str) and res in residueOrder:
            return residues[res]
        # If it is a string requesting a subset of Residues
        elif isinstance(res, str) and res.lower() in ['all', 'ends', 'no ends']:
            # Store the residues in this list
            reses = []
            # Get all residues
            if res.lower() == "all":
                for residue in self:
                    reses.append(residue)
            # If only the first and last Residues are wanted
            elif res.lower() == "ends":
                reses.append(self.__getitem__(0))
                if self.length > 1:
                    reses.append(self.__getitem__(-1))
            # If all Residues but the N and C terminus are wanted
            else:
                for i in range(self.length):
                    if i == 0 or i == self.length - 1:
                        pass
                    else:
                        reses.append(self.__getitem__(i))
            return reses
        # If the input is a request for a subset of Residues
        elif isinstance(res, list):
            reses = []
            for obj in res:
                reses.append(self.__getitem__(obj))
            return reses
        # Otherwise raise an error
        else:
            text = "The following input is not valid:\n" + str(res).strip()
            text += "\nTo retrieve a Residue from this Molecule:\n"
            text += self.__str__()
            raise MoleculeError(text)

    def __setitem__(self, resName, res):
        """Store the specified Residue in the Molecule."""
        # Confirm it is a Residue
        if not isinstance(res, Residue):
            text = "Only Residue class objects may be stored in Molecules."
            raise MoleculeError(text)
        # Have a unique method for each supported file format
        if self.fileFormat == "PDB":
            # Make sure the Residue is being stored using its name
            if resName != res.name:
                text ="A Residue's name must be used to store it in a Molecule."
                raise MoleculeError(text)
            # Confirm that the Residue matches the Molecule's file format, force
            # field, and dimensions
            for label in ['fileFormat', 'forceField', 'dimensions']:
                # Get the values of the attributes and compare them.
                v1 = res.__getattribute__(label)
                v2 = self.__getattribute__(label)
                if v1 != v2:
                    text = "All Residues in a Molecule must use the Molecule's "
                    text += label + " attribute."
                    raise MoleculeError(text)
            # Get the current numbers of the first Atom and Residue in the
            # Molecule
            atomNum = self[0][0].number
            resNum = self[0].number
            # Make sure the Residue knows it belongs to this Molecule
            res.molecule = self
            # If the Residue is replacing a current Residue
            if resName in self:
                object.__getattribute__(self, "_residues")[resName] = res
            # Otherwise, it must be inserted before the first residue with a
            # higher integer portion of its name
            else:
                # Get the integer portion of the Residue's name
                try:
                    n1 = int(res.name)
                except ValueError:
                    n1 = int(res.name[:-1])
                # Identify where to insert the Residue
                index = self.length
                for residue in self:
                    # Get the Residue's name's integer portion
                    try:
                        n2 = int(residue.name)
                    except ValueError:
                        n2 = int(residue.name[:-1])
                    # Determine if it is higher than n1
                    if n2 > n1:
                        index = object.__getattribute__(self, \
                                       "_residueOrder").index(residue.name)
                        break
                # Insert the Residue's name in the molecule
                object.__getattribute__(self, "_residueOrder").insert(index, \
                                                                      resName)
                # Store the Residue
                object.__getattribute__(self, "_residues")[resName] = res
                # Update the Molecule's length
                self.length = 0
            # Number everything sequentially
            self.renumber(atomNum, resNum)
        # If the file format is not supported
        else:
            text = "The __setitem__ function of the Molecule class does not "
            text += "support the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)

    def __eq__(self, other):
        """The equality check of the Molecule class."""
        # Make sure the other is a Molecule with the same number of Residues and
        # the same design status
        if not isinstance(other, Molecule) or other.length != self.length or \
        other.design != self.design:
            return False
        # Compare every Residue, which covers the equality checking properly
        for res in self:
            if res.name not in other or res.__ne__(other[res.name]):
                return False
        return True

    def __ne__(self, other):
        """The inequality check of the Molecule class."""
        return not self.__eq__(other)

    def __format__(self, how = ''):
        """Create a formatted string of text for the Molecule's Atoms."""
        # Store the text in this string
        text = ''
        # Loop through the Molecule's residues
        for residue in self:
            # If non-bonded energy calculations are being done
            if how == "energy":
                # If the Residue is receiving a rotamer, only include its
                # backbone atoms
                if residue.permission in ['ROTAMER', 'ISOMER']:
                    text += format(residue, "backbone - 1")
                else:
                    text += format(residue, "all - 1")
            # Otherwise use the specified format
            else:
                text += format(residue, how)
        return text

    def _summarize(self, exclude = []):
        """Create a summary of the Molecule's Residues"""
        # Store the summary here
        summary = ''
        # If there are Residues to not include, say that
        if len(exclude) > 0:
            summary += "Eligible "
        summary += "Residues in Molecule " + self.name + ":"
        # Get the Residues that will be summarized
        residues = []
        for residue in self:
            if residue.name not in exclude:
                residues.append(residue)
        # If there are no Residues, say that
        if len(residues) == 0:
            summary += " NONE"
            return summary
        # Figure out the maximum length of a Residue's kind or name
        max = 0
        for residue in residues:
            if len(residue.name) > max:
                max = len(residue.name)
            if len(residue.kind) > max:
                max = len(residue.kind)
        # Increment this to account for a space between Residues
        max += 1
        # Calculate how many can fit in a line (SPACES is from STANDARDS)
        N = (SPACES-5) / max
        # Separate the Residues into groups of N
        groups = []
        group = []
        for residue in residues:
            if len(group) == N:
                groups.append(group)
                group = []
            group.append(residue)
        groups.append(group)
        # Summarize each Group
        for group in groups:
            summary += "\nNAME:"
            for residue in group:   
                summary += residue.name.rjust(max)
            summary += "\nKIND:"
            for residue in group:
                summary += residue.kind.rjust(max)
        return summary

    def __str__(self):
        """Create a summary of the contents of the Molecule."""
        return self._summarize()

    def __repr__(self):
        """Call the Molecule's __str__ function."""
        return self._summarize()

    def duplicate(self):
        """Create a new Molecule with this one's attributes."""
        # Create the formatted text
        text = self.__format__(self.fileFormat)
        # Make a new Molecule class object
        new = Molecule(text, self[0][0].number, self[0].number, self.name, \
                       self.design, self.forceField, self.fileFormat)
        # Match up the permission, freedom, fixedAtoms, and permittedKinds
        # attributes of the Residues
        for i in range(len(new)):
            new[i].design = self[i].design
            new[i].freedom = self[i].freedom
            new[i].permission = self[i].permission
            new[i].permittedKinds = list(self[i].permittedKinds)
            new[i].fixedAtoms = list(self[i].fixedAtoms)
        # Return the new molecule
        return new

    def _move(self, coordinates, sign):
        """Move the Molecule.

        This function should ONLY be called by the general move function at the
        beginning of the MOLECULES module for error-checking reasons."""
        for residue in self:
            residue._move(coordinates, sign)

    def _rotate(self, rmatrix):
        """Rotate the Molecule.

        This function should ONLY be called by the general rotate function at
        the start of the MOLECULES module for error-checkng reasons."""
        for residue in self:
            residue._rotate(rmatrix)

    def center(self, residues, backboneOnly = False):
        """Move the Molecule so that its average position is the origin."""
        # Keep track of the movements used to center the Molecule
        movements = []
        for i in range(self.dimensions):
            movements.append(0.0)
        # Keep track of how many Atoms are used in the centering
        atomCount = 0
        # Loop through the Residues specified by the residues input
        for residue in self[residues]:
            # If specified and appropriate, only use non-hydrogen atoms in
            # the backbones of amino acids
            if backboneOnly and residue.kind in aminoAcids[self.fileFormat]:
                # Have different lists depending on the file format
                if self.fileFormat == "PDB":
                    atoms = ['N', 'CA', 'C', 'O']
                else:
                    text = "The center function of the Molecule class does not "
                    text +="support the " + str(self.fileFormat)+" file format."
                    raise MoleculeError(text)
            elif isinstance(backboneOnly, list):
                atoms = backboneOnly
            else:
                atoms = object.__getattribute__(residue, "_atomOrder")
            # Loop through those Atoms
            for name in atoms:
                # Try to get the Atom
                try:
                    atom = residue[name]
                    atomCount += 1
                    # Include the Atom's coordinates in the centering
                    for i in range(atom.dimensions):
                        movements[i] += atom[i]
                except MoleculeError:
                    pass
        # Make an average of the movements list
        for i in range(self.dimensions):
            movements[i] /= atomCount
        # Move the Molecule
        self._move(movements, '-')
        # Return the movements
        return movements

    def generate_name(self, procedure = None, fileFormat=None):
        """Create the name of a file to store the Molecule's structure in."""
        # If no file format has been specified, use the Molecule's
        if not fileFormat:
            fileFormat = self.fileFormat
        # Have a unique method for each file format
        if fileFormat == "PDB":
            # Store the name in this string
            name = ''
            # Try to access a group number
            try:
                groupNum = self.groupNumber
            except MoleculeError:
                groupNum = None
            # If there is a group number
            if groupNum != None:
                # Make sure the group number (name) is acceptable
                if not isinstance(groupNum, str):
                    groupNum = str(groupNum)
                if ' ' in groupNum or '\t' in groupNum or '\n' in groupNum:
                    text = "This is not an acceptable name for a DesignGroup in"
                    text += " the output function of the Molecule class:\n"
                    text += groupNum
                    raise MoleculeError(text)
                name += "Group" + str(groupNum) + "_"
            # Include information about the Molecule
            name += "Molecule"
            if self.name != ' ':
                name += self.name
            else:
                name += "_"
            # If there is a procedure, include information about that
            if isinstance(procedure, str):
                # Remove whitespace from the procedure
                items = procedure.split()
                for item in items:
                    name += "_" + item
            name += ".pdb"
            return name
        # CHARMM needs the name to be lower cased.
        elif fileFormat == "CHARMM":
            name = self.generate_name(procedure, self.fileFormat)
            return name.lower()
        # If the file format is not supported
        else:
            text = "The generate_name function of the Molecule class does not "
            text += "support the " + str(fileFormat) + " file format."
            raise MoleculeError(text)

    def output(self, name = None, fileFormat = None, user = defaultUser):
        """Write the Molecule's structure to a file."""
        # If no file format has been specified, use the Molecule's default
        if not fileFormat:
            fileFormat = self.fileFormat
        # if no name has been given for the structure, create one
        if not isinstance(name, str) or ' ' in name or '\t' in name or \
        '\n' in name:
            name = self.generate_name(None, fileFormat)
        # Have a different method for each supported file format
        if fileFormat in ['PDB', 'CHARMM']:
            # Create a remark to head the file
            text = ''
            # Identify who is making this structure
            if user == defaultUser:
                try:
                    user = os.getlogin()
                except OSError:
                    pass
            text += "REMARK  Created by " + str(user) + " on "
            # Include the date and time
            a = datetime.datetime.now()
            text += str(a.month) + "/" + str(a.day) + "/" + str(a.year) + " at "
            text += str(a.hour).rjust(2, '0') + ":" + str(a.minute).rjust(2,'0')
            text += ":" + str(a.second).rjust(2, '0') + "\n"
            # Include the formatted text of the Molecule
            text += format(self, fileFormat)
            # End the Molecule's structure
            text += "END\n"
            # Write this to a file
            f = open(name, "w")
            f.write(text)
            f.close()
        # If the file format is not supported, raise an error
        else:
            text = "The output function of the Molecule class does not support "
            text += "the " + str(fileFormat) + " file format."
            raise MoleculeError(text)

    def calculate_dihedrals(self):
        """Calculate the dihedral angles of a Molecule's Residues."""
        # List the required Atoms by file format
        if self.fileFormat == "PDB":
            A1 = "C"
            A2 = "N"
            A3 = "CA"
            A4 = "C"
            A5 = "N"
        else:
            text = "The calculate_dihedrals function of the Molecule class "
            text += "does not support the " + str(self.fileFormat)
            text += " file format."
            raise MoleculeError(text)
        # Loop through the Molecule's Residues
        for i in range(self.length):
            # Get the current residue
            current = self[i]
            # If it is not an amino acid, skip it
            if current.kind not in aminoAcids[self.fileFormat]:
                current.phi = None
                current.psi = None
                continue
            # Extract the relevant backbone atoms
            try:
                atom2 = current[A2]
                atom3 = current[A3]
                atom4 = current[A4]
            # If those Atoms can't be gathered, no dihedral angles can be
            # calculated
            except MoleculeError:
                current.phi = None
                current.psi = None
                continue
            # Try to calculate the phi dihedral angle by retrieving the Atom
            # from the previous Residue
            try:
                # The first Residue doesn't have a phi dihedral angle
                if i == 0:
                    raise MoleculeError
                atom1 = self[i-1][A1]
                # Don't check to see if the previous is an amino acid, so
                # that non-canonical amino acids can be used to calculate
                # the dihedral angles of neighbors
                current.phi = calculate_dihedral_angle(atom1, atom2, atom3,\
                                                       atom4)
            except MoleculeError:
                current.phi = None
            # Do the same thing for psi
            try:
                if i == self.length - 1:
                    raise MoleculeError
                atom5 = self[i+1][A5]
                current.psi = calculate_dihedral_angle(atom2, atom3, atom4,\
                                                       atom5)
            except MoleculeError:
                current.psi = None

    def identify_design_positions(self, spots = {}):
        """Identify the Residues in the Molecule that may mutate."""
        # The default setting of the Residues is that they are not design
        # positions. Start by identifying if there are any design positions for
        # this Molecule in the spots dictionary
        if self.name in spots and len(spots[self.name]) > 0:
            any = True
        else:
            any = False
        # If there are any design positions and this is not a design Molecule,
        # raise an error
        if any and not self.design:
            text = "Only DESIGN Molecules may contain DESIGN Residues."
            raise MoleculeError(text)
        # Loop through the Molecule's Residues
        for residue in self:
            # if this Residue is a design position
            if any and residue.name in spots[self.name]:
                # If the Residue is not an amino acid, raise an error
                if residue.kind not in aminoAcids[self.fileFormat]:
                    text = "Only Residues that are amino acids are permitted "
                    text += "to be DESIGN Residues."
                    raise MoleculeError(text)
                # If the Residue has any side chain atoms that should always be
                # fixed in place, don't permit it to be a design position
                for atom in residue.fixedAtoms:
                    if atom not in backboneAtoms[self.fileFormat]:
                        text = "Residues that have side chain atoms that must "
                        text += "always be fixed in place are not permitted to "
                        text += "be DESIGN positions."
                        raise MoleculeError(text)
                # If those tests were passed, allow the position to be mutated
                residue.design = True
            else:
                residue.design = False

    def identify_fixedAtoms(self, spots = {}):
        """Identify the Atoms that are never permitted to be altered."""
        # Determine if there are any such Atoms in this Molecule
        if self.name in spots and len(spots[self.name]) > 0:
            any = True
        else:
            any = False
        # Loop through the Molecule's Residues
        for residue in self:
            # If this Residue has any fixed in place Atoms
            if any and residue.name in spots[self.name]:
                # Check the Atoms
                for atom in spots[self.name][residue.name]:
                    # The __setattr__ function of the Residue class confirms
                    # that the Atoms are in the Residue.
                    # If the Residue is a DESIGN residue and the atom is not a
                    # backbone atom
                    if residue.design and atom not in \
                    backboneAtoms[self.fileFormat]:
                        text = "DESIGN Residues may not have side chain Atoms "
                        text += "that are permanently fixed in place."
                        raise MoleculeError(text)
                # Store the Atoms if no error was raised
                residue.fixedAtoms = spots[self.name][residue.name]
            else:
                residue.fixedAtoms = []

    def identify_permittedKinds(self, permitted = {}):
        """Identify what amino acids DESIGN Residues may mutate to."""
        # Loop through the Molecule's Residues
        for residue in self:
            # The default is a blank list, so only consider DESIGN Residues
            if residue.design:
                # If there are listed permissions
                if self.name in permitted and residue.name in \
                permitted[self.name]:
                    residue.permittedKinds = permitted[self.name][residue.name]
                # If the attribute is a blank list, store the entire amino acids
                # dictionary. Remember, the default setting is a blank list, so
                # this is a valid check for EVERY residue
                if len(residue.permittedKinds) == 0:
                    residue.permittedKinds = aminoAcids[self.fileFormat]
            # If the Residue can't be mutated, store a blank list
            else:
                residue.permittedKinds = []

class MoleculeFile(object):
    """A class for parsing files containing Molecules."""
    def __init__(self, fileName, path = "./", forceField = defaultField, \
                 fileFormat = defaultFormat):
        """The initialization routine of the MoleculeFile class."""
        # Try to open the file
        try:
            f = IO_CHECK.for_file(fileName, path, True)
            lines = f.readlines()
            f.close()
        except IPRO_Error as error:
            raise MoleculeError(str(error))
        # Gather the Atoms and other lines from the file
        try:
            atomList, self.lines = gather_Atoms(lines, forceField, fileFormat,\
                                                True)
        except TypeError:
            text = fileName + " did not contain any Atoms."
            raise MoleculeError(text)
        # Store the File's file format and force field
        self.fileFormat = fileFormat
        self.forceField = forceField
        # Store the file's name
        self.name = fileName
        # Have a different initialization for each file format
        if self.fileFormat == "PDB":
            self._PDB_init(atomList)
        # If the file format is not supported
        else:
            text = "The __init__ function of the MoleculeFile class does not "
            text += "support the " + str(self.fileFormat) + " file format."
            raise MoleculeError(text)

    def _PDB_init(self, atomList):
        """Initialize a Molecule File that is using the PDB file format."""
        # Start by separating the Atoms by Molecule name
        names = []
        atoms = {}
        for atom in atomList:
            if atom.moleculeName not in names:
                names.append(atom.moleculeName)
                atoms[atom.moleculeName] = []
            atoms[atom.moleculeName].append(atom)
        # NMR structures will frequently have the same protein multiple times
        # using the same name in a file. We want to name each of them
        # independently. Here are the possible and permitted names
        letters = string.ascii_uppercase + string.digits
        # Keep track of the original names for each Molecule
        self.originals = {}
        # Separate the Atoms into Molecules
        molecules = {}
        self.moleculeNames = []
        # Loop through the different Molecule Names
        for I in names:
            # Store the current name that is being worked with
            name = I
            self.originals[I] = [name]
            self.moleculeNames.append(name)
            molecules[name] = []
            # Get the NUMBER of the first Atom with this Molecule Name. Because
            # no values were given to the Atoms this corresponds to the
            # numerical portion of the Atoms' Residues' names.
            old = atoms[I][0].residueNumber
            # Loop through all of the Atoms with this Molecule Name
            for atom in atoms[I]:
                # If this Atom's Residue comes before the last Atom's, then it
                # is the start of a new molecule
                if atom.residueNumber < old:
                    # Get a new name for the Molecule
                    i = 0
                    name = letters[i]
                    while name in self.moleculeNames or name in names:
                        i += 1
                        if i == len(letters):
                            text = "The _PDB_init function of the MoleculeFile "
                            text += "class has run out of optional Molecule "
                            text += "names. Please edit the function."
                            raise MoleculeError(text)
                        name = letters[i]
                    # Update the storage containers for the new name
                    self.originals[I].append(name)
                    self.moleculeNames.append(name)
                    molecules[name] = []
                # update the old number
                old = atom.residueNumber
                # Store the Atom
                molecules[name].append(atom)
        # Store the Atoms as Molecules
        self.molecules = {}
        for name in self.moleculeNames:
            self.molecules[name] = Molecule(molecules[name], None, None, \
                            name, False, self.forceField, self.fileFormat)

    def __len__(self):
        """The number of Molecules in the Molecule File."""
        return len(self.moleculeNames)

    def __iter__(self):
        """An iterator for the MoleculeFile class."""
        # Create a new list of the Molecules for iteration
        new = []
        for moleculeName in self.moleculeNames:
            new.append(self.molecules[moleculeName])
        return MolIter(new)

    def __contains__(self, label):
        """Determine if the specified Molecule is in the File."""
        return label in self.moleculeNames

    def __getitem__(self, label):
        """Retrieve the specified Molecule."""
        if isinstance(label, str) and label in self.moleculeNames:
            return self.molecules[label]
        elif isinstance(label, int) and -len(self) <= label < len(self):
            return self.molecules[self.moleculeNames[label]]
        else:
            text = str(label) + " is not a valid label to retrieve a Molecule "
            text += "from this MoleculeFile:\n" + str(self)
            raise MoleculeError(text)

    def _summarize(self, exclude = []):
        """Summarize the contents of a Molecule File"""
        # Store the summary here
        summary = ''
        # If there are Molecules to not include, indicate that
        if len(exclude) > 0:
            summary += "Eligible "
        summary += "Molecules in Molecule File " + self.name + ":"
        # Identify the Molecules that may be used
        molecules = []
        for molecule in self:
            if molecule.name not in exclude:
                molecules.append(molecule)
        # If there are no Molecules, say that
        if len(molecules) == 0:
            summary += " NONE"
            return summary
        # Include each Molecule
        for molecule in molecules:
            summary += "\nMolecule " + molecule.name + " contains " + \
                       str(len(molecule)) + " Residues"
            # Figure out the Molecules original name
            for I in self.originals:
                if molecule.name in self.originals[I]:
                    summary += "\n    Was originally Molecule " + I
                    break
        return summary

    def __str__(self):
        """Create a summary of the MoleculeFile's contents."""
        return self._summarize()

    def __format__(self, how = ''):
        """Call __str__"""
        return self.__str__()

    def __repr__(self):
        """Call __str__."""
        return self.__str__()

class DesignGroup(object):
    """A container of Molecules that are simultaneously present in a system."""
    def __init__(self, number, molecules, forceField = defaultField, \
                 fileFormat = defaultFormat):
        """The initialization routine of the DesignGroup class."""
        # Store the Group's file format, force field, and name
        self.fileFormat = fileFormat
        self.forceField = forceField
        self.number = number
        # A Design group should be initialized from a list of Molecules. Make
        # sure that's what is being done
        error = False
        if not isinstance(molecules, list) or len(molecules) == 0:
            error = True
        else:
            for obj in molecules:
                if not isinstance(obj, Molecule):
                    error = True
                    break
        if error:
            text = "A DesignGroup must be initialized using a list of Molecules"
            text += ", not:\n" + str(molecules)
            raise MoleculeError(text)
        # Store the Group's dimensions attribute from the first molecule
        self.dimensions = molecules[0].dimensions
        # The Group's moleculeOrder and molecules attributes are protected, so
        # create them in an unusual way
        object.__setattr__(self, "_molecules", {})
        object.__setattr__(self, "_moleculeOrder", [])
        # Check that each Molecule is compatible with the group, then store it.
        for molecule in molecules:
            for label in ['fileFormat', 'forceField', 'dimensions']:
                # Get and compare the values of the attributes
                v1 = molecule.__getattribute__(label)
                v2 = self.__getattribute__(label)
                if v1 != v2:
                    text = "Each Molecule stored in a DesignGroup must use the "
                    text += "Group's " + label + " attribute."
                    raise MoleculeError(text)
            # Make sure each has a unique name
            if molecule.name in object.__getattribute__(self, "_moleculeOrder"):
                text = "Each Molecule stored in a DesignGroup must have a "
                text += "unique name."
                raise MoleculeOrder(text)
            # Store the Molecule's name
            object.__getattribute__(self, \
                                    "_moleculeOrder").append(molecule.name)
            # Store a duplicate of the molecule so that DesignGroups definitely
            # contain different copies of Molecules. Make sure the duplicate
            # knows it belongs to this Design Group
            mol = molecule.duplicate()
            mol.group = self
            object.__getattribute__(self, "_molecules")[mol.name] = mol
        # Store the number of Molecule's in the Group. __setattr__ makes sure
        # the proper value is always used
        self.length = 0
        # Sort the Molecules so that the DESIGN Molecules come first
        # Start by getting the _molecules and _moleculeOrder attributes
        mols = object.__getattribute__(self, "_molecules")
        molOrder = object.__getattribute__(self, "_moleculeOrder")
        # Use a while loop to do the sort
        sorted = False
        while not sorted:
            # Set sorted to True
            sorted = True
            # Loop through the Molecule names
            for i in range(len(molOrder) - 1):
                # If this Molecule isn't a Design Molecule and the next is,
                # reverse them and indicate that everything isn't sorted yet
                if not mols[molOrder[i]].design and mols[molOrder[i+1]].design:
                    sorted = False
                    temp = molOrder[i]
                    molOrder[i] = molOrder[i+1]
                    molOrder[i+1] = temp
        # It should be fine, but store the moleculeOrder anyway
        object.__setattr__(self, "_moleculeOrder", molOrder)
        # Number all Molecules sequentially
        self.renumber("all")
        # The default is that the Design Group will try to improve binding
        self.objective = "improve"

    def renumber(self, which = "all"):
        """Number the DesignGroup's Molecules sequentially."""
        # Start the numberings at 1
        atomNum = 1
        resNum = 1
        # Loop through the Molecules
        for molecule in self[which]:
            molecule.renumber(atomNum, resNum)
            atomNum = molecule.lastAtom + 1
            resNum = molecule.lastResidue + 1

    def __setattr__(self, name, value):
        """Control attribute assignment in the DesignGroup class."""
        # If the attribute name is a method of the class
        if name in ['__init__', '_PDB_init', 'renumber', '__setattr__', \
        '__getattribute__', '__len__', '__iter__', '__contains__', \
        '__getitem__', '__setitem__', '__format__', '__str__', '__repr__', \
        'duplicate', '__eq__', '__ne__', '_move', '_rotate', '_summarize', \
        'output']:
            text = name + " is a method of the DesignGroup class and does not "
            text += "support value assignment."
            raise MoleculeError(text)
        # If the attribute is a private attribute of the class
        elif name in ['_molecules', '_moleculeOrder']:
            text = name + " is a private attribute of the DesignGroup class and"
            text += " does not support value assignment."
            raise MoleculeOrder(text)
        # File format must be checked
        elif name == "fileFormat":
            check_fileFormat(value)
            # If no file format has been stored before, store this one
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            # If it has, but the value isn't being changed
            elif value == self.__dict__[name]:
                pass
            # Otherwise raise an error
            else:
                text = "The DesignGroup class does not support converting the "
                text += str(self.fileFormat) + " file format to the "
                text += str(value) + " file format."
                raise MoleculeError(text)
        # Do the same for the force field
        elif name == "forceField":
            check_forceField(value)
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            elif value == self.__dict__[name]:
                pass
            else:
                text = "The __setattr__ function of the DesignGroup class does "
                text += "not support converting the " + str(self.forceField) 
                text += " force field to the " + str(value) + " force field."
                raise MoleculeError(text)
        # Attributes that are calculated properties of the Group
        elif name == "length":
            N = len(object.__getattribute__(self, "_moleculeOrder"))
            object.__setattr__(self, name, N)
        # Attributes that are associated with IPRO
        elif name in ["number", "objective"]:
            check_IPRO_attribute("DesignGroup", name, value)
            object.__setattr__(self, name, value)
        # Attributes that are semi-associated with a file format
        elif name == "dimensions":
            # If there are no molecules in the Design Group, use the provided
            # value
            if "_molecules" not in self.__dict__:
                object.__setattr__(self, name, value)
            # Otherwise use the value of the first Molecule in the group
            else:
                object.__setattr__(self, name, self[0].dimensions)
        # Otherwise raise an error
        else:
            text = "The " + str(self.fileFormat) + " file format and the "
            text += str(self.forceField) + " force field do not support the "
            text += "existence of the " + str(name) + " attribute of the "
            text += "DesignGroup class."
            raise MoleculeError(text)

    def __getattribute__(self, name):
        """Access attributes of the Design Group class."""
        # Control access to the _molecules and _moleculeOrder attributes
        if name in ["_molecules", "_moleculeOrder"]:
            error = "The " + name + " attribute of the DesignGroup class is "
            error += "protected from direct access. Please use __getitem__."
            raise MoleculeError(error)
        # Simply call the default get attribute function of the object class.
        return object.__getattribute__(self, name)

    def __len__(self):
        """The number of Molecules in the design group."""
        # Use __setattr__ to recalculate that value
        self.length = 0
        return self.length

    def __iter__(self):
        """The iterator of the DesignGroup class."""
        # Create a new list of the Molecules
        new = []
        # Access the molecules and their order directly
        molecules = object.__getattribute__(self, "_molecules")
        moleculeOrder = object.__getattribute__(self, "_moleculeOrder")
        # Store each Molecule in the list
        for moleculeName in moleculeOrder:
            new.append(molecules[moleculeName])
        return MolIter(new)

    def __contains__(self, mol):
        """Determine if the specified Molecule is in the DesignGroup."""
        # if the input is an integer
        if isinstance(mol, int):
            return -self.length <= mol < self.length
        # Or a string of a Molecule's name
        elif isinstance(mol, str):
            return mol in object.__getattribute__(self, "_moleculeOrder")
        # Or a Molecule
        elif isinstance(mol, Molecule):
            for molecule in self:
                if molecule.__eq__(mol):
                    return True
            return False
        # Or a Residue
        elif isinstance(mol, Residue):
            for molecule in self:
                for residue in molecule:
                    if residue.__eq__(mol):
                        return True
            return False
        # Or an Atom
        elif isinstance(mol, Atom):
            for molecule in self:
                for residue in molecule:
                    for atom in residue:
                        if atom.__eq__(mol):
                            return True
            return False
        # Or a list of any of those things
        elif isinstance(mol, list):
            # If it is empty, return false
            if len(mol) == 0:
                return False
            for obj in mol:
                if not self.__contains__(obj):
                    return False
            return True
        # Otherwise raise an error
        else:
            text = "This is not a valid object to check and see if it is "
            text += "contained in a DesignGroup:\n" + str(mol)
            raise MoleculeError(text)

    def __getitem__(self, mol):
        """Retrieve the specified Molecule(s) from the DesignGroup."""
        # Get the Molecules and the order they appear
        molecules = object.__getattribute__(self, "_molecules")
        moleculeOrder = object.__getattribute__(self, "_moleculeOrder")
        # if it is an integer input
        if isinstance(mol, int) and -self.length <= mol < self.length:
            return molecules[moleculeOrder[mol]]
        # Or a string of a Molecule's name
        elif isinstance(mol, str) and mol in moleculeOrder:
            return molecules[mol]
        # Or a string requesting a set of Molecules
        elif isinstance(mol, str) and mol.upper() in ['ALL', 'DESIGN','TARGET']:
            # Store the molecules in this list
            mols = []
            # Get all molecules
            if mol.upper() == "ALL":
                for molecule in self:
                    mols.append(molecule)
            # Get only the DESIGN molecules
            elif mol.upper() == "DESIGN":
                for molecule in self:
                    if molecule.design:
                        mols.append(molecule)
            # Get only the TARGET Molecules
            else:
                for molecule in self:
                    if not molecule.design:
                        mols.append(molecule)
            return mols
        # Or a list requesting a set of Molecules
        elif isinstance(mol, list):
            mols = []
            for obj in mol:
                mols.append(self.__getitem__(obj))
            return mols
        # Otherwise raise an error
        else:
            text = "The following is not an acceptable label to retrieve "
            text += "Molecule(s) from this DesignGroup:\n" + str(mol)
            raise MoleculeError(text)

    def __setitem__(self, molName, mol):
        """Store the specified Molecule in the DesignGroup."""
        # First, make sure the Molecule is a Molecule
        if not isinstance(mol, Molecule):
            text = "Only Molecule class objects may be stored in a DesignGroup."
            raise MoleculeError(text)
        # Make sure the Molecule's name is being used to store it
        if molName != mol.name:
            text = "A Molecule's name must be used to store it in a "
            text += "DesignGroup."
            raise MoleculeError(text)
        # confirm that the file format, force field, and dimensions attributes
        # match up
        for label in ['fileFormat', 'forceField', 'dimensions']:
            # Get and compare the values of the attributes
            v1 = mol.__getattribute__(label)
            v2 = self.__getattribute__(label)
            if v1 != v2:
                text = "All Molecules stored in a DesignGroup must use the "
                text += "Group's " + label + " attribute."
                raise MoleculeError(text)
        # Make a duplicate of the molecule and make sure it knows to belong to
        # this Design Group
        MOL = mol.duplicate()
        MOL.group = self
        # If the Molecule's name is not yet in the DesignGroup's list
        if MOL.name not in object.__getattribute__(self, "_moleculeOrder"):
            object.__getattribute__(self, "_moleculeOrder").append(MOL.name)
        # Store the Molecule in the Design Group
        object.__getattribute__(self, "_molecules")[MOL.name] = MOL
        # Update the length of the Design Group
        self.length = 0
        # Sort the Molecules so the DESIGN Molecules come first. Start by
        # getting the molecules and moleculeOrder attributes
        mols = object.__getattribute__(self, "_molecules")
        molOrder = object.__getattribute__(self, "_moleculeOrder")
        # Use a while loop to do the sorting
        sorted = False
        while not sorted:
            # Set sorted to true
            sorted = True
            # Loop through the Molecule names
            for i in range(len(molOrder) - 1):
                # If this molecule isn't a Design Molecule and the next is,
                # reverse them and say the sorting isn't done
                if not mols[molOrder[i]].design and mols[molOrder[i+1]].design:
                    sorted = False
                    temp = molOrder[i]
                    molOrder[i] = molOrder[i+1]
                    molOrder[i+1] = temp
        # make CERTAIN the Molecule Order is correct, even though pointers
        # dictate that it should be
        object.__setattr__(self, "_moleculeOrder", molOrder)
        # Number everything sequentially
        self.renumber("ALL")

    def __format__(self, how = ''):
        """Create a formatted string of text for the Group's Molecules."""
        # If how hasn't been specified, use the Group's file format
        if how == '':
            how = self.fileFormat
        # Store the formatted text in this string
        text = ''
        # Loop through the Molecules
        for molecule in self:
            # Add the formatting for that Molecule
            text += format(molecule, how)
            # If it is PDB related, add END to the end of the text
            if how == "PDB" or (how == "CHARMM" and self.fileFormat == "PDB"):
                text += "END\n"
        return text
    
    def _summarize(self, exclude = []):
        """Summarize the contents of a Design Group"""
        # Store the summary here
        summary = ''
        # If there are Molecules to not include, indicate that
        if len(exclude) > 0:
            summary += "Eligible "
        summary += "Molecules"
        # If this isn't Design Group 0, say what DG it is
        if self.number > 0:
            summary += " in Design Group " + str(self.number)
        summary += ":"
        # Figure out the Molecules to use
        molecules = []
        for molecule in self:
            if molecule.name not in exclude:
                molecules.append(molecule)
        # If there are no Molecules, say that
        if len(molecules) == 0:
            summary += " NONE"
            return summary
        # Include information for each Molecule
        for molecule in molecules:
            summary += "\nMolecule " + molecule.name + " contains " + \
                       str(len(molecule)) + " Residues"
        return summary

    def __str__(self):
        """Create a summary of the contents of the Design Group."""
        return self._summarize()

    def __repr__(self):
        """Call the __str__ function of the DesignGroup class."""
        return self._summarize()

    def duplicate(self):
        """Make a new Design Group that is a match for the current one."""
        # Make a new list of unique Molecules
        molecules = []
        for molecule in self:
            molecules.append(molecule.duplicate())
        # Return the new DesignGroup
        return DesignGroup(self.number, molecules, self.forceField, \
                           self.fileFormat)

    def __eq__(self, other):
        """The equality check of the Design Group class."""
        # Confirm that the other is a Design Group with an equal number of
        # Molecules
        if not isinstance(other, DesignGroup) or other.length != self.length:
            return False
        # Check each Molecule in the Design Groups
        for molecule in self:
            if molecule.name not in other or \
            molecule.__ne__(other[molecule.name]):
                return False
        return True

    def __ne__(self, other):
        """The not equals check of the Design Group class."""
        return not self.__eq__(other)

    def _move(self, coordinates, sign):
        """Move the DesignGroup."""
        for molecule in self:
            molecule._move(coordinates, sign)

    def _rotate(self, rmatrix):
        """Rotate the DesignGroup."""
        for molecule in self:
            molecule._rotate(rmatrix)

    def output(self, molecules = "ALL", path = "./", procedure = None, \
               fileFormat = None):
        """Output the specified Molecules in the specified folder."""
        # Number the specified Molecules sequentially
        self.renumber(molecules)
        # Make sure the path is usable
        if not isinstance(path, str):
            path = str(path)
        if ' ' in path or '\t' in path or '\n' in path:
            text = "This is not an acceptable path for outputting a "
            text += "DesignGroup:\n" + path
            raise MoleculeError(text)
        # If no file format has been specified, use the Group's
        if fileFormat == None:
            fileFormat = self.fileFormat
        # Output each Molecule
        for molecule in self[molecules]:
            # Generate a name for the molecule
            name = path + molecule.generate_name(procedure, fileFormat)
            # output the Molecule to that file
            molecule.output(name, fileFormat)

class OptzymeGroup(object):
    """A container of design groups that are involved in a single catalytic
    property optimization."""
    def __init__(self, number, assemblies, forceField = defaultField, \
                 fileFormat = defaultFormat):
        """The initialization routine of the OptzymeGroup class."""
        # Store the assembly's file format, force field, and name
        self.fileFormat = fileFormat 
        self.forceField = forceField
        self.number = number
        # If a single BindingAssembly is provided, convert it to a list
        if isinstance(assemblies, DesignGroup):
            assemblies = [assemblies]    
        # An Optzyme Group should be initialized from a list of design groups.
        # Ensure that is what is being done
        error = False
        if not isinstance(assemblies, list) or len(assemblies) == 0:
            error = True
        else:
            for obj in assemblies:
                if not isinstance(obj, DesignGroup):
                    error = True
                    break
        if error:
            text = "An OptzymeGroup must be initialized using a list of "
            text += "Design Groups, not:\n" + str(assemblies)
            raise MoleculeError(text)
        # Initialize the assemblies and assemblyOrder attributes
        object.__setattr__(self, "_assemblies", {})
        object.__setattr__(self, "_assemblyOrder", [])
        # Make sure the assemblies input is a list of BindingAssemblies
        if not isinstance(assemblies, list):
            text = "The initialization of an OptzymeGroup class object requires"
            text += " a list of BindingAssembly class objects as input. The "
            text += "provided input is:\n" + str(assemblies)
            raise MoleculeError
        # Make an initial copy of assemblies and assemblyOrder
        i = -1
        for assembly in assemblies:
            for obj in object.__getattribute__(self, "_assemblyOrder"):
                if obj == assembly:
                    text = "The same BindingAssembly was used twice in "
                    text += "the OptzymeGroup. This is not allowed."
                    raise MoleculeError(text)
            object.__getattribute__(self, "_assemblyOrder").append( \
            str(assembly.number))
            i += 1
            object.__getattribute__(self, "_assemblies")[i] = assembly
            object.__getattribute__(self, "_assemblies")[ \
            str(assembly.number)] = assembly
        # Reorder the BindingAssemblies in a consistent manner
        self.reorder()
        # Store the number of design groups within the Optzyme Group. __setattr__
        # makes sure the proper value is always used
        self.length = 0
        # Store the Optzyme Group's dimensions from the first molecule of the
        # first group
        self.dimensions = assemblies[0][0].dimensions
        # The default is that the Optzyme Group will try to improve binding
        self.objective = "improve"
        # By default, the Optzyme Group's relevant catalytic property is kcat/Km
        self.property = "kcat/km"
        # Store the value of (RT)TSA/(RT)S used for kcat calculations. The default
        # value is that of GUS (see documentation)
        self.weight = 0.03957

    def reorder(self):
        """Order the BindingAssemblies in a consistent manner."""
        # List binding assemblies by their objective. OptzymeGroups will use the
        # objective to separate the TSA binding assembly (objective: improve)
        # from the substrate binding assembly (objective: eliminate). Thus, the
        # TSA binding assembly will be listed before substrate binding
        # assemblies
        unsorted = []
        # Access the assemblies and their order directly
        assemblies = object.__getattribute__(self, "_assemblies")
        assemblyOrder = object.__getattribute__(self, "_assemblyOrder")
        # Store each BindingAssembly in the "unsorted" list
        for assemblyName in assemblyOrder:
            unsorted.append(assemblies[assemblyName]) 
        i = -1
        sorted = []
        for objective in ['improve', 'eliminate']:
            for assembly in unsorted:
                if not isinstance(assembly, DesignGroup):
                    text = "One of the entries in the assemblies list is not a "
                    text += "BindingAssembly class object:\n" + str(assembly)
                    raise MoleculeError
                if objective != assembly.objective:
                    continue
                # Store the BindingAssembly's number in the assemblyOrder
                sorted.append(str(assembly.number))
                # Store the BindingAssembly under the assemblies attribute using
                # either the BindingAssembly's number by integer
                i += 1
                object.__getattribute__(self, "_assemblies")[i] = assembly
                object.__getattribute__(self, \
                "_assemblies")[str(assembly.number)] = assembly
        object.__setattr__(self, "_assemblyOrder", sorted)

    def __setattr__(self, name, value):
        """Control attribute assignment in the OptzymeGroup class."""
        # If the attribute name is a method of the class
        if name in ['__init__', '__setattr__', '__getattribute__', '__len__', \
        '__iter__', '__contains__', '__getitem__', '__setitem__', '__format__',\
        '__str__', '__repr__', 'duplicate', '__eq__', '__ne__', '_move', \
        '_rotate', '_summarize', 'output', 'reorder']:
            text = name + " is a method of the OptzymeGroup class and does "
            text += "not support value assignment."
            raise MoleculeError(text)
        elif name in ['_assemblies', '_assemblyOrder']:
            text += name + " is a private attribute of the OptzymeGroup "
            text += "class and does not support value assignment."
            raise MoleculeError(text)
        elif name == "fileFormat":
            check_fileFormat(value)
            # If no file format has been stored before, store this one
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            # Raise an error if the fileFormat is being changed 
            elif value != self.__dict__[name]:
                text = "The OptzymeGroup class does not support converting "
                text += "the " + str(self.fileFormat) + " file format to the "
                text += str(value) + " file format."
                raise MoleculeError(text)
        # Do the same for the force field
        elif name == "forceField":
            check_forceField(value)
            if name not in self.__dict__:
                object.__setattr__(self, name, value)
            elif value != self.__dict__[name]:
                text = "The OptzymeGroup class does not support converting "
                text += "the " + str(self.forceField) + " force field to the "
                text += str(value) + " force field."
                raise MoleculeError(text)
        # Attributes that are calculated properties of the Optzyme Group
        elif name == "length":
            N = len(object.__getattribute__(self, "_assemblyOrder"))
            object.__setattr__(self, name, N)
        # Attributes that are associated with OptZyme
        elif name in ["property", "objective", "weight", "number"]:
            # Store the value of strings using lowercase characters
            if name not in ["weight", "number"]:
                value = value.lower()
            check_IPRO_attribute("OptzymeGroup", name, value)
            object.__setattr__(self, name, value)
        # Attributes that are semi-associated with a file format
        elif name == "dimensions":
            # If there are no Design Groups in the Optzyme Group, use the
            # provided value
            if "_assemblies" not in self.__dict__:
                object.__setattr__(self, name, value)
            # Otherwise, use the value of the first Design Group in the Optzyme
            # Group
            else:
                object.__setattr__(self, name, self[0].dimensions)
        # Otherwise, this name and/or value are invalid so raise an error
        else:            
            text = "The " + str(self.fileFormat) + " file format and the "
            text += str(self.forceField) + " force field do not support the "
            text += "existence of the " + str(name) + " attribute of the "
            text += "Optzyme Group class."
            raise MoleculeError(text) 

    def __getattribute__(self, name):
        """Access attributes of the Optzyme Group class."""
        # Control access to the _assemblies and _assemblyOrder attributes
        if name in ["_assemblies", "_assemblyOrder"]:
            error = "The " + name + " attribute of the OptzymeGroup class "
            error += "is protected from direct access. Please use __getitem__."
            raise MoleculeError(error)
        # Otherwise, simply call the default "get attribute" function of the 
        # object class.
        return object.__getattribute__(self, name)

    def __len__(self):
        """The number of Binding Assemblies in the Optzyme Group."""
        self.length = 0
        return self.length

    def __iter__(self):
        """The iterator of the OptzymeGroup class."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder()
        # Create a new list of the assemblies
        new = []
        # Access the assemblies and their order directly
        assemblies = object.__getattribute__(self, "_assemblies")
        assemblyOrder = object.__getattribute__(self, "_assemblyOrder")
        # Store each Binding Assembly in the list
        for assemblyName in assemblyOrder:
            new.append(assemblies[assemblyName])
        return MolIter(new)

    def __contains__(self, ba):
        """Determine if the specified BindingAssembly is in the OptzymeGroup."""
        # If the input is an integer
        if isinstance(ba, int):
            return -self.length <= ba < self.length
        # Or if it is a string of the BindingAssembly's number
        elif isinstance(ba, str):
            return ba in object.__getattribute__(self, "_assemblyOrder")
        # Or a BindingAssembly
        elif isinstance(ba, DesignGroup):
            for assembly in self:
                if assembly.__eq__(ba):
                    return True
            return False
        # Or a Molecule
        elif isinstance(ba, Molecule):
            for assembly in self:
                for molecule in assembly:
                    if molecule.__eq__(ba):
                        return True
            return False
        # Or a Residue
        elif isinstance(ba, Residue):
            for assembly in self:
                for molecule in assembly:
                    for residue in molecule:
                        if residue.__eq__(ba):
                            return True
            return False
        # Or an Atom
        elif isinstance(ba, Atom):
            for assembly in self:
                for molecule in assembly:
                    for residue in molecule:
                        for atom in residue:
                            if atom.__eq__(ba):
                                return True
            return False
        # Or if it is a list of any of the above
        elif isinstance(ba, list):
            # If it is empty, return false
            if len(ba) == 0:
                return False
            for obj in ba:
                if not self.__contains__(obj):
                    return False
            return True
        # Otherwise, raise an error
        else:
            text = "This is not a valid object to check and see if it is "
            text ++ "contained in an OptzymeGroup:\n" + str(ba)
            raise MoleculeError(text)

    def __getitem__(self, ba):
        """Retrieve the specified BindingAssemblies from the OptzymeGroup."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder() 
        # Get the BindingAssemblies and the order they appear
        assemblies = object.__getattribute__(self, "_assemblies")
        assemblyOrder = object.__getattribute__(self, "_assemblyOrder")
        # If the input is an integer
        if isinstance(ba, int) and -self.length <= ba < self.length:
            return assemblies[assemblyOrder[ba]]
        # Or if the input is a string of the BindingAssembly's name
        elif isinstance(ba, str) and ba in assemblyOrder:
            return assemblies[ba]
        # Or a string requesting a set of Molecules
        elif isinstance(ba, str) and ba.upper() in ['ALL', 'DESIGN', 'TARGET']:
            # Store the molecules in this list
            mols = []
            # Go though the molecules
            for assembly in self: 
                for molecule in assembly:
                    # Gather design molecules if desired
                    if molecule.design:
                        if ba.upper() in ['ALL', 'DESIGN']:
                            mols.append(molecule)
                    # Gather target molecules if desired 
                    else:
                        if ba.upper() in ['ALL', 'TARGET']:
                            mols.append(molecule)
            return mols
        # Or a list requesting a set of Binding Assemblies       
        elif isinstance(ba, list):
            assemblies = []
            for obj in ba:
                assemblies.append(self.__getitem__(obj))
            return assemblies
        # Otherwise, raise an error
        else:
            text = "The following is not an acceptable label to retrieve "
            text += "Binding Assemblies from this OptzymeGroup:\n" 
            text += str(ba)
            raise MoleculeError(text)

    def __setitem__(self, baName, ba):
        """Store the specified BindingAssembly in the OptzymeGroup."""
        # First, make sure ba is a BindingAssembly
        if not isinstance(ba, DesignGroup):
            text = "Only BindingAssembly class objects may be stored "
            text += "in an OptzymeGroup."
            raise MoleculeError(text)
        # Make sure the BindingAssembly's number is being used to store it
        if str(baName) != str(ba.number):  
            text = "A BindingAssembly's name must be used to store it in a "
            text += "OptzymeGroup."
            raise MoleculeError
        # Confirm that the file format, force field, and dimensions attributes
        # match up
        for label in ['fileFormat', 'forceField', 'dimensions']:
            # Retrieve and compare the values of the attributes
            v1 = ba.__getattribute__(label)
            v2 = self.__getattribute__(label)
            if v1 != v2:
                text = "All BindingAssemblies stored in an OptzymeGroup must "
                text += "use the Group's " + label + " attribute."
                raise MoleculeError(text)
        # Make a duplicate of the BindingAssembly and ensure it knows to belong
        # to this opt
        BA = ba.duplicate()
        # If the BindingAssembly's name is not yet in the assemblyOrder list,
        # add it
        if str(BA.number) not in object.__getattribute__(self, \
        "_assemblyOrder"):
            object.__getattribute__(self, \
                                    "_assemblyOrder").append(str(BA.number))
        # Store the BindingAssembly in the OptzymeGroup
        object.__getattribute__(self, "_assembly")[str(BA.number)] = BA
        # Use the reorder function to sort the BindingAssemblies
        self.reorder()
        # Update the length of the OptzymeGroup
        self.length = 0
        
    def __format__(self, how = ''):
        """Create a formatted string of text for the BindingAssemblies within
        the OptzymeGroup."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder()
        # If 'how' hasn't been specified, use the OptzymeGroup's file format
        if how == '':
            how = self.fileFormat
        # Store the formatted text in this string
        text = ''
        # Loop through the Molecules
        for assembly in self:
            for molecule in assembly:
                # Add the formatting for that Molecule
                text += format(molecule, how)
                # If it is PDB related, add END to the end of the string
                if how == "PDB" or (how == "CHARMM" and \
                                    self.fileFormat == "PDB"):
                    text += "END\n"
        return text

    def _summarize(self, exclude = []):
        """Summarize the contents of an OptzymeGroup."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder()
        # Store the summary in this string
        summary = ''
        # Go through each BindingAssembly in the OptzymeGroup
        for i, assembly in enumerate(self):
            # Store all of the Molecules in the assembly
            names = []
            for molecule in assembly:
                names.append(molecule.name)
            # Ensure that at least one of the molecules is eligible for
            # a summary
            include = False
            all = True
            for name in names:
                if name not in exclude:
                    include = True
                else:
                    all = False
            # If there are no Molecules to include in this BindingAssembly, skip
            # it
            if not include:
                summary += "All Molecules in Binding Assembly "
                summary += str(assembly.number) + " of OptzymeGroup " 
                summary += str(self.number) + " are INELIGIBLE.\n\n"
                continue
            # If some of the Molecules are excluded, indicate that
            if not all:
                summary += "Eligible "
            summary += "Molecules"
            # State the OptzymeGroup number and BindingAssembly number
            summary += " in BindingAssembly " + str(i+1)
            summary += " of OptzymeGroup " + str(self.number) + ":"
            # Figure out the Molecules to use in the summary
            molecules = []
            for name in names:
                if name not in exclude:
                    molecules.append(assembly[name])
            # Include information for each Molecule
            for molecule in molecules:
                summary += "\nMolecule " + molecule.name + " contains " + \
                            str(len(molecule)) + " Residues"
            summary += "\n"
        return summary

    def __str__(self):
        """Create a summary of the contents for the OptzymeGroup."""
        return self._summarize()

    def __repr__(self):
        """Call the __str__ function of the OptzymeGroup class."""
        return self._summarize()

    def duplicate(self):
        """Make a new OptzymeGroup that is a match of the current one."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder()
        # Make a new list of BindingAssemblies
        assemblies = []
        for assembly in self:
            assemblies.append(assembly.duplicate())
        # Return the new OptzymeGroup
        new = OptzymeGroup(self.number, assemblies, self.forceField, \
                            self.fileFormat)
        # Set the proper BindingAssembly objectives
        for assembly in self:
            new[str(assembly.number)].objective = \
            self[str(assembly.number)].objective
        return new       

    def __eq__(self, other):
        """Test for equality between two OptzymeGroup class objects."""
        # Check to make sure the assemblyOrder does not need to be resorted
        self.reorder()
        # Test if the other OptzymeGroup has an identical number of
        # BindingAssemblies. If it does not, then they are unequal.
        if not isinstance(other, OptzymeGroup) or other.length != self.length:
            return False
        # Check each BindingAssembly in the OptzymeGroup
        for assembly in self:
            if str(assembly.number) not in other or \
            assembly != other[str(assembly.number)] or \
            assembly.objective != other[str(assembly.number)].objective:
                return False
        return True

    def __ne__(self, other):
        """The inequality check of the OptzymeGroup class."""
        return not self.__eq__(other)

    def _move(self, coordinates, sign):
        """Move the OptzymeGroup."""
        for assembly in self:
            for molecule in assembly:
                molecule._move(coordinates, sign)

    def _rotate(self, rmatrix):
        """Rotate the OptzymeGroup."""
        for assembly in self:
            for molecule in assembly:
                molecule._rotate(rmatrix)

    def output(self, assemblies = "ALL", molecules = "ALL", path = "./", \
               procedure = None, fileFormat = None):
        """Output the specified Molecules in the specified folder."""
        # Create a copy of the assemblies attribute of the OptzymeGroup class
        # object
        Assemblies = object.__getattribute__(self, "_assemblies")
        # If desired, include output molecules from all assemblies
        if assemblies == "ALL":
            for assembly in self:
                assembly.output(molecules, path, procedure, fileFormat)
        # If a list is provided, go through each entry of the list
        elif isinstance(assemblies, list):              
            for assemblyName in assemblies:
                Assemblies[assemblyName].output(molecules, path, procedure, \
                                                fileFormat)
        # If a number or string are provided, use that specific BindingAssembly
        elif isinstance(assemblies, (int, str)):
            Assemblies[assemblies].output(molecules, path, procedure, \
                                          fileFormat)
        # Otherwise, raise an error
        else:   
            text = "The OptzymeGroup class object requires the designated "
            text += "BindingAssemblies to be provided as 'ALL', a list of "
            text += "integers and/or strings, a string, or a number."
            raise MoleculeError(text)

        
