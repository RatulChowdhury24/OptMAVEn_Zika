#!/usr/bin/env python

# The name of the file
__name__ = "IPRO Suite Experiment Module"
# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains the functions needed to create ane run an experiment using
the IPRO Suite of Programs."""

# Import relevant PYTHON modules
import os
import sys
import copy
import math
# Include all contents of the STANDARDS module, which contains all 'default'
# variables, the 'supported' lists, the aminoAcids and backboneAtoms
# dictionaries, the screen_formatting function and the IPRO_Error class
from STANDARDS import *
# Include the I/O modules
import IO_CHECK
import IO_LOADING
import IO_OUTPUT
import IO_VALIDATE
import IO_ASK
# And the sharing and MOLECULES modules
import CHARMM
import SHARING
import ROTAMERS
import MOLECULES
import IPRO_FUNCTIONS
import DOCKING_FUNCTIONS

# Create an error for problems unique to this module
class ExperimentError(IPRO_Error):
    """An error for problems in the EXPERIMENT module."""
    def __init__(self, error = ''):
        """Initialization of the ExperimentError class."""
        IPRO_Error.__init__(self, error)

def translate_attribute_name(name, computer = True):
    """Standardize attribute names"""
    # This function should only be called by other functions that use it
    # correctly. Don't include error checking so that it runs faster.
    # Separate the words in the name by white space
    spacingWords = name.split()
    # Separate words internally by capitalization
    words = []
    for entry in spacingWords:
        # Store the first letter of the word
        word = entry[0]
        # Loop through the characters in entry
        for i in range(1, len(entry)):
            # If this letter is capitalized and the previous letter was not,
            # start a new word
            if entry[i].isupper() and not entry[i-1].isupper():
                words.append(word)
                word = ''
            # Add the character to the word
            word += entry[i]
        # Store the word
        words.append(word)
    # Assemble a new attribute name depending on what method is needed
    attribute = ''
    if computer:
        # have a single string that uses camel case
        for i in range(len(words)):
            if i == 0:
                attribute += words[i].lower()
            else:
                attribute += words[i].capitalize()
    # If the attribute is for the user to read
    else:
        # Loop through the words
        for i in range(len(words)):
            # Put spaces before follow-up words
            if i != 0:
                attribute += " "
            # Make sure some words are capitalized
            if words[i].lower() in ['ipro', 'charmm', 'sa', 'sd', 'lk', 'noe']:
                attribute += words[i].upper()
            else:
                attribute += words[i].capitalize()
    return attribute

def appropriateness(experiment, nonbonded = True):
    """Check for force field and implicit solvation information validity"""
    # Make a list of the Molecules
    molecules = []
    for data in experiment["Molecules"]:
        # This shouldn't be the case, but the function is being written so it
        # can be used as a generic checker, in which case "Molecules" would just
        # be a list of Molecules
        if isinstance(data, MOLECULES.Molecule):
            molecules.append(data)
        # the normal format is a list with the 3rd entry being a Molecule
        else:
            molecules.append(data[2])
    # Add missing Atoms, depending on the Force Field
    if experiment["Force Field"] == "CHARMM":
        try:
            CHARMM.Missing_Atoms(molecules, experiment)
        except CHARMM.CHARMM_Error as error:
            text = "There was an error when CHARMM was tested:\n" + str(error)
            raise ExperimentError(text)
    # If the force field isn't supported
    else:
        text = "The " + str(experiment["Force Field"]) + " is not supported by "
        text += "the appropriateness function"
        raise ExperimentError(text)
    # Try to parameterize the system
    if nonbonded:
        try:
            ROTAMERS.parameterize(molecules, experiment)
        except ROTAMERS.RotamerError as error:
            text = "There was an error when attempting to assign non-bonded "
            text += "energy parameters:\n" + str(error)
            raise ExperimentError(text)
    # Return the list of Molecule structures
    return molecules

def input_validation(experiment):
    """Make sure that the force field and rotamer calculations will work"""
    # Clear the screen
    os.system("clear")
    # Tell the user what is happening
    message = """
The inputs you have provided are now being validated. This should only take a
moment, so please be patient."""
    print screen_formatting(message[1:])
    # Make the Experiment's folder
    os.mkdir(experiment["Folder"])
    # Get the current folder, so we can move back to it when finished
    current = os.getcwd()
    # Move to the Experiment's folder
    os.chdir(experiment["Folder"])
    # copy in the force field and solvation files
    SHARING.copy_standard_files(experiment, current + "/input_files/", False, \
                                True, True)
    # Check the appropriateness of the Molecules, force field, and non-bonded
    # energy parameters
    molecules = appropriateness(experiment, True)
    # If everything worked correctly, the folder can be prepped for the
    # experiment
    # Delete all files, as they're not needed anymore
    names = os.listdir("./")
    for name in names:
        i = os.remove(name)
    # Make an input files folder
    os.mkdir("input_files")
    # Move into that folder
    os.chdir("input_files")
    # Copy in the force field and solvation files
    SHARING.copy_standard_files(experiment, current + "/input_files/", False, \
                                True, True)
    # Move back to the Experiment's folder
    os.chdir(experiment["Folder"])
    # Create a results folder
    os.mkdir("results")
    # Create a structures folder
    os.mkdir("structures")
    # Move into that folder
    os.chdir("structures")
    # Output each Molecule in the structures folder
    for molecule in molecules:
        molecule.output(None, experiment["File Format"], experiment["User"])
    # Move back to the starting folder
    os.chdir(current)

# Create the Experiment class
class Experiment(object):
    """A class for storing the information of an IPRO Suite experiment"""
    def __init__(self, load = True):
        """The initialization of the Experiment class."""
        # Initialize the groups and groupOrder attributes
        object.__setattr__(self, "_groups", {})
        object.__setattr__(self, "_groupOrder", [])
        # Extra section
        print "Are you currently testing?"
        #answer = raw_input()
        answer = False
        if answer in ["yes", "y", "YES", "Y"]:
            self.miracle_load()
        else:    
            # Gather the information about the experiment using the appropriate
            # method
            if load:
                self.load_from_file()
            else:
                self.user_creation()
                
            # Carry out final calculations
            self.finish_creation(load)

    def __setattr__(self, name, value):
        """Control attribute assignment in an Experiment."""
        # If the name is the name of a function, raise an error
        if name in ['__init__', '__setattr__', '__len__', '__iter__', \
        '__contains__', '__getitem__', '__setitem__', '__format__', '__str__', \
        '__repr__', 'output', 'load_from_file', 'user_creation', 'validate', \
        'finish_creation', 'make_DesignGroups']:
            text = name + " is a function of the Experiment class and does not "
            text += "support value assignment."
            raise ExperimentError(text)
        # If the name is for the Design Group information, raise an error
        elif name in ['_groups', '_groupOrder']:
            text = name + " is a private attribute of the Experiment class "
            text += "and does not support value assignment."
            raise ExperimentError(text)
        # If we've gotten past that, generate the internal and external names
        # for the attribute
        internal = translate_attribute_name(name, True)
        external = translate_attribute_name(name, False)
        # Certain attributes should just be accepted and stored
        if external in ['Nonbonded Parameters', "Current", "Energies", \
        "Superimpose", "Target Energies", "Summary", "Last Update", \
        "Docking Groups", "Activity", "Optzyme Information"]:
            object.__setattr__(self, internal, value)
        # The length is the number of Design Groups in the Experiment
        elif external == "Length":
            object.__setattr__(self, internal, len(self._groupOrder))
        # If the attribute is a piece of basic information every IPRO Suite
        # Experiment needs to have
        elif external in ['User', 'Type', 'Name', 'Force Field', 'File Format',\
                          'Folder']:
            IO_CHECK.basic_info(external, value)
            object.__setattr__(self, internal, value)
        # if the attribute is related to the use of the CHARMM force field
        elif external in ['CHARMM Topology Files', 'CHARMM Parameter Files', \
                          'CHARMM Energy Terms', 'CHARMM Iterations', 'Lj Phi']:
            IO_CHECK.CHARMM_info(external, value)
            object.__setattr__(self, internal, value)
        # If the attribute is related to how to use implicit solvation
        elif external in ['Use Solvation', 'Relaxation Solvation', \
        'Perturbation Solvation', 'Energy Solvation', 'Solvation Type', \
        'LK Solvation Files']:
            IO_CHECK.solvation_info(external, value)
            object.__setattr__(self, internal, value)
        # If the attribute relates to the use of Rotamers
        elif external in ['Rotamer Library', 'Max Rotamer Number', \
        'Rotamer Window', 'Packing Method', 'Packing Selection', \
        'Packing Cutoff']:
            IO_CHECK.rotamer_info(external, value)
            object.__setattr__(self, internal, value)
        # If the attribute is for structure restraints
        elif external == "Restraints":
            IO_CHECK.Restraints(value, self)
            object.__setattr__(self, internal, value)
        # If the attribute relates to how to run docking
        elif external in ['Docking Frequency', 'Docking Iterations', \
        'Docking SD Movements', 'Docking SD Rotations', 'Docking Start Temp', \
        'Docking End Temp']:
            IO_CHECK.docking_info(external, value, self)
            object.__setattr__(self, internal, value)
        # If the attribute is involved with how to do structure refinements
        elif external in ['Do Refinement', 'Refinement Iterations', \
                          'Ensemble Number']:
            IO_CHECK.refinement_info(external, value)
            object.__setattr__(self, internal, value)
        # If this is the list of Molecules used in the Experiment
        elif external == "Molecules":
            IO_CHECK.Molecules(value)
            object.__setattr__(self, internal, value)
        # If this is the dictionary of Design Positions
        elif external == "Design Positions":
            IO_CHECK.DesignPositions(value, self)
            object.__setattr__(self, internal, value)
        # If this is the dictionary of Epitope Positions
        elif external == "Epitope Positions":
            IO_CHECK.EpitopePositions(value, self)
            object.__setattr__(self, internal, value)
        # If this is the list of Design Groups
        elif external == "Design Groups":
            IO_CHECK.DesignGroups(value, self)
            object.__setattr__(self, internal, value)
        # If this is the list of Optzyme Groups
        elif external in ["Optzyme Groups"]:
            IO_CHECK.OptzymeGroups(value, self)
            object.__setattr__(self, internal, value)
        # If this is the dictionary of permitted kinds of amino acids for design
        # positions
        elif external == "Permitted Kinds":
            IO_CHECK.PermittedKinds(value, self)
            object.__setattr__(self, internal, value)
        # If the information is related to how to run IPRO
        elif external in ['IPRO Iterations', 'IPRO Annealing Temperature', \
                          'Annealing Sharing', 'Energy Calculation']:
            IO_CHECK.IPRO_info(external, value)
            object.__setattr__(self, internal, value)
        # If the information is about the mutants in a Mutator experiment
        elif external == "Mutants":
            IO_CHECK.Mutants(value, self)
            object.__setattr__(self, internal, value)
        # If the information is about pairs of Molecules that are dimers
        elif external == "Dimers":
            IO_CHECK.Dimers(value, self)
            object.__setattr__(self, internal, value)
        # If the attribute is anything else, raise an error
        else:
            text = str(name) + " is not an attribute that may be stored in an "
            text += "IPRO Suite Experiment class object."
            raise ExperimentError(text)

    def __len__(self):
        """The number of Design Groups in the Experiment."""
        # The length attribute is controlled by __setattr__
        self.length = 0
        return self.length

    def __iter__(self):
        """The iterator of the Experiment class."""
        # Make sure the Experiment's length is correct
        self.length = 0
        # Create a new list of the Design Groups
        new = []
        for i in self._groupOrder:
            new.append(self._groups[i])
        return MOLECULES.MolIter(new)

    def __contains__(self, name):
        """Determine if the specified attribute is in the Experiment."""
        # If the name is a string
        if isinstance(name, str):
            # Get the computer version of the name
            attribute = translate_attribute_name(name, True)
            # Check the Experiment's dictionary
            return attribute in self.__dict__
        # If the name is an integer
        elif isinstance(name, int):
            try:
                return name in self._groupOrder
            except AttributeError:
                return False
        # If the name is anything else
        else:
            return False

    def __getitem__(self, name):
        """Retrieve the specified Experiment attribute"""
        # If the name is a string
        if isinstance(name, str):
            # Convert the name to the format the computer expects
            attribute = translate_attribute_name(name, True)
            # Try to get it
            try:
                return object.__getattribute__(self, attribute)
            # If it doesn't exist
            except AttributeError:
                text = name + " is not an attribute of the Experiment."
                raise ExperimentError(text)
        # If it is an integer, it is an attempt to access a Design Group
        elif isinstance(name, int):
            if name in range(len(self) + 1):
                return self._groups[name]
            else:
                text = str(name) + " is not a Design Group in this Experiment."
                raise ExperimentError(text)
        # Otherwise, raise an error
        else:
            text = str(name) + " is not an attribute of the Experiment."
            raise ExperimentError(text)

    def __setitem__(self, name, value):
        """Store an attribute in an experiment."""
        # If it is a string
        if isinstance(name, str):
            # attempt to store it using __setattr__
            self.__setattr__(name, value)
        else:
            text = str(name) + " is not an attribute that may be stored in an "
            text += "IPRO Suite Experiment class object."
            raise ExperimentError(text)
    
    def validate(self):
        """Make sure an Experiment has all information that it needs."""
        # Start by validating the basic information
        IO_VALIDATE.basic_info(self)
        # If all of that information was found, other information can be
        # checked successfully, too. Store other errors here
        errors = ''
        # Force field
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            if self["Force Field"] == "CHARMM":
                errors += IO_VALIDATE.CHARMM_info(self)
            else:
                errors += "\nThe validate function does not support the "
                errors += str(self["Force Field"]) + " force field."
        # How to run Docking
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.docking_info(self)
        # How to run IPRO
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.IPRO_info(self)
        # How to do structure refinements
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.refinement_info(self)
        # How to use Rotamers
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.rotamer_info(self)
        # how to use implicit solvation
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.solvation_info(self)
        # What the Molecules are
        if self["Type"] in ['IPRO', 'Mutator', 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.Molecules(self)
        # The remaining checks are all dependent on the Molecule information
        # being correct
        if errors != '':
            raise ExperimentError(errors[1:])
        # If there are Dimers
        if self["Type"] in ["IPRO", "Mutator", 'OptMAVEn', 'OptZyme']:
            errors += IO_VALIDATE.Dimers(self)
        # Design Groups
        if self["Type"] in ["IPRO", "Mutator", 'OptMAVEn']:
            errors += IO_VALIDATE.DesignGroups(self)
        # Optzyme Groups
        if self["Type"] in ["OptZyme"]:
            errors += IO_VALIDATE.OptzymeGroups(self)
        # Design Positions
        if self["Type"] in ["IPRO", "OptZyme"]:
            errors += IO_VALIDATE.DesignPositions(self)
        # The permitted kinds of amino acids is dependent on the design
        # positions
        if errors != '':
            raise ExperimentError(errors[1:])
        if self["Type"] in ["IPRO", "OptZyme"]:
            errors += IO_VALIDATE.PermittedKinds(self)
        # The mutants
        if self["Type"] in ["Mutator"]:
            errors += IO_VALIDATE.Mutants(self)
        # Any structure restraints
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            errors += IO_VALIDATE.Restraints(self)
        if errors != '':
            raise ExperimentError(errors[1:])

    def miracle_load(self):
        """Load the contents of an Experiment from a file"""
        # Start by loading the data
        data = IO_LOADING.load_Experiment_Details()
        # Now go through the different types of information
        IO_LOADING.basic_info(self, data)
        # Store any remaining errors here
        errors = ''
        # Get the force field information (those functions check for the
        # experiment's force field at their beginning, so that doesn't have to
        # happen here)
        if self["Type"] in ['IPRO', 'Mutator', "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.CHARMM_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.docking_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.IPRO_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.refinement_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.rotamer_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.solvation_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.Molecules(self, data)
        # The following are dependent on the Molecules being loaded
        # successfully, so if there are errors be done
        if errors != '':
            raise ExperimentError(errors[1:])
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            errors += IO_LOADING.Dimers(self, data)
        if self["Type"] in ["OptZyme"]:
            IO_ASK.OptzymeGroups(self)
            self.make_OptzymeGroups()
        if self["Type"] in ["OptZyme"]:
            text = IO_OUTPUT.OptzymeGroups(self)
        file = open("output.txt", "a")
        file.write(text)
        file.close()
        sys.exit(0)

    def load_from_file(self):
        """Load the contents of an Experiment from a file"""
        # Start by loading the data
        data = IO_LOADING.load_Experiment_Details()
        # Now go through the different types of information
        IO_LOADING.basic_info(self, data)
        # Store any remaining errors here
        errors = ''
        # Get the force field information (those functions check for the
        # experiment's force field at their beginning, so that doesn't have to
        # happen here)
        if self["Type"] in ['IPRO', 'Mutator', "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.CHARMM_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.docking_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.IPRO_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.refinement_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.rotamer_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.solvation_info(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            errors += IO_LOADING.Molecules(self, data)
        # The following are dependent on the Molecules being loaded
        # successfully, so if there are errors be done
        if errors != '':
            raise ExperimentError(errors[1:])
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            errors += IO_LOADING.Dimers(self, data)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn"]:
            errors += IO_LOADING.DesignGroups(self, data)
            if errors == '':
                self.make_DesignGroups()
        if self["Type"] in ["OptZyme"]:
            errors += IO_LOADING.OptzymeGroups(self, data)
            if errors == '':
                self.make_OptzymeGroups()
        if self["Type"] in ["IPRO", "OptZyme"]:
            errors += IO_LOADING.DesignPositions(self, data)
        #if self["Type"] in ["OptMAVEn"]:
        #    errors += IO_LOADING.EpitopePositions(self, data)
        # The following are dependent on having accurate residue info
        if errors != '':
            raise ExperimentError(errors[1:])
        if self["Type"] in ["IPRO", "OptZyme"]:
            errors += IO_LOADING.PermittedKinds(self, data)
        if self["Type"] in ["Mutator"]:
            errors += IO_LOADING.Mutants(self, data)
        if errors != '':
            raise ExperimentError(errors[1:])
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            errors += IO_LOADING.Restraints(self, data)
        if errors != '':
            raise ExperimentError(errors[1:])
        # Validate that the experiment has all of the information that it should
        self.validate()

    def output(self):
        """Write the Experiment's information to an Experiment Details file"""
        # There's no need to validate the information - that happens in the
        # OUTPUT module
        text = IO_OUTPUT.basic_info(self)
        # Include the other information in a sensible order
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]: 
            text += IO_OUTPUT.IPRO_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.CHARMM_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.solvation_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.rotamer_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.refinement_info(self)
        if self["Type"] in ['IPRO', 'Mutator', "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.docking_info(self)
        # Now start including information about Molecules and how they're used
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            text += IO_OUTPUT.Molecules(self)
            # Molecules has the dimer info output at the same time
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn"]:
            text += IO_OUTPUT.DesignGroups(self)
        if self["Type"] in ["OptZyme"]:
            text += IO_OUTPUT.OptzymeGroups(self)
        if self["Type"] in ["IPRO", "OptZyme"]:
            text += IO_OUTPUT.DesignPositions(self)
        if self["Type"] in ["IPRO", "OptZyme"]:
            text += IO_OUTPUT.PermittedKinds(self)
        if self["Type"] in ["Mutator"]:
            text += IO_OUTPUT.Mutants(self)
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            text += IO_OUTPUT.Restraints(self)
        #if self["Type"] in ["OptMAVEn"]:
        #    text += IO_OUTPUT.EpitopePositions(self)
        # Write this to an appropriate file
        try:
            f = open(self["Folder"] + "Experiment_Details.txt", "w")
        except IOError:
            f = open("Experiment_Details.txt", "w")
        f.write(text)
        f.close()

    def user_creation(self):
        """Ask the user for the information about how to run an experiment."""
        # Start by asking for the basic information
        IO_ASK.basic_info(self)
        # Start getting other information based on the type of the Experiment
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IO_ASK.Standard_Molecules(self)
        #if self["Type"] in ["OptMAVEn"]:
        #    IO_ASK.antigen_molecules(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            if self["Force Field"] == "CHARMM":
                IO_ASK.CHARMM_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IO_ASK.solvation_info(self)
        # Make sure those inputs will all work together
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            input_validation(self)
        # Now that we know they do, start getting more information about how to
        # run the Experiment
        if self["Type"] in ["IPRO", "OptZyme"]:
            IO_ASK.DesignPositions(self)
            IO_ASK.PermittedKinds(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn"]:
            IO_ASK.DesignGroups(self)
            self.make_DesignGroups()
        if self["Type"] in ["OptZyme"]:
            IO_ASK.OptzymeGroups(self)
            self.make_OptzymeGroups()
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IO_ASK.IPRO_info(self)
        if self["Type"] in ["Mutator"]:
            IO_ASK.mutator_info(self)
            IO_ASK.Mutants(self)
        if self["Type"] in ["IPRO", "OptMAVEn", "OptZyme"]:
            IO_ASK.refinement_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IO_ASK.rotamer_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IO_ASK.docking_info(self)
        if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
            IO_ASK.Restraints(self)
        #if self["Type"] in ["OptMAVEn"]:
            #IO_ASK.OptMaven_info(self)
        #    IO_ASK.Epitope(self)

    def make_DesignGroups(self, load = True):
        """Make the Experiment's Design Groups"""
        # Tell the user what is happening
        if not load:
            os.system("clear")
            text = "The actual MOLECULES class objects that make up the Design "
            text += "Groups are now being generated."
            print screen_formatting(text)
        # Make an initial Design Group that contains every Molecule
        molecules = []
        for data in self["Molecules"]:
            molecules.append(data[2])
        group = MOLECULES.DesignGroup(0, molecules, self["Force Field"], \
                                      self["File Format"])
        # Store the Design Group in the Experiment
        object.__setattr__(self, "_groups", {})
        object.__setattr__(self, "_groupOrder", [])
        self._groups[0] = group
        # Go through the Design Groups
        for i, group in enumerate(self["Design Groups"]):
            new = []
            for molecule in molecules:
                if molecule.design or molecule.name in group[1:]:
                    new.append(molecule)
            # Make a new Design Group out of this
            design = MOLECULES.DesignGroup(i+1, new, self["Force Field"], \
                                           self["File Format"])
            design.objective = group[0]
            # Store the Design Group
            self._groupOrder.append(i+1)
            self._groups[i+1] = design

    def make_OptzymeGroups(self, load = True):
        """Make the Experiment's Optzyme Groups"""
        # Tell the user what is happening
        if not load:
            os.system("clear")
            text = "The actual MOLECULES class objects that make up the "
            text += "Optzyme Groups are now being generated."
            print screen_formatting(text)
        # Construct an initial Binding Assembly that contains every Molecule
        molecules = []
        for data in self["Molecules"]:
            molecules.append(data[2])
        assembly = MOLECULES.DesignGroup(0, molecules, self["Force Field"], \
                                         self["File Format"])
        # Store the Binding Assembly in the Experiment
        object.__setattr__(self, "_groups", {})
        object.__setattr__(self, "_groupOrder", [])
        self._groups[0] = assembly
        # Store the unique Binding Assemblies listed within the Optzyme Groups
        unique = []
        # Go through the Optzyme Groups
        for group in self["Optzyme Groups"]:
            for assembly in group[3]:
                # If this is a new Binding Assembly, construct it
                if assembly not in unique:
                    unique.append(assembly)
                    # If this is a new Binding Assembly, store the Molecules
                    new = []
                    for molecule in molecules:
                        if molecule.design or molecule.name in assembly[1:]:
                            new.append(molecule)
                    # Make the new Binding Assembly
                    bindingAssembly = MOLECULES.DesignGroup(len(unique), new, \
                                       self["Force Field"], self["File Format"])
                    bindingAssembly.objective = assembly[0]
                    # Store the Binding Assembly
                    self._groupOrder.append(len(unique))
                    self._groups[len(unique)] = bindingAssembly
        # Store the Optzyme Groups in a list
        groups = []
        # Go through the Optzyme Groups
        for i, group in enumerate(self["Optzyme Groups"]):
            # Find the corresponding Binding Assemblies
            assemblies = []
            # Go through each Assembly within the Optzyme Group
            for ba in group[3]:
                # State whether the matching Assembly has been found
                found = False
                # Go through the stored Binding Assemblies 
                for assembly in self:
                    # Get the list of Target Molecules in the assembly
                    target = []
                    for mol in assembly:
                        if not mol.design:
                            target.append(mol.name)
                    # See if the list matches the one provided from the Optzyme
                    # Group
                    if target.sort() == ba[1:].sort():
                        found = True
                        break
                # Raise an error if the Binding Assembly was not found        
                if not found:           
                    text = "A matching Binding Assembly for the inputs "
                    text += "provided was not found"
                    raise ExperimentError(text)
                # Otherwise, store the information in the 'assemblies' list
                else:
                    assemblies.append(assembly)
            # Generate the Optzyme Group
            optzymeGroup = MOLECULES.OptzymeGroup(i+1, assemblies, \
                            self["Force Field"], self["File Format"])
            # Store the attributes of the Optzyme Group
            optzymeGroup.property = group[0]
            optzymeGroup.objective = group[1]
            if group[0] != "kcat":
                optzymeGroup.weight = group[2]
            # Append this Optzyme Group to the list of Groups
            groups.append(optzymeGroup)
        # Store the list of Optzyme Groups
        self["Optzyme Information"] = groups

    def finish_creation(self, load = True):
        """Finish creating an Experiment"""
        # If this is not being loaded, tell the user what is happening
        if not load:
            os.system("clear")
            text = "Some final calculations are being carried out to finish "
            text += "preparing the Experiment for use. Please wait a moment."
            print screen_formatting(text)
        # Initialize all of the attributes that got ignored in the __setattr__
        # function
        # Start with the ones that are simple values
        self["Summary"] = ''
        self["Last Update"] = -1
        self["Activity"] = "Standard"
        # Store the Design Position and permitted kinds
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            for group in self:
                for molecule in group:
                    if self["Type"] in ["IPRO", "Mutator", "OptZyme"]:
                        molecule.identify_design_positions(self["Design Positions"])
                        # if there are permitted kinds of amino acids
                        if "Permitted Kinds" in self:
                            molecule.identify_permittedKinds(self.permittedKinds)
                        else:
                            molecule.identify_permittedKinds()
                    if self["Type"] in ["OptMAVEn"] and molecule.design:
                    # Set the designs attribute of all the residues in the molecule are True and permittedKinds attribute to 20 AAs
                        spots = {}
                        names = []
                        #permitted = {}
                        for r in molecule:
                            names.append(r.name)
                            #permitted[molecule.name][r.name] = aminoAcids["PDB"]
                        spots[molecule.name] = names
                        molecule.identify_design_positions(spots)
                        molecule.identify_permittedKinds()
        # Store Atoms that may NEVER move in the Molecules
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            if "Restraints" in self and "Fixed Atoms" in self["Restraints"]:
                data = self["Restraints"]["Fixed Atoms"]
                # Go through the design group specifications
                for gn in data:
                    groups = []
                    for group in self:
                        if gn in ['all', group.number]:
                            groups.append(group)
                    for group in groups:
                        for molecule in group:
                            molecule.identify_fixedAtoms(data[gn])
        # Initialize the Docking Groups
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            DOCKING_FUNCTIONS.calculate_docking_groups(self)
        # Create the Current dictionary and its structures
        self["Current"] = {}
        for group in self:
            self["Current"][group.number] = {}
            for molecule in group:
                new = molecule.duplicate()
                self["Current"][group.number][new.name] = new
        # Initialize the Optzyme Groups
        if self["Type"] in ["OptZyme"]:
            IPRO_FUNCTIONS.initialize_optzyme_groups(self)
        # Create the Energies dictionary
        self["Energies"] = {}
        for group in self:
            self["Energies"][group.number] = {}
        # If Molecules should be superimposed during relaxations
        if self["Type"] in ["IPRO", "Mutator", "OptMAVEn", "OptZyme"]:
            IPRO_FUNCTIONS.do_superimpose(self)
