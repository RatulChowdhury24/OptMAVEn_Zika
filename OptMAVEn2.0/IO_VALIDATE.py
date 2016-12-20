#!/usr/bin/env python

# The name of the file
__name__ = "IPRO Suite Input / Output Attribute Validation"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for validating IPRO Suite Experiment attributes"""

# Include the standard python modules
import os
import sys
# Include the CHECK and FUNCTIONS I/O module
import IO_CHECK as CHECK
import IO_FUNCTIONS as FUNCTIONS

# Have functions for each of the various types of information
def basic_info(experiment):
    """Check that the Experiment has all basic information."""
    # Store any generated errors here
    errors = ''
    # Loop through the attributes
    for name in ['User', 'Type', 'Name', 'File Format', 'Force Field','Folder']:
        # If the attribute is missing
        if name not in experiment:
            text += "\nThe " + name + " attribute is missing."
        else:
            # If it's not missing, check its value
            try:
                CHECK.basic_info(name, experiment[name])
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # Unlike in other validation functions, if there is a problem with basic
    # information the experiment should immediately die instead of searching for
    # other errors
    if errors != '':
        raise FUNCTIONS.IPRO_IOError(errors[1:])

def CHARMM_info(experiment):
    """Make sure the Experiment has all needed CHARMM information."""
    # Store any errors here
    errors = ''
    # Only do the search if the experiment is using the CHARMM force field
    if experiment["Force Field"] == "CHARMM":
        # Search through the attributes
        for name in ['CHARMM Topology Files', 'CHARMM Parameter Files', \
                     'CHARMM Energy Terms', 'CHARMM Iterations', 'Lj Phi']:
            # If the attribute is missing
            if name not in experiment:
                errors += "\nThe " + name + " attribute is missing"
            # Otherwise, check the values
            else:
                try:
                    CHECK.CHARMM_info(name, experiment[name])
                except FUNCTIONS.IPRO_IOError as error:
                    errors += str(error)
    # Return the errors to where the validation function was called from
    return errors

def docking_info(experiment):
    """Make sure all Docking data for an Experiment is valid."""
    # Store the errors here
    errors = ''
    # Go through the attributes
    for attribute in ['Docking Frequency', 'Docking Iterations', \
    'Docking SD Movements', 'Docking SD Rotations', 'Docking Start Temp', \
    'Docking End Temp']:
        # If the attribute is missing
        if attribute not in experiment:
            errors += "\nThe " + attribute + " attribute is missing"
        # Otherwise, check its value
        else:
            try:
                CHECK.docking_info(attribute, experiment[attribute], experiment)
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # Return the errors to the user
    return errors

def IPRO_info(experiment):
    """Check the values of attributes for how to run IPRO."""
    # Store errors here
    errors = ''
    # Go through the IPRO attributes
    for name in ['IPRO Iterations', 'IPRO Annealing Temperature', \
                 'Annealing Sharing', 'Energy Calculation']:
        # If the attribute is missing
        if name not in experiment:
            errors += "\nThe " + name + " attribute is missing"
        # Otherwise, try to check it
        else:
            try:
                CHECK.IPRO_info(name, experiment[name])
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    return errors

def refinement_info(experiment):
    """Make sure the Experiment has valid refinement information."""
    # Store errors here
    errors = ''
    # Go through the attributes
    for attribute in ['Do Refinement', 'Refinement Iterations', \
                      'Ensemble Number']:
        # If the attribute is missing
        if attribute not in experiment:
            errors += "\nThe " + attribute + " attribute is missing."
        # Otherwise, check the value
        else:
            try:
                CHECK.refinement_info(attribute, experiment[attribute])
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    # Return the error information to the user
    return errors

def rotamer_info(experiment):
    """Confirm that all needed Rotamer information is available and correct"""
    # Store errors here
    errors = ''
    # Go through the attributes
    for attribute in ['Rotamer Library', 'Rotamer Window','Max Rotamer Number',\
                      'Packing Cutoff', 'Packing Method', 'Packing Selection']:
        # If the attribute is missing
        if attribute not in experiment:
            errors += "\nThe " + attribute + " attribute is missing."
        # Otherwise, check the value
        else:
            try:
                CHECK.rotamer_info(attribute, experiment[attribute])
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
    return errors

def solvation_info(experiment):
    """Check the values of all needed implicit solvation attributes"""
    # Store the errors here
    errors = ''
    # And whether or not a problem has been encountered here
    problem = False
    # If there is no information about how to use solvation
    if "Use Solvation" not in experiment:
        errors += "\nThe Use Solvation attribute is missing"
        problem = True
    # Otherwise check it, then move on
    else:
        try:
            CHECK.solvation_info("Use Solvation", experiment["Use Solvation"])
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
            problem = True
    # If there's no problem and solvation is being used
    if not problem and experiment["Use Solvation"]:
        # Check for the type of solvation
        if "Solvation Type" not in experiment:
            errors += "\nThe Solvation Type attribute is missing"
            problem = True
        else:
            try:
                CHECK.solvation_info("Solvation Type", \
                                     experiment["Solvation Type"])
            except FUNCTIONS.IPRO_IOError as error:
                errors += str(error)
                problem = True
        # if there's no problem and LK implicit solvation is being used
        if not problem and experiment["Solvation Type"] == "Lazaridis-Karplus":
            # Check the LK solvation files
            name = "LK Solvation Files"
            if name not in experiment:
                errors += "\nThe " + name + " attribute is missing"
            else:
                try:
                    CHECK.solvation_info(name, experiment[name])
                except FUNCTIONS.IPRO_IOError as error:
                    errors += str(error)
    # Return the error information
    return errors

def Molecules(experiment):
    """Check whether or not the Experiment has valid Molecule information"""
    # If the Molecules information is missing
    if "Molecules" not in experiment:
        return "\nThe Molecules information is missing"
    # Otherwise, just check them
    errors = ''
    try:
        CHECK.Molecules(experiment["Molecules"])
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors

def Dimers(experiment):
    """If there is Dimer information, check it"""
    errors = ''
    # There is no guarantee that there will be dimer information
    if "Dimers" in experiment:
        try:
            CHECK.Dimers(experiment["Dimers"], experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    return errors

def DesignPositions(experiment):
    """Confirm that the Design Position information is acceptable."""
    if "Design Positions" not in experiment:
        return "\nThe Design Positions information is missing"
    errors = ''
    try:
        CHECK.DesignPositions(experiment["Design Positions"], experiment)
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors

def EpitopePositions(experiment):
    """Confirm that the Epitope Position information is acceptable."""
    if "Epitope Positions" not in experiment:
        return "\nThe Epitope Positions information is missing"
    errors = ''
    try:
        CHECK.EpitopePositions(experiment["Epitope Positions"], experiment)
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors

def DesignGroups(experiment):
    """Make sure the declaration of the Design Groups is acceptable"""
    if "Design Groups" not in experiment:
        return "\nThe Design Groups information is missing"
    errors = ''
    try:
        CHECK.DesignGroups(experiment["Design Groups"], experiment)
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors

def OptzymeGroups(experiment):
    """Make sure the information about the Optzyme Groups is acceptable"""
    if "Optzyme Groups" not in experiment:
        return "\nThe Optzyme Groups information is missing"
    errors = ''
    try:
        CHECK.OptzymeGroups(experiment["Optzyme Groups"], experiment)
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors      

def PermittedKinds(experiment):
    """If there is information about how Residues may mutate, check it"""
    errors = ''
    if "Permitted Kinds" in experiment:
        try:
            CHECK.PermittedKinds(experiment["Permitted Kinds"], experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    return errors

def Mutants(experiment):
    """Confirm the usefulness of the Mutation information"""
    if "Mutants" not in experiment:
        return "\nThe Mutants information is missing"
    errors = ''
    try:
        CHECK.Mutants(experiment["Mutants"], experiment)
    except FUNCTIONS.IPRO_IOError as error:
        errors += str(error)
    return errors

def Restraints(experiment):
    """If an Experiment has restraints, check each type."""
    # Store errors here
    errors = ''
    # Only do this if there ARE restraints
    if "Restraints" in experiment:
        try:
            CHECK.Restraints(experiment["Restraints"], experiment)
        except FUNCTIONS.IPRO_IOError as error:
            errors += str(error)
    # Return the errors
    return errors
