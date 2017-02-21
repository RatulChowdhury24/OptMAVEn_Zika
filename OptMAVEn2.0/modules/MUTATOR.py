#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Mutation Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions specific to running the Mutator program - it
creates mutants."""

# Import normal PYTHON modules
import os
import sys
# Include IPRO Suite modules
from STANDARDS import *
import IPRO_FUNCTIONS
import ROTAMERS
import SHARING

# An error class for problems specific to this module
class MutatorError(IPRO_Error):
    """An error for problems in the MUTATOR module"""
    def __init__(self, error = ''):
        """The initialization of the MutatorError class."""
        IPRO_Error.__init__(self)

def permission_setter(experiment, gn, mn):
    """Set the permissions for a Design Group before an initial mutation"""
    # First, go through every Residue in the Design Group and freeze it in place
    for molecule in experiment[gn]:
        for residue in molecule:
            residue.freedom = "FIXED"
            residue.permission = "FIXED"
    # Go through the mutations in this mutant
    for mutation in experiment["Mutants"][mn - 1]:
        residue = experiment[gn][mutation[0]][mutation[1]]
        residue.freedom = "RESTRAINED"
        residue.permission = "ROTAMER"
        residue.permittedKinds = [mutation[2]]
    # Repack things around that based on distance
    text = experiment["Packing Method"]
    experiment["Packing Method"] = "Distance"
    IPRO_FUNCTIONS.distance_setting(experiment, experiment[gn])
    experiment["Packing Method"] = text
    # Make sure that dimer information matches
    IPRO_FUNCTIONS.Dimers(experiment, experiment[gn], {})

def mutate_DesignGroups(experiment, mn):
    """Make the initial mutants of the Design Groups"""
    # Move into the folder to do the initial calculations in 
    folder = "initial_mutant" + str(mn)
    os.chdir(folder)
    # Loop through the Design Groups
    for group in experiment:
        # The calculations will be stored in this folder
        folder = "Group" + str(group.number)
        # Try to claim the calculations
        do = SHARING.claim_calculations(folder)
        # If this processor is doing those calculations
        if do:
            # Time stamp when this started
            experiment["Summary"] = "Started" + SHARING.time_stamp()
            # Move into the folder
            os.chdir(folder)
            # Copy in the C++ and force field files
            SHARING.copy_standard_files(experiment)
            # Use the Current structures
            for molecule in experiment[group.number]:
                text =format(experiment["Current"][group.number][molecule.name])
                molecule.load(text)
            # Set the permissions for the Molecules
            permission_setter(experiment, group.number, mn)
            # Mutate the Residues
            refinement = IPRO_FUNCTIONS.Optimal_Rotamers(experiment, \
                                                         group.number)
            refinement = IPRO_FUNCTIONS.Relaxation(experiment, group.number, \
                                                   True)
            energies, refinement = IPRO_FUNCTIONS.Calculate_Energy(experiment, \
                                                                   group.number) 
            # Start sharing
            SHARING.Start(experiment)
            # Write a brief summary file
            name = SHARING.summary_name(SHARING.get_current())
            f = open(name, "w")
            f.write(experiment["Summary"])
            f.close()
            # Move up a folder
            os.chdir("../")
            # Store the structures in the Current dictionary
            IPRO_FUNCTIONS.store_structures(experiment, group.number)
            IPRO_FUNCTIONS.store_energies(experiment, energies, group.number)
            # Write the structures to the Current folder
            SHARING.output_Current(experiment, "./Current/", group.number)
            SHARING.output_Energies(experiment, "./Current/", group.number)
            # End sharing
            SHARING.End(experiment)

def check_finish(experiment):
    """Check to see if the initial structures for a mutant have been made yet"""
    # Start sharing
    SHARING.Start(experiment)
    # Keep track of whether everything is finished
    finished = True
    # Check every Design Group for a summary file
    for group in experiment:
        label = "Group" + str(group.number)
        if label + "_Summary.txt" not in os.listdir(label):
            finished = False
            break
    return finished

def Finish(experiment, mn):
    """Finish the initial creation of a mutant"""
    # Sharing was started in check finish
    # Load the structures and energies
    SHARING.update_Current(experiment, "./Current/")
    SHARING.update_Energies(experiment, "./Current/")
    # Create a brief summary of the information
    experiment["Summary"] = '\nMutant ' + str(mn) + " Creation\n"
    f = open("Group1/Group1_Summary.txt", "r")
    for line in f:
        if line.startswith("Started"):
            experiment["Summary"] += line
            break
    f.close()
    # Include information about the energies
    for group in experiment:
        experiment["Summary"] += \
            SHARING.format_energies(experiment["Energies"][group.number], \
                                    group.number, False)
    # Put this in the overall summary
    os.chdir(experiment["Folder"])
    name = SHARING.summary_name(SHARING.get_current())
    f = open(name, "a")
    f.write(experiment["Summary"])
    f.close()
    # Make a folder to do the structure refinements in
    folder = "mutant" + str(mn)
    os.mkdir(folder)
    folder += "/Current/"
    os.mkdir(folder)
    SHARING.output_Current(experiment, folder)
    SHARING.output_Energies(experiment, folder)
    f = open(folder + "iteration.txt", "w")
    f.write(str(mn))
    f.close()
    # Delete the current folder
    name = "initial_mutant" + str(mn)
    os.system("rm -rf " + name)
    # End sharing
    SHARING.End(experiment)

def initial_check(experiment, mn):
    """Do an initial check about doing calculations for a mutant"""
    # Start sharing
    SHARING.Start(experiment)
    # This is the name of the folder we're looking for
    folder = "mutant" + str(mn)
    # If the results are already completed
    if folder in os.listdir(experiment["Folder"] + "results/"):
        SHARING.End(experiment)
        return False
    # Or if the refinement is ongoing
    elif folder in os.listdir(experiment["Folder"]):
        SHARING.End(experiment)
        return False
    # If the initial folder for that mutant does not yet exist, make it
    folder = "initial_" + folder
    if folder not in os.listdir(experiment["Folder"]):
        os.mkdir(folder)
        os.mkdir(folder + "/Current")
    # Load the wildtype structures
    SHARING.update_Current(experiment, experiment["Folder"] + \
                           "results/wildtype/")
    SHARING.update_Energies(experiment, experiment["Folder"] + \
                            "results/wildtype/")
    # End sharing
    SHARING.End(experiment)
    return True

def DO(experiment, mn):
    """Make a particular Mutant"""
    # Determine if the mutant should be made
    do = initial_check(experiment, mn)
    # If it should be, do so
    if do:
        # Mutate the Design Groups
        mutate_DesignGroups(experiment, mn)
        # Check to see if everything is finished
        finished = check_finish(experiment)
        # If everything is done, start the refinement calculations for the
        # mutant
        if finished:
            Finish(experiment, mn)
        else:
            # Otherwise, go back to the Experiment's folder, end sharing, and
            # wait for the initial folder to be deleted
            os.chdir(experiment["Folder"])
            SHARING.End(experiment)
            SHARING.Wait("initial_mutant" + str(mn))
