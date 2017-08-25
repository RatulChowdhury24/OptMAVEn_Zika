#!/usr/bin/env python

# The name
__name__ = "IPRO Suite Refinement Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions to carry out the extensive, ensemble-based
structure refinements that are part of the IPRO Suite of Programs."""

# Include standard python modules
import os
import sys
import math
# And IPRO Suite information
from STANDARDS import *
import SHARING
import IPRO_FUNCTIONS
import MOLECULES

class RefinementError(IPRO_Error):
    """An error for problems in the Refinement Module"""
    def __init__(self, error = ''):
        """The initialization of the RefinementError class"""
        IPRO_Error.__init__(self, error)

def ensemble_refinement(experiment, gn, en):
    """Do a refinement of the specified Design Group"""
    # The refinement folder is here
    folder = "Ensemble" + str(en)
    # Move into the folder
    os.chdir(folder)
    # Do an appropriate number of iterations
    iteration = 0
    experiment["Last Update"] = -1
    while iteration <= experiment["Refinement Iterations"]:
        print "Doing refinement"
        iteration = IPRO_FUNCTIONS.IPRO_ITERATION(experiment, gn)
    # Move back to the initial folder
    os.chdir("../")

def first_group_refinement(experiment):
    """Start all ensembles for the Design Groups"""
    # Start sharing (for file creation reasons)
    SHARING.Start(experiment)
    # Loop through the design groups
    print "Check1"
    for group in experiment:
        gn = group.number
        # Generate the name for the Design Group's calculations
        folder = "Group" + str(gn)
        try:
            os.mkdir(folder)
        except OSError:
            pass
        os.chdir(folder)
        # Start each ensemble
        for en in range(1, experiment["Ensemble Number"] + 1):
            eFolder = "Ensemble" + str(en)
            do = SHARING.claim_calculations(eFolder)
            # If this processor has started the calculations for this ensemble
            print en, do
            if do:
                # Make sure the experiment has the right structures and energies
                SHARING.update_Current(experiment, "../Current/", gn)
                SHARING.update_Energies(experiment, "../Current/", gn)
                # Move into the folder
                os.chdir(eFolder)
                # Make a copy of the best structures and energies for the group
                os.mkdir("Current")
                SHARING.output_Current(experiment, "./Current/", gn, 0)
                os.mkdir("Best")
                SHARING.output_Current(experiment, "./Best/", gn)
                SHARING.output_Energies(experiment, "./Best/", gn)
                # Move out of the e Folder
                os.chdir("../")
                SHARING.End(experiment)
                # Do the ensemble refinement
                ensemble_refinement(experiment, gn, en)
                SHARING.Start(experiment)
        # Move out of the folder
        os.chdir("../")
    SHARING.End(experiment)

def second_group_refinement(experiment):
    """Finish all ensembles of all Design Groups"""
    # Loop through the groups
    for group in experiment:
        gn = group.number
        # Move into the folder
        os.chdir("Group" + str(gn))
        # Loop through the ensembles
        for en in range(1, experiment["Ensemble Number"] + 1):
            ensemble_refinement(experiment, gn, en)
        # Move out of the folder
        os.chdir("../")

def finish_check(experiment, mn):
    """Check that all ensembles are entirely finished"""
    # NOTE THAT THIS DOES NOT END SHARING
    SHARING.Start(experiment)
    # Figure out what folder we should be looking for
    if experiment["Type"] == "Mutator":
        if mn in [None, 0]:
            folder = "wildtype"
        else:
            folder = "mutant" + str(mn)
    else:
        folder = "refinement"
    # Move to the Experiment's folder
    os.chdir(experiment["Folder"])
    # If the appropriate folder no longer exists, be done but return False so
    # the Finish function isn't called
    if folder not in os.listdir("./"):
        return False
    os.chdir(folder)
    finished = True
    # Loop through the groups
    for group in experiment:
        os.chdir("Group" + str(group.number))
        # and the ensembles
        for en in range(1, experiment["Ensemble Number"] + 1):
            os.chdir("Ensemble" + str(en))
            # Get the number of the last completed iteration
            I = SHARING.iteration_counter(SHARING.get_current(), False)
            # If they aren't all done yet
            if I < experiment["Refinement Iterations"]:
                finished = False
                os.chdir("../../")
                break
            os.chdir("../")
        if not finished:
            break
        os.chdir("../")
    return finished

def make_extra_group(experiment):
    """Create an extra Design Group with no Target Molecules"""
    # Only do this if Binding energy calculations are being done
    if experiment["Energy Calculation"] == "Binding":
        # Get all of the Design Molecules from Design Group 1
        molecules = []
        for molecule in experiment[1]:
            if molecule.design:
                molecules.append(molecule)
        # Make and store a new Design Group (which duplicates these Molecules
        # for run independence reasons)
        N = len(experiment) + 1
        group = MOLECULES.DesignGroup(N, molecules, experiment["Force Field"], \
                                      experiment["File Format"])
        experiment._groupOrder.append(N)
        experiment._groups[N] = group
        n = len(experiment)
        # Update the Current dictionary, too
        experiment["Current"][n] = {}
        for molecule in experiment[n]:
            new = molecule.duplicate()
            experiment["Current"][n][new.name] = new
        # Those structures are essentially place holders at this point. However,
        # we do need to run an energy minimization and calculate an initial
        # energy for that group
        name = "Group" + str(n) + "_Energies.txt"
        SHARING.Start(experiment)
        # If another processor already did the calculation, we're fine
        if name not in os.listdir("./Current/"):
            # Try to make a temp directory to do the calculations in
            if "temp" not in os.listdir("./"):
                # Make the directory and move into it
                os.mkdir("temp")
                os.chdir("temp")
                # Stop sharing
                SHARING.End(experiment)
                # Copy in the relevant files
                SHARING.copy_standard_files(experiment, False)
                # Relax the Design Group
                refinement = IPRO_FUNCTIONS.Relaxation(experiment, n, True)
                # And calculate the energies
                energies, refinement = \
                          IPRO_FUNCTIONS.Calculate_Energy(experiment, n)
                # Move back up a folder
                os.chdir("../")
                # Start sharing
                SHARING.Start(experiment)
                # Store the energy
                text = SHARING.format_energies(energies[n])
                f = open("./Current/" + name, "w")
                f.write(text)
                f.close()
                SHARING.output_Current(experiment, "./Current/", n)
                # Delete the temp folder
                os.system("rm -rf temp")
            # Otherwise, just wait
            else:
                SHARING.End(experiment)
                SHARING.Wait("temp", "./")
                SHARING.Start(experiment)
        # Store that complex energy
        f = open("./Current/" + name, "r")
        experiment["Energies"][n] = {"Complex":float(f.readline().split()[2])}
        f.close()
        SHARING.End(experiment)

def delete_extra_group(experiment):
    """Delete the extra design group made for binding calculations"""
    if experiment["Energy Calculation"] == "Binding":
        N = len(experiment)
        del experiment._groupOrder[N-1]
        del experiment._groups[N]
        del experiment["Current"][N]
        if N in experiment["Energies"]:
            del experiment["Energies"][N]
        n = len(experiment)

def Start(experiment, iteration, gn = None):
    """Create a folder to run a refinement in."""
    # Only do this when appropriate
    if not experiment["Do Refinement"] or gn != None or experiment["Activity"] \
    != "Standard":
        return False
    # Try to make the directory
    if experiment["Type"] == "Mutator":
        folder = "wildtype"
    else:
        folder = "refinement"
    try:
        os.mkdir(folder)
    except OSError:
        return False
    # Assuming that was successful, prep that folder
    folder += "/Current/"
    os.mkdir(folder)
    SHARING.output_Current(experiment, folder)
    SHARING.output_Energies(experiment, folder)
    f = open(folder + "iteration.txt", "w")
    f.write(str(iteration))
    f.close()
    return True

def change_iterations(experiment, mn = None):
    """Do 2x the calculations for initial / wildtype calcs"""
    # If there is a mutant number that's not 0, don't do this
    if mn not in [None, 0]:
        return False
    # If it is None, check the number of completed iterations
    if mn == None:
        n = SHARING.iteration_counter(SHARING.get_current(), False)
        if n > 0:
            return False
    # multiply the number by 2
    experiment["Refinement Iterations"] *= 2
    return True

def initialize(experiment, mn):
    """Set up an Experiment to work in a refinement"""
    # Start sharing
    SHARING.Start(experiment)
    # Modify the number of iterations
    modify = change_iterations(experiment, mn)
    # Identify the proper folder for the refinement
    if experiment["Type"] == "Mutator":
        if mn in [None, 0]:
            folder = "./wildtype/"
        else:
            folder = "./mutant" + str(mn) + "/"
    else:
        folder = "./refinement/"
    # Move into the refinement folder
    os.chdir(folder)
    # Load the current structures and energies
    SHARING.update_Current(experiment, "./Current/")
    SHARING.update_Energies(experiment, "./Current/")
    # Stop sharing
    SHARING.End(experiment)
    # Make the extra Design Group
    make_extra_group(experiment)
    # Modify the 'Activity' of the Experiment
    experiment['Activity'] = "Refinement"
    return modify

def gather_energies(experiment):
    """Collect all calculated energies"""
    # Store them here
    energies = {}
    # Go through the Design Groups
    for group in experiment:
        # initialize the energies dictionary
        energies[group.number] = {}
        # Go through the ensembles
        for en in range(1, experiment["Ensemble Number"] + 1):
            energies[group.number][en] = {}
            # Access the appropriate file
            name = "Group"+str(group.number)+"/Ensemble"+str(en)+"/Best/"
            name += "Group" + str(group.number) + "_Energies.txt"
            f = open(name, "r")
            # Store the energies
            for line in f:
                items = line.split(":")
                if len(items) != 2:
                    continue
                type = items[0].split()[0].capitalize()
                energies[group.number][en][type] = float(items[1])
            f.close()
    return energies

def calculate_average(values):
    """Calculate the average and standard deviations of a list of values"""
    # Store the calculated average here
    average = 0
    for value in values:
        average += value
    average /= float(len(values))
    # And the calculated standard deviation here
    sd = 0
    for value in values:
        sd += math.pow(average - value, 2)
    sd = math.sqrt(sd/float(len(values)))
    return average, sd

def calculate_IEs(energies):
    """Calculate interaction energy information"""
    # Store the calculated values here
    IEs = {}
    # Go through each Design Group
    for gn in energies:
        # If there isn't interaction energies
        if "Interaction" not in energies[gn][1]:
            continue
        # Make a list of the interaction energies from each Design Group
        values = []
        for en in energies[gn]:
            values.append(energies[gn][en]["Interaction"])
        # Calculate the average and standard deviation of the interaction
        # energies
        average, sd = calculate_average(values)
        # Store that info for the Design Group
        IEs[gn] = {"Average":[average, sd]}
        # Find the Design Group that is closest to the average value
        best = None
        difference = 5000000
        for en in energies[gn]:
            value = math.fabs(average - energies[gn][en]["Interaction"])
            if value < difference:
                difference = value
                best = en
        # Store that best choice
        IEs[gn]["Choice"] = [best, energies[gn][en]]
    return IEs

def calculate_CEs(energies):
    """Calculate the average complex energies of the Design Groups"""
    # This is pretty similar to the IE function, so I'm not going to comment
    CEs = {}
    for gn in energies:
        values = []
        for en in energies[gn]:
            values.append(energies[gn][en]["Complex"])
        average, sd = calculate_average(values)
        CEs[gn] = {"Average":[average, sd]}
    return CEs

def calculate_BEs(experiment, energies):
    """Calculate Binding Energies"""
    # Store the energies here
    BEs = {}
    # Only do this if there are binding energies to calculate
    if experiment["Energy Calculation"] == "Binding":
        # Collect the Complex energies from the last Design Group
        N = len(experiment)
        designs = []
        for en in energies[N]:
            designs.append(energies[N][en]["Complex"])
        # Go through the other Design Groups
        for gn in energies:
            # Skip the last group
            if gn == N:
                continue
            # Calculate the different energies for each ensemble
            ensembles = {}
            allValues = []
            for en in energies[gn]:
                values = []
                for value in designs:
                    values.append(energies[gn][en]["Complex"] - value - \
                                  experiment["Target Energies"][gn])
                # Calculate this average and SD
                average, sd = calculate_average(values)
                # Store the average
                ensembles[en] = average
                # And all of the values
                allValues.extend(values)
            # Calculate the overall average
            average, sd = calculate_average(allValues)
            # Store that
            BEs[gn] = {"Average":[average, sd]}
            # Determine which ensemble is closest to that
            best = None
            difference = 5000000
            for en in energies[gn]:
                value = math.fabs(ensembles[en] - average)
                if value < difference:
                    difference = value
                    best = en
            # Store the best ensemble choice
            BEs[gn]["Choice"] = [best, ensembles[en]]
    return BEs

def make_refinement_choice(experiment, IEs, BEs):
    """Make a decision regarding the results of a structure refinement"""
    # Always keep the results of mutator experiments
    if experiment["Type"] == "Mutator":
        return True
    # Or if this is the initial refinement
    f = open("./Current/iteration.txt", "r")
    n = int(f.readline())
    f.close()
    if n == 0:
        return True
    # Otherwise, make sure the Experiment is using its reference energies
    SHARING.load_reference_Energies(experiment)
    # Get the type of energy calculation
    type = experiment["Energy Calculation"]
    # Go through the Design groups
    for group in experiment:
        if type == "Binding" and group.number == len(experiment):
            continue
        # Get the energy to compare
        if type == "Binding":
            energy = BEs[group.number]["Average"][0]
        else:
            energy = IEs[group.number]["Average"][0]
        energy -= experiment["Energies"][group.number][type]
        # Make the choice based on the objective
        if group.objective in ['improve', 'maintain']:
            if energy > 0:
                return False
        elif group.objective in ['reduce', 'eliminate']:
            if energy < 0:
                return False
    # if all Design Groups passed the test, return True
    return True

def store_results(experiment, IEs, BEs, CEs, mn = None):
    """Store the results of the refinement and share them"""
    # If this isn't a mutator, get the current iteration
    if experiment["Type"] != "Mutator":
        f = open("./Current/iteration.txt", "r")
        mn = int(f.readline())
        f.close()
    # Go through the Design Groups, identify the chosen ensemble, and store the
    # structures and relevant energies
    for group in experiment:
        if group.number == len(experiment) and BEs != {}:
            continue
        # Store the energies
        experiment["Energies"][group.number] = {}
        experiment["Energies"][group.number]["Complex"] = \
                CEs[group.number]["Average"][0]
        experiment["Energies"][group.number]["Interaction"] = \
                IEs[group.number]["Average"][0]
        if BEs != {}:
            experiment["Energies"][group.number]["Binding"] = \
                BEs[group.number]["Average"][0]
            # Identify the best ensemble structure
            en = BEs[group.number]["Choice"][0]
        # Otherwise get the best ensemble structure from the IEs
        else:
            en = IEs[group.number]["Choice"][0]
        # Load the structures from that ensemble
        folder = "Group" + str(group.number) + "/Ensemble" + str(en) + "/Best/"
        SHARING.update_Current(experiment, folder, group.number)
        # If this isn't a mutator, share them in the Experiment's folder
        if experiment["Type"] != "Mutator":
            # Share the structures in the Experiment's Current folder
            SHARING.output_Current(experiment, experiment["Folder"]+"Current/",\
                                   group.number, mn)
        # Share the results
        SHARING.output_best(experiment, mn, group.number)

def Summarize(experiment, choice, IEs, BEs, CEs, mn = None):
    """Summarize the results of a structure refinement"""
    # Come up with an appropriate header
    if experiment["Type"] == "Mutator":
        if mn == 0:
            experiment["Summary"] = "\nWildtype Refinement Results\n"
        else:
            experiment["Summary"] = "\nMutant " + str(mn)+" Refinement Results\n"
    else:
        f = open("./Current/iteration.txt", "r")
        n = int(f.readline())
        f.close()
        if n == 0:
            experiment["Summary"] = "\nInitial Structures Refinement\n"
        else:
            experiment["Summary"] = "\nRefinement of Iteration " + str(n) + "\n"
    # Find out when the refinement started
    f = open("./Group1/Ensemble1/Ensemble1_Summary.txt", "r")
    for line in f:
        if line.startswith("Started"):
            experiment["Summary"] += line
            break
    f.close()
    # Store the information about the calculated energies
    for group in experiment:
        if experiment["Energy Calculation"] == "Binding" and group.number == \
        len(experiment):
            continue
        # Summarize Complex energy
        experiment["Summary"] += "Design Group " + str(group.number) + "\n"
        experiment["Summary"] += "Average Complex Energy: " + \
            format(CEs[group.number]["Average"][0], '.3f') + " +/- " + \
            format(CEs[group.number]["Average"][1], '.3f') + " kcal/mol\n"
        # IE
        experiment["Summary"] += "Average Interaction Energy: " + \
            format(IEs[group.number]["Average"][0], '.3f') + " +/- " + \
            format(IEs[group.number]["Average"][1], '.3f') + " kcal/mol\n"
        # BE
        if experiment["Energy Calculation"] == "Binding":
            experiment["Summary"] += "Average Binding Energy: " + \
            format(BEs[group.number]["Average"][0], '.3f') + " +/- " + \
            format(BEs[group.number]["Average"][1], '.3f') + " kcal/mol\n"
    if experiment["Type"] != "Mutator" and n != 0:
        if choice:
            experiment["Summary"] += "These were the BEST results so far\n"
        else:
            experiment["Summary"] += "These results were DISCARDED\n"
    # Say when the calculations finished
    experiment["Summary"] += "Ended" + SHARING.time_stamp()

def Finish(experiment, mn = None):
    """Finish a Refinement"""
    # Sharing will have been started elsewhere (in finish_check)
    # Collect the calculated energies
    energies = gather_energies(experiment)
    # Get interaction, complex, and binding energies
    IEs = calculate_IEs(energies)
    CEs = calculate_CEs(energies)
    BEs = calculate_BEs(experiment, energies)
    # Evaluate the results
    choice = make_refinement_choice(experiment, IEs, BEs)
    # If the choice is to keep the results, make a permanent copy
    if choice:
        store_results(experiment, IEs, BEs, CEs, mn)
    # Summarize the results
    Summarize(experiment, choice, IEs, BEs, CEs, mn)
    # Move out of the refinement's folder
    os.chdir(experiment["Folder"])
    # Store the refinement summary
    name = SHARING.summary_name(experiment["Folder"])
    f = open(name, "a")
    f.write(experiment["Summary"])
    f.close()
    # Delete the folder the refinement happened in
    if experiment["Type"] == "Mutator":
        if mn in [None, 0]:
            name = "wildtype"
        else:
            name = "mutant" + str(mn)
    else:
        name = "refinement"
    os.system("rm -rf " + name)
    # End sharing
    SHARING.End(experiment)

def DO(experiment, mn = None):
    """Do a Refinement"""
    # Only do a refinement when there is a refinement folder to do the
    # refinement IN
    if experiment["Type"] == "Mutator":
        if mn in [None, 0]:
            name = "wildtype"
        else:
            name = "mutant" + str(mn)
    else:
        name = "refinement"
    # If the folder doesn't exist, don't do anything
    if name in os.listdir("./"):
        # Initialize the experiment for the refinement
        modify = initialize(experiment, mn)
        # do the first set of refinement calculations - which exist to try and
        # keep the processors working on separate jobs as much as possible
        first_group_refinement(experiment)
        # And the second, which exists to get everything finished ASAP
        second_group_refinement(experiment)
        # Determine if the calculations are finished yet
        finished = finish_check(experiment, mn)
        # If the refinement is not finished, move out of the refinement's folder
        # and wait for it to be deleted
        if not finished:
            os.chdir(experiment["Folder"])
            SHARING.End(experiment)
            # There's no need to wait during mutator experiments, except for the
            # Wildtype values to finish
            if experiment["Type"] != "Mutator" or mn in [None, 0]:
                SHARING.Wait(name)
        # If the refinement is finished, store the information
        else:
            Finish(experiment, mn)
        # The refinement folder is gone, the structures and energies are stored,
        # so we can be done with the refinement after a few more things
        delete_extra_group(experiment)
        if modify:
            experiment["Refinement Iterations"] /= 2
        experiment["Activity"] = "Standard"
