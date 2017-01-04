#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Functions for File Sharing"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for the sharing of information between different
processors during IPRO Suite Experiments."""


# Standard PYTHON modules
import os
import shutil #mfa#
import sys
import time
import datetime
# Include all contents from the STANDARDS file
from STANDARDS import *
import REFINEMENT

class SharingError(IPRO_Error):
    """An error for problems sharing files between processors"""
    def __init__(self, error = ''):
        """The initialization of the SharingError class"""
        IPRO_Error.__init__(self, error)

def Start(experiment):
    """Create a folder that says that a processor is sharing files"""
    # Know the time everything started
    startTime = time.time()
    did = False
    # Do this for 15 minutes
    print "=" * 60
    while time.time() - startTime < 900:
        # Make a folder called 'sharing'
        try:
            os.mkdir(os.path.join(experiment["Folder"], "sharing")) #mfa#
            # If successful, break the while loop
            did = True
            break
        except NameError:
            # Take a 1 second break before trying again
            print os.listdir(".")
            time.sleep(1)
    # If this never succeeded, raise an error
    if not did:
        text = "The processor waited more than 15 minutes to start sharing "
        text += "files. This indicates a major problem."
        raise SharingError(text)

def End(experiment):
    """Delete the sharing folder"""
    # This folder
    folderName = experiment["Folder"] + "sharing"
    # Get rid of it
    try:
        print "deleted"
        os.rmdir(folderName)
    except OSError:
        print "deleted2"
        # Its weird if that failed. Do a hard delete of the folder
        if "sharing" in os.listdir(experiment["Folder"]):
            os.system('rm -rf ' + folderName)

def get_current():
    """Get the current directory.

    This function exists because I've run into problems with os.getcwd()"""
    # Use a while loop
    accept = False
    while not accept:
        try:
            current = os.getcwd() + "/"
            accept = True
        except OSError:
            pass
    return current

def time_stamp():
    """Declare the date and time"""
    # Use datetime
    a = datetime.datetime.now()
    # Store the information here
    text = " on " + str(a.month) + "/" + str(a.day) + "/" + str(a.year) + " at "
    text += str(a.hour).rjust(2, '0') + ":" + str(a.minute).rjust(2, '0') + ":"
    text += str(a.second).rjust(2, '0') + "\n"
    return text

def summary_update(startTime, procedure, extra = ''):
    """Create a string saying how long a procedure took"""
    # Calculate the number of seconds
    s = int(time.time() - startTime)
    # minutes and seconds
    m = s/60
    if m > 0:
        s = s % (m*60)
    # Here's the text
    text = procedure + " took "
    if m > 0:
        text += str(m) + " minutes and "
    text += str(s) + " seconds\n" + extra
    return text

def summary_name(folder):
    """Generate the name of a summary file based on the folder"""
    # Split the folder on /
    folders = folder.split("/")
    # Do this differently depending on whether or not there was a slash at the
    # end
    if folders[-1] == '':
        name = folder + folders[-2] + "_Summary.txt"
    else:
        name = folder + "/" + folders[-1] + "_Summary.txt"
    return name

def iteration_counter(folder, ongoing = True):
    """count completed and ongoing iterations of IPRO"""
    # Store the number of completed iterations
    iteration = 0
    # Get the name of the file
    fileName = summary_name(folder)
    # Try to read it 
    try:
        f = open(fileName, "r")
        for line in f:
            if line.startswith("Iteration"):
                iteration = int(line.split()[1])
        f.close()
    # If you can't, that's OK. It may not exist yet
    except IOError:
        pass
    # If the number of ongoing iterations should be counted
    if ongoing:
        names = os.listdir(folder)
        for name in names:
            if name.startswith("iteration"):
                iteration += 1
    return iteration

def claim_calculations(folder):
    """Try to create the specified folder. Return success or failure"""
    try:
        os.mkdir(folder)
        return True
    except OSError:
        return False

def copy_cpp_files(rotamer = True, docking = True, folder = "./"):
    """Copy the CPP files from the installation folder to the current one."""
    # Store the names of the files here
    files = []
    # If the rotamer files are needed
    if rotamer:
        files.extend(['rotamer_constant.out', 'rotamer_rotamer.out'])
    # If the docking program is needed
    if docking:
        files.append("docking.out")
    # copy each file from the installation folder
    for fileName in files:
        shutil.copy(os.path.join(InstallFolder, "modules", "CPP", fileName),
                os.getcwd()) #mfa#
    
def copy_forceField_files(experiment, folder):
    """Copy the force field files to the current folder."""
    # Store the file names here
    files = []
    if experiment["Force Field"] == "CHARMM":
        files.extend(experiment["CHARMM Topology Files"])
        files.extend(experiment["CHARMM Parameter Files"])
    else:   
        text = "The copy_forceField_files function does not support the "
        text += str(experiment["Force Field"]) + " force field."
        raise IPRO_IO_Error(text)
    # Copy each file
    for fileName in files:
        shutil.copy(os.path.join(InstallFolder, "input_files", fileName),
                os.getcwd()) #mfa#
        
def copy_solvation_files(experiment, folder):
    """Copy the implicit solvation files to the current folder."""
    # Store the file names here
    files = []
    # Go through the different types of possible file lists (currently only LK)
    for attribute in ["LK Solvation Files"]:
        if attribute in experiment:
            files.extend(experiment[attribute])
    # Copy each file
    for fileName in files:
        shutil.copy(os.path.join(InstallFolder, "input_files", fileName),
                os.getcwd()) #mfa#

def copy_human_sequences_files(experiment, folder):
    """ Copy the human 9mer sequences file to the current folder. """
    # Store the file names here
    files = []
    # Go through the different types of possible file lists
    for attribute in ["Human Sequences Files"]:
        if attribute in experiment:
            files.extend(experiment[attribute])
    # Copy each file
    for fileName in files:
        command = "cp " + folder + fileName + " ."
        i = os.system(command)
        if i != 0:
            text = "Copying " + fileName + " has failed."
            raise SharingError(text)

def copy_integer_cuts_files(experiment, folder):
    """ Copy the MAPs integer cuts file to the current folder. """
    # Store the file names here
    files = []
    # Go through the different types of possible file lists (currently only LK)
    for attribute in ["Maps Integer Cuts Files"]:
        if attribute in experiment:
            files.extend(experiment[attribute])
    # Copy each file
    for fileName in files:
        command = "cp " + folder + fileName + " ."
        i = os.system(command)
        if i != 0:
            text = "Copying " + fileName + " has failed."
            raise SharingError(text)

def copy_antibody_parts_files(experiment, folder):
    """ Copy the antibody parts files to the current folder. """
    # Store the file names here
    files = []
    # Go through the different types of possible file lists (currently only LK)
    for attribute in ["Antibody Parts Files"]:
        if attribute in experiment:
            files.extend(experiment[attribute])
    # Copy each file
    for fileName in files:
        command = "cp " + folder + fileName + " ."
        i = os.system(command)
        if i != 0:
            text = "Copying " + fileName + " has failed."
            raise SharingError(text)

def copy_optmaven_cpp_files():
    """Copy the OptMAVEn CPP files from the installation folder to the current one."""
    # Store the names of the files here
    files = []
    # If the rotamer files are needed
    files.extend()
    if rotamer:
        files.extend(['rotamer_constant.out', 'rotamer_rotamer.out'])
    # If the docking program is needed
    if docking:
        files.append("docking.out")
    # copy each file from the installation folder
    for fileName in files:
        # Generate the command to copy the file
        command = "cp " + InstallFolder + "modules/CPP/" + fileName
        command += " ." 
        # Copy it
        i = os.system(command)
        # If that didn't work, raise an error
        if i != 0:
            text = "Copying " + fileName + " has failed."
            raise SharingError(text)

def copy_standard_files(experiment, folder = None, cpp = True, ff = True, solv=\
                        False):
    """Copy the standard IPRO Suite files to a folder."""
    # If no folder has been specified, use the Experiment's input files folder
    if not isinstance(folder, str):
        folder = os.path.join(experiment["Folder"], "input_files") #mfa#
    # Copy the expected files
    if cpp:
        copy_cpp_files(True, True, experiment["Folder"])
    if ff:
        copy_forceField_files(experiment, folder)
    # The default is to NOT copy the solvation files, because the experiment
    # really doesn't need them often.
    if solv:
        copy_solvation_files(experiment, folder)

def output_Current(experiment, folder, gn = None, last = None):
    """Share the contents of the Current dictionary"""
    # Loop through the relevant Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Output every Molecule
        for molecule in group:
            name = folder + molecule.generate_name()
            experiment["Current"][group.number][molecule.name].output(name, \
                                     molecule.fileFormat, experiment["User"])
    # If appropriate, make a 'last update' file
    if last != None:
        f = open(folder + "last_update.txt", "w")
        f.write(str(last))
        f.close()

def identify_last_update(folder):
    """Identify the last time structures and energies were updated"""
    last = 0
    try:
        f = open(folder + "last_update.txt", "r")
        last = int(f.readline())
        f.close()
    except IOError:
        pass
    return last

def format_energies(energies, gn = None, simple = True):
    """Create text listing energy values"""
    # A list of expected energy terms
    terms = ['Complex', 'Interaction', 'Binding']
    # Include any other terms in alphabetical order
    all = energies.keys()
    all.sort()
    for term in all:
        if term not in terms:
            terms.append(term)
    # Create the formatted text
    text = ''
    for term in terms:
        # Skip irrelevant terms
        if term not in energies:
            continue
        # Format the energy
        energy = format(energies[term], '.3f')
        # Create the text differently, depending on method
        if simple or gn == None:
            label = term + " Energy:"
            text += label.ljust(20) + energy + "\n"
        else:
            text += "The " + term + " Energy of Design Group " + str(gn) +" is "
            text += energy + " kcal / mol\n"
    return text

def output_Energies(experiment, folder, gn = None):
    """Output the contents of the Energies dictionary"""
    # Loop through relevant Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # The file for the energies
        name = folder + "Group" + str(group.number) + "_Energies.txt"
        # The text
        text = format_energies(experiment["Energies"][group.number])
        f = open(name, "w")
        f.write(text)
        f.close()

def update_Current(experiment, folder, gn = None):
    """Update the Current structures of the Experiment"""
    # Loop through the relevant Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Go through every Molecule
        for molecule in group:
            # Identify the name of the file
            name = folder + molecule.generate_name()
            # Load those contents
            f = open(name, "r")
            lines = f.readlines()
            f.close()
            experiment["Current"][group.number][molecule.name].load(lines)

def update_Energies(experiment, folder, gn = None):
    """Update the Energies dictionary of the Experiment"""
    # Loop through the relevant Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Load the file that contains the energies
        f = open(folder + "Group" + str(group.number) + "_Energies.txt", "r")
        # Reset the entry for this Design Group
        experiment["Energies"][group.number] = {}
        for line in f:
            # Split on the expected :
            items = line.split(":")
            if len(items) != 2:
                continue
            # Get the type and the energy
            type = items[0].split()[0].capitalize()
            energy = float(items[1])
            experiment["Energies"][group.number][type] = energy
        f.close()

def load_Target_Energies(experiment, folder):
    """Load the Target Molecule energies"""
    # Only do this if Binding energy calculations are being done
    if experiment["Energy Calculation"] == "Binding":
        # Store the energies of the Target Molecules here
        targets = {}
        # And the sums of their energies here
        experiment["Target Energies"] = {}
        # Go through each Design Group
        for group in experiment:
            # The total energy for each of this Group's Target Molecules
            energy = 0
            # Find each Target Molecule
            for molecule in group:
                if not molecule.design:
                    mn = molecule.name
                    # if the energy hasn't been found yet, look it up
                    if mn not in targets:
                        try:
                            name = "Molecule" + mn + "_Energy.txt"
                            f = open(folder + name, "r")
                            targets[mn] = float(f.readline().split()[-1])
                            f.close()
                        except IOError:
                            text = "Could not find " + name 
                            raise SharingError(text)
                    # modify the energy sum
                    energy += targets[mn]
            # Store the energy for the Design Group
            experiment["Target Energies"][group.number] = energy

def load_reference_Energies(experiment):
    """Load the energies of the Best results so far"""
    # Determine where the best results are located
    iteration = 0
    folders = os.listdir(os.path.join(experiment["Folder"], "results")) #mfa#
    for folder in folders:
        if folder.startswith("iteration"):
            n = int(folder[9:])
            if n > iteration:
                iteration = n
    # The folder where initial results should be stored
    if experiment["Type"] == "Mutator":
        initial = experiment["Folder"] + "results/wildtype/"
    else:
        initial = experiment["Folder"] + "results/initial/"
    # Go through the Design Groups
    for group in experiment:
        # If this is during a refinement and binding energy calculations are
        # being used, don't look for reference energies for the last Design
        # Group. They're meaningless and useless.
        if experiment["Activity"] == "Refinement" and \
        experiment["Energy Calculation"] == "Binding" and group.number == \
        len(experiment):
            continue
        # If the initial energies should be used for comparison
        if iteration == 0 or group.objective in ['maintain', 'reduce']:
            folder = initial
        else:
            folder = experiment["Folder"] + "results/iteration" + str(iteration)
            folder += "/"
        # Update the Design Group's energies from that folder
        update_Energies(experiment, folder, group.number)

def output_best(experiment, iteration, gn = None):
    """Permanently store the Experiment's contents in results"""
    folder = experiment["Folder"] + "results/"
    # Create the appropriate folder name
    if experiment["Type"] == "Mutator":
        if iteration == 0:
            folder += "wildtype/"
        else:
            folder += "mutant" + str(iteration) + "/"
    else:
        if iteration == 0:
            folder += "initial/"
        else:
            folder += "iteration" + str(iteration) + "/"
    # Make that folder
    try:
        os.mkdir(folder)
    # If it already exists, that can be ignored
    except OSError:
        pass
    # Write the Experiment's Current and Energies dictionaries to that folder
    output_Current(experiment, folder, gn)
    output_Energies(experiment, folder, gn)

def Results(experiment, iteration, best, keep, gn = None):
    """Share results between processors if that is appropriate"""
    # First, dump the contents of the experiment to the Current folder
    if best or (experiment["Annealing Sharing"] and keep):
        output_Current(experiment, "./Current/", gn, iteration)
    # If these are the best results so far, more needs to be done
    if best:
        # If this is not during a standard portion of IPRO
        if experiment["Activity"] != "Standard":
            output_Current(experiment, "./Best/", gn)
            output_Energies(experiment, "./Best/", gn)
        # If a refinement should be started, do so
        elif experiment["Do Refinement"]:
            did = REFINEMENT.Start(experiment, iteration, gn)
        # Otherwise, output these results as the best so far
        else:
            output_best(experiment, iteration, gn)

def Wait(folder1, folder2 = "./", seconds = 7200):
    """Wait for up to 2 hours for a folder to be deleted"""
    # Know when this started
    startTime = time.time()
    # Use a flag
    gone = False
    while time.time() - startTime <= seconds:
        if folder1 not in os.listdir(folder2):
            gone = True
            break
        else:
            # Take a 1 second break
            time.sleep(1)
    # If this has taken too long, throw an error
    if not gone:
        os.system("rm -r /gpfs/scratch/mjg5185/requested")
        os.mkdir("/gpfs/scratch/mjg5185/requested")
        os.system("cp -r " + folder2 + " /gpfs/scratch/mjg5185/requested") 
        text = "The processor waited over " + str(seconds) + " seconds for the "
        text += folder1 + " folder to be deleted. This indicates a major issue."
        raise SharingError(text)
