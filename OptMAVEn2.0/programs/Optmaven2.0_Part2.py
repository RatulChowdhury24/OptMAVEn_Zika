#!/usr/bin/env python

# Where the IPRO Suite is installed
InstallFolder = "/home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0"

# Name of this file
__name__ = "The Iterative Protein Redesign & Optimization Program"

# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This program runs IPRO to design proteins to have modified binding behaviors
(typically improved) to target molecules.

This is the SECOND part of OptMAVEn. This part is run when the interaction
energy calculations, which occupy the bulk of the time spent running OptMAVEn,
finish.
"""

# Load the needed PYTHON Modules
import math
import os
import sys

import numpy as np

# Import the IPRO modules.
sys.path.append(os.path.join(InstallFolder, "modules"))
import EXPERIMENT
import IPRO_FUNCTIONS
import OPTMAVEN2
import REFINEMENT


# Get the experiment directory from the arguments.
if len(sys.argv) != 2:
    raise OSError("Usage: args: experiment/directory")
directory = sys.argv[1]
if not os.path.isdir(directory):
    raise OSError("Experiment directory does not exist: {}".format(directory))
os.chdir(directory)

# Load an experiment from the Experiment Details file.
experiment = EXPERIMENT.Experiment()
# Ensure that the experiment is an OptMAVEn experiment.
if experiment["Type"] != "OptMAVEn":
    raise TypeError('Cannot run OptMAVEn on an experiment of type "{}."'
            .format(experiment["Type"]))

# Check if all of the interaction energy calculations have finished.
energy_directory = os.path.join(directory, "energies")
for molecule in os.listdir(energy_directory):
    mol_directory = os.path.join(energy_directory, molecule)
    for part in os.listdir(mol_directory):
        part_directory = os.path.join(mol_directory, part)
        # If at least one energy calculation has not finished, abort and wait
        # until they are all finished.
        if not os.path.isfile(os.path.join(part_directory, "finished")):
            sys.exit(0)

# If all of the calculations have finished, then reorganize the energies, which are grouped first by part and second by position, to be grouped first by position and then by part.
OPTMAVEN2.regroup_energies(experiment)


# Pick the antibody parts combination for each position



# Refine the antigen position and antibody parts combination

# IPRO to redesign the selected antibody parts against antigen
# Declare an iteration variable and set it to 0
iteration = 0
# Do the appropriate number of iterations
while iteration <= experiment["IPRO Iterations"]:
    # If appropriate, do a refinement
    REFINEMENT.DO(experiment)
    # And if not, do an iteration of IPRO. It will prematurely terminate if
    # another processor calls for a refinement while it is running
    iteration = IPRO_FUNCTIONS.IPRO_ITERATION(experiment)
# Wait for IPRO to be finished, assisting with any Refinements that get started
# while waiting
IPRO_FUNCTIONS.Wait(experiment)
