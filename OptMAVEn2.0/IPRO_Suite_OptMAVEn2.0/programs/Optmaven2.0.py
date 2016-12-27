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
(typically improved) to target molecules."""

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

# Start by loading an Experiment
experiment = EXPERIMENT.Experiment()
# Ensure that the experiment is an OptMAVEn experiment.
if experiment["Type"] != "OptMAVEn":
    raise TypeError('Cannot run OptMAVEn on an experiment of type "{}."'
            .format(experiment["Type"]))
# Initialize the Experiment so that it can run calculations
#mfa# dummy = IPRO_FUNCTIONS.INITIALIZE(experiment)

# Move the antigen to the initial position and generate a PSF.
OPTMAVEN2.initialize_antigen(experiment)

exit()
# Lay out the grid search.
GRID_LEVELS_FILE = "grid_levels.dat"
grid_dimensions = ["z angle", "x translation", "y translation", "z translation"]
grid_levels = list()
"""
for dimension in grid_dimensions:
	dmin = input("What is the minimum {}? ".format(dimension))
	dmax = input("What is the maximum {}? ".format(dimension))
	dstep = input("What is the {} step? ".format(dimension))
	grid_levels.append(np.arange(dmin, dmax, dstep))
"""
grid_levels = [[0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340], [-5.0, -3.75, -2.5, -1.25, 0.0, 1.25, 2.5, 3.75, 5.0, 6.25, 7.5, 8.75], [-10.0, -8.75, -7.5, -6.25, -5.0, -3.75, -2.5, -1.25, 0.0, 1.25, 2.5, 3.75], [3.0, 4.25, 5.5, 6.75, 8.0, 9.25, 10.5, 11.75, 13.0, 14.25, 15.5]]
# Determine the optimal way to delegate the grid points to various processes.
# The optimal way to delegate points will minimize the number of points given to each process while keeping the total number of processes under the max_processes limit.
# For each dimension, compute the maximum number of grid points that a process will need to compute at that dimension given that the level is divided a certain number of times.
max_points = [[int(math.ceil(float(len(dimension)) / num_divisions)) for num_divisions in range(1, 1 + len(dimension))] for dimension in grid_levels]
raw_input(max_points)

open(GRID_LEVELS_FILE, "w").write("\n".join([" ".join(map(str, dimension)) for dimension in grid_levels]))

# Initialize the antigen position and generate a PSF.
os.system("python test_initial_antigen_steps.py")

# Determine which positions do not cause clashes.
os.system("vmd test_vmd_clashes.tcl")

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
