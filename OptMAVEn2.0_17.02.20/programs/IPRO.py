#!/usr/bin/env python

# Where the IPRO Suite is installed
InstallFolder = "/gpfs/work/m/mfa5147/OptMAVEn2.0"

# Name of this file
__name__ = "The Iterative Protein Redesign & Optimization Program"

# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This program runs IPRO to design proteins to have modified binding behaviors
(typically improved) to target molecules."""

# Load the needed PYTHON Modules
import os
import sys
sys.path.append(InstallFolder + "modules/")
import EXPERIMENT
import IPRO_FUNCTIONS
import REFINEMENT

# Start by loading an Experiment
experiment = EXPERIMENT.Experiment()
# Initialize the Experiment so that it can run calculations
dummy = IPRO_FUNCTIONS.INITIALIZE(experiment)
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
