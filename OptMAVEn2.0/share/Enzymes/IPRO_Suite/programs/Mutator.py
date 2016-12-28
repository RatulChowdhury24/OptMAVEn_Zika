#!/usr/bin/env python

# Where the IPRO Suite is installed
InstallFolder = "/gpfs/group/c/cdm/IPRO_Suite/"

# name of this file
__name__ = "IPRO Suite Program to Make Mutants"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This program creates and evaluates user-specified mutants of proteins."""

# Load the PYTHON modules
import os
import sys
sys.path.append(InstallFolder + "modules/")
import EXPERIMENT
import IPRO_FUNCTIONS
import REFINEMENT
import MUTATOR

# Load the experiment
experiment = EXPERIMENT.Experiment()
# Initialize it
dummy = IPRO_FUNCTIONS.INITIALIZE(experiment)
# Do a Refinement of those structures
REFINEMENT.DO(experiment, 0)

# Make each mutant
for mn in range(1, len(experiment["Mutants"]) + 1):
    # Make the initial structures of the Mutant
    MUTATOR.DO(experiment, mn)
    # And then refine them
    REFINEMENT.DO(experiment, mn)
