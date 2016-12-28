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
#FIXME OPTMAVEN2.initialize_antigen(experiment)

# Determine which positions do not cause clashes.
#FIXME OPTMAVEN2.cull_clashes(experiment)

# Calculate the interaction energies with the MAPs database.
OPTMAVEN2.interaction_energies(experiment)
