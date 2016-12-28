#!/usr/bin/env python

InstallFolder = "/gpfs/scratch/aml52/IPRO_Suite/"

# The name of this program
__name__ = "The IPRO Suite Experiment Starter"

# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This program creates an Experiment class object, which asks the user for all
information needed to run an experiment. Once that has been gathered, it goes to
the Experiment's folder, puts the proper program there, and optionally starts
the Experiment."""

# Include PYTHON information
import os
import sys
sys.path.append(InstallFolder + "modules/")
# Include the EXPERIMENT module
import EXPERIMENT
# The STANDARDS module is needed for answering questions
from STANDARDS import *
# The source of the functions that can generate scripts to run the experiment
import SUBMITTER

# Make the experiment, using user input
experiment = EXPERIMENT.Experiment(False)
# Move to the experiment's folder and output it's information
os.chdir(experiment["Folder"])
experiment.output()
# Copy the appropriate program to this folder
command = "cp " + InstallFolder + "programs/" + experiment["Type"] + ".py ."
i = os.system(command)
# Generate scripts to run this experiment and possibly submit them
dummy = SUBMITTER.script_generator(experiment)
