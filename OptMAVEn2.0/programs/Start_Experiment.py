#!/usr/bin/env python

InstallFolder = "/home/matthew/Maranas_Lab/Mini_OptMAVEn2_Test"

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
sys.path.append(os.path.join(InstallFolder, "modules"))
# Include the EXPERIMENT module
import EXPERIMENT

# Make the experiment, using user input
experiment = EXPERIMENT.Experiment(False)
# Move to the experiment's folder and output it's information
os.chdir(experiment["Folder"])
experiment.output()
# Start running OptMAVEn from the experiment's folder.
os.system("python {}".format(os.path.join(InstallFolder, "programs",
        "OptMAVEn2.0.py")))

