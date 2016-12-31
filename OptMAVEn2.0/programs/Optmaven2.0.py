#!/usr/bin/env python

# Where the IPRO Suite is installed
InstallFolder = "/gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0"

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
directory = experiment["Folder"]

# Determine the status of the experiment.
""" Status codes:
None: The experiment has not yet started.
   1: The MAPs parts interaction energy calculations are being performed and may
      have finished.
   2: The MILP is being solved and may have finished.
"""
status = experiment.get_OptMAVEn_status()

if status is None:
    # Move the antigen to the initial position and generate a PSF.
    OPTMAVEN2.initialize_antigen(experiment)

    # Determine which positions do not cause clashes.
    OPTMAVEN2.cull_clashes(experiment)

    # Calculate the interaction energies with the MAPs database.
    experiment.incr_OptMAVEn_status()
    OPTMAVEN2.interaction_energies(experiment)

elif status == 1:
    # Check if all of the interaction energy calculations have finished.
    MAPs_parts = list()
    MAPs_dir = os.path.join(InstallFolder, "databases", "MAPs")
    for part_type in os.listdir(MAPs_dir):
        MAPs_parts.extend([os.path.splitext(part)[0] for part in os.listdir(
                os.path.join(MAPs_dir, part_type))])
    energy_directory = os.path.join(directory, "energies")
    for molecule in os.listdir(energy_directory):
        mol_directory = os.path.join(energy_directory, molecule)
        for part in MAPs_parts:
            part_directory = os.path.join(mol_directory, part)
            # If at least one energy calculation has not finished, abort and
            # wait until they are all finished.
            part_file = os.path.join(part_directory, "finished")
            if not os.path.isfile(part_file):
                print "{} is not yet finished.".format(part_file)
                sys.exit(0)

    # If all of the calculations have finished, then reorganize the energies,
    # which are grouped first by part and second by position, to be grouped
    # first by position and then by part.
    OPTMAVEN2.regroup_energies(experiment)

    # Pick the antibody parts combination for each position.
    experiment.incr_OptMAVEn_status()
    OPTMAVEN2.select_all_parts(experiment)

elif status == 2:
    # Combine the files of selected parts into one file for each molecule.
    energy_directory = os.path.join(directory, "energies")
    for molecule in os.listdir(energy_directory):
        mol_directory = os.path.join(energy_directory, molecule)
        lines = list()
        header = False
        for parts_file in [f for f in os.listdir(mol_directory) if f.startswith(
                "parts_")]:
            parts_path = os.path.join(mol_directory, parts_file)
            lines.extend(open(parts_path).readlines()[int(header):])
            header = True
        all_file = os.path.join(mol_directory, "parts_all.csv")
        open(all_file, "w").write("".join(lines))

