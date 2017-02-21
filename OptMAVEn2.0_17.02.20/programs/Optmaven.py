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
import math
import datetime
import os
import sys
import numpy as np

# Import the IPRO modules.
sys.path.append(os.path.join(InstallFolder, "modules"))
import EXPERIMENT
import IPRO_FUNCTIONS
import OPTMAVEN
import REFINEMENT
import PERFORMANCE

DTIME_FORMAT = "on %Y %B %d at %H:%M:%S"

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

if experiment.get_OptMAVEn_status() == "initialize":
    # Record the start time.
    open(os.path.join(directory, "performance.txt"), "a").write("Started " + 
            datetime.datetime.now().strftime(DTIME_FORMAT))

    # Move the antigen to the initial position and generate a PSF.
    OPTMAVEN.initialize_antigen(experiment)
    experiment.incr_OptMAVEn_status()

if experiment.get_OptMAVEn_status() == "cull clashes":
    # Determine which positions do not cause clashes.
    OPTMAVEN.cull_clashes(experiment)
    experiment.incr_OptMAVEn_status()

if experiment.get_OptMAVEn_status() == "MAPs energies":
    # Calculate the interaction energies with the MAPs database.
    experiment.incr_OptMAVEn_status()
    OPTMAVEN.interaction_energies(experiment)

if experiment.get_OptMAVEn_status() == "MAPs energies callback":
    # Determine disk usage.
    du_out = "du.out"
    while os.path.isfile(du_out):
        du_out = "1" + du_out
    os.system("du -sh {} > {}".format(directory, du_out))
    du = open(du_out).read()
    try:
        os.remove(du_out)
    except OSError:
        pass
    open(os.path.join(directory, "performance.txt"), "a").write("\nDisk usage "
            "during energy calculations is {}".format(du))
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
    timer = PERFORMANCE.Timer(os.path.join(directory, "performance.txt"))
    timer.start("Regrouping MAPs energies")
    OPTMAVEN.regroup_energies(experiment)
    timer.stop()
    timer.dump()

    # Pick the antibody parts combination for each position.
    experiment.incr_OptMAVEn_status()
    OPTMAVEN.select_all_parts(experiment)

if experiment.get_OptMAVEn_status() == "select parts callback":
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
    # Determine disk usage.
    du_out = "du.out"
    while os.path.isfile(du_out):
        du_out = "1" + du_out
    os.system("du -sh {} > {}".format(directory, du_out))
    du = open(du_out).read()
    try:
        os.remove(du_out)
    except OSError:
        pass
    open(os.path.join(directory, "performance.txt"), "a").write("\nDisk usage "
            "after all energy calculations and MILP is {}".format(du))
    # Record the finish time and total time.
    now = datetime.datetime.now()
    started = None
    with open(os.path.join(directory, "performance.txt")) as f:
        for line in f:
            if line.startswith("Started on"):
                try:
                    started = datetime.datetime.strptime(line[8:], DTIME_FORMAT)
                except ValueError:
                    pass
    if started is not None:
        duration = now - started
        open(os.path.join(directory, "performance.txt"), "a").write("\nFinished {}"
                "\nTotal Time: {}".format(now.strftime(DTIME_FORMAT), duration))
    

