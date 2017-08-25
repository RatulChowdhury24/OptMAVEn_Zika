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
import itertools #mfa#
import os
import sys
import numpy as np
import tempfile

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
performance_path = os.path.join(directory, "performance.txt")
performance_file = PERFORMANCE.File(performance_path)

# Determine the status of the experiment.

if experiment.get_OptMAVEn_status() == "initialize":
    # Record the start time.
    performance_file.start(experiment)
    performance_file.disk_usage("initial", directory)
    # Move the antigen to the initial position and generate a PSF.
    OPTMAVEN.initialize_antigen(experiment)
    performance_file.disk_usage("after_antigen_initialization", directory)
    experiment.set_OptMAVEn_status("cull clashes starting")

if experiment.get_OptMAVEn_status() == "cull clashes starting":
    experiment.set_OptMAVEn_status("cull clashes running")
    # Determine which positions do not cause clashes.
    OPTMAVEN.cull_clashes(experiment)
    performance_file.disk_usage("after_culling_clashes", directory)
    experiment.set_OptMAVEn_status("MAPs energies starting")

if experiment.get_OptMAVEn_status() == "MAPs energies starting":
    # Calculate the interaction energies with the MAPs database.
    experiment.set_OptMAVEn_status("MAPs energies submission")
    OPTMAVEN.interaction_energies(experiment)
    performance_file.disk_usage("beginning_MAPs_energy_calculations", directory)
    experiment.set_OptMAVEn_status("MAPs energies callback")

if experiment.get_OptMAVEn_status() == "MAPs energies callback":
    # Check if all of the interaction energy calculations have finished.
    MAPs_parts = list()
    MAPs_dir = os.path.join(InstallFolder, "databases", "MAPs")
    for domain, region in itertools.product("HLK", ("V", "J", "CDR3")):
        part_type = domain + region
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
    experiment.set_OptMAVEn_status("regroup energies starting")

if experiment.get_OptMAVEn_status() == "regroup energies starting":
    experiment.set_OptMAVEn_status("regroup energies running")
    performance_file.disk_usage("after_MAPs_energy_calculations", directory)
    # If all of the calculations have finished, then reorganize the energies,
    # which are grouped first by part and second by position, to be grouped
    # first by position and then by part.
    timer = PERFORMANCE.Timer()
    timer.start()
    OPTMAVEN.regroup_energies(experiment)
    performance_file.sub_time("regrouping_interaction_energies", timer)
    performance_file.disk_usage("after_regrouping_energies", directory)
    # Pick the antibody parts combination for each position.
    experiment.set_OptMAVEn_status("select parts starting")

if experiment.get_OptMAVEn_status() == "select parts starting":
    experiment.set_OptMAVEn_status("select parts submission")
    OPTMAVEN.select_all_parts(experiment)
    experiment.set_OptMAVEn_status("select parts callback")

if experiment.get_OptMAVEn_status() == "select parts callback":
    # Determine the antigen positions to use.
    positions = [tuple(map(float, pos.split()[0: 4])) for pos in open(
            os.path.join(directory, "input_files", "positions.dat"))]
    energy_directory = os.path.join(directory, "energies")
    parts_files = ["parts_zr{}x{}y{}z{}.csv".format(*pos) for pos in positions]
    for molecule in os.listdir(energy_directory):
        mol_directory = os.path.join(energy_directory, molecule)
        parts = set(os.listdir(mol_directory))
        unfinished = [parts_file for parts_file in parts_files if parts_file
                not in parts]
        if len(unfinished) == 0:
            experiment.set_OptMAVEn_status("collate designs starting")
        else:
            print "MAPs energies for {} part(s) is/are not yet finished."\
                    .format(len(unfinished))
            exit()


if experiment.get_OptMAVEn_status() == "collate designs starting":
    experiment.set_OptMAVEn_status("collate designs running")
    # Determine disk usage.
    performance_file.disk_usage("after_MILP", directory)
    positions = [tuple(map(float, pos.split()[0: 4])) for pos in open(
            os.path.join(directory, "input_files", "positions.dat"))]
    parts_files = ["parts_zr{}x{}y{}z{}.csv".format(*pos) for pos in positions]
    # Combine the files of selected parts into one file for each molecule.
    energy_directory = os.path.join(directory, "energies")
    for molecule in os.listdir(energy_directory):
        mol_directory = os.path.join(energy_directory, molecule)
        lines = list()
        unfinished = list()
        header = False
        for parts_file in parts_files:
            parts_path = os.path.join(mol_directory, parts_file)
            lines.extend([line.strip() for line in open(
                    parts_path).readlines()[int(header):]])
            header = True
        # Write a combined parts file if all part selections are complete.
        all_file = os.path.join(mol_directory, "parts_all.csv")
        open(all_file, "w").write("\n".join(lines))
        # Remove the individual parts files to save disk space.
        for pos in positions:
            files = ["parts_zr{}x{}y{}z{}.csv".format(*pos),
                     "zr{}x{}y{}z{}.dat".format(*pos)]
            for f in files:
                try:
                    os.remove(os.path.join(mol_directory, f))
                except OSError:
                    pass
    # Determine disk usage.
    performance_file.disk_usage("after_removing_MILP_files", directory)
    experiment.set_OptMAVEn_status("clustering starting")

 
if experiment.get_OptMAVEn_status() == "clustering starting":
    experiment.set_OptMAVEn_status("clustering running")
    timer = PERFORMANCE.Timer()
    timer.start()
    OPTMAVEN.cluster_designs(experiment)
    # Record the finish time and total time.
    performance_file.sub_time("clustering", timer)
    performance_file.disk_usage("final", directory)
    performance_file.finish(experiment)
    experiment.set_OptMAVEn_status("finished")
    # Email the user when the job is finished.
    # FIXME: have the user enter the address somewhere else
    email_address = "quizmatthew118@gmail.com"
    handle, email_file = tempfile.mkstemp()
    with open(email_file, "w") as f:
        f.write("""Dear """ + experiment["User"] + """,

Your OptMAVEn experiment, """ + experiment["Name"] + """, has finished.

Sincerely,
The Maranas Group
http://www.maranasgroup.com""")
    os.system('cat {} | mail -s "OptMAVEn: {}" {}'.format(email_file, experiment["Name"], email_address))
    os.remove(email_file)


if experiment.get_OptMAVEn_status() == "finished":
    print "{} is finished.".format(experiment["Name"])
