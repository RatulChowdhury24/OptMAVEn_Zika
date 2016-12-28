#!/usr/bin/env python

# The name of this file
__name__ = "General IPRO Suite Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions fro running IPRO. Each of the main functions
requires an Experiment class object as its input."""

# Include PYTHON modules
import os
import sys, pdb
import math
import copy
import time
import random
# Include the IPRO Suite STANDARDS module
from STANDARDS import *
# And include other IPRO Suite modules
import MOLECULES
import CHARMM
import EXPERIMENT
import DOCKING_FUNCTIONS
import SHARING
import ROTAMERS
import REFINEMENT

from pore_area import *

def do_superimpose(experiment):
    """Calculate whether or not Molecule superpositions should occur"""
    # Store that information here
    superimpose = {}
    # Initialize that dictionary to False
    for group in experiment:
        superimpose[group.number] = {}
        for molecule in group:
            superimpose[group.number][molecule.name] = False
    # Go through the Restraints to find out
    if "Restraints" in experiment and "Position" in experiment["Restraints"]:
        for restraint in experiment["Restraints"]["Position"]:
            # If the last entry is an integer, the restraint connects between
            # Design Groups, not to the initial positions, and can be skipped
            if isinstance(restraint[-1], int):
                continue
            # Make a list of the Design Groups the restraint applies to
            gns = []
            if restraint[0] == 'all':
                for group in experiment:
                    gns.append(group.number)
            else:
                gns.append(restraint[0])
            # Loop through those Design Groups
            for gn in gns:
                # Determine which Molecules it applies to
                mns = []
                if restraint[1] == 'all':
                    for molecule in experiment[gn]:
                        mns.append(molecule.name)
                else:
                    mns.append(restraint[1])
                # Superimpose those Molecules if they aren't Design Molecules
                for mn in mns:
                    if not experiment[gn][mn].design:
                        superimpose[gn][mn] = True
    # Store this information in the Experiment
    experiment["Superimpose"] = superimpose

def Closest_Rotamers(experiment, gn = None):
    """Assign the Closest Rotamers to amino acids in Design Molecules"""
    # Make sure the experiment is an experiment
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Closest Rotamers function requires an Experiment class "
        text += "object as its input."
        raise IPRO_Error(text)
    # Start timing this
    startTime = time.time()
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # Go through the amino acids of the Design Molecules
        for molecule in group:
            for residue in molecule:
                if molecule.design and residue.kind in \
                aminoAcids[residue.fileFormat]:
                    residue.permission = "ISOMER"
                else:
                    residue.permission = "FIXED"
        # use the ROTAMERS function
        ROTAMERS.Closest_Rotamer(group, experiment)
        # Update the summary
        experiment["Summary"] += SHARING.summary_update(groupTime, \
        "Assigning the Closest Rotamers to Design Group " + str(group.number))
    if gn == None:
        experiment["Summary"] += SHARING.summary_update(groupTime, \
                                 "The Closest Rotamers function")

def standard_position_selection(experiment):
    """The standard method to choose a Design Position for perturbation"""
    # Make a list of all of the Design Positions
    spots = []
    for mn in experiment["Design Positions"]:
        for rn in experiment["Design Positions"][mn]:
            spots.append([mn, rn])
    # Choose the perturbation position
    return spots[random.randint(0, len(spots) - 1)]

def optmaven_position_selection(experiment):
    """ The optmaven method to choose a Design Position for perturbation """
    # Make a list of all of the Design Positions
    spots = []
    for entry in experiment["Molecules"]:
        if entry[2].design:
            moleculeName = entry[1]
            molecule = entry[2]
            for residue in molecule:
        #Extract each residue name and molecule they belong to in an antibody and save them to a list
        #which can be used to random selecting the perturbed residues for OPTMAVEN refinement
        #If a residue belongs to the CDRs, this residue is saved three times
        #which will make sure CDRs residues can be bias randomly selected with three times higher probablity
        #than the non-CDRs residues.
                try:
                    r = int(residue.name)
                except ValueError:
                    r = int(residue.name[:-1])
                if r in range(27, 39) or r in range(56, 66) or r in range(105, 118):
                    N = 3
                else:
                    N = 1
                for i in range(N):
                    spots.append([moleculeName, residue.name])
        # Choose the perturbation position
    return spots[random.randint(0, len(spots) - 1)]

def select_perturbation_position(experiment):
    """Select a Design Position for perturbation"""
    # Do this differently depending no the Experiment's type
    if experiment["Type"] in ['IPRO', 'Mutator']:
        return standard_position_selection(experiment)
    elif experiment["Type"] in ['OptMAVEn']:
        return optmaven_position_selection(experiment)
    else:
        text = "The select perturbation postion function does not support the "
        text += experiment["Type"] + " IPRO Suite Experiment type."
        raise IPRO_Error(text)

def standard_perturbation_selection(group, spot, number = 4):
    """The standard method for choosing a perturbation in IPRO"""
    # Choose how many extra spots there will be
    more = random.randint(0, number)
    # And how many come before
    before = random.randint(0, more)
    after = more - before
    # Get the selected Molecule
    molecule = group[spot[0]]
    # Get the index of the Residue in the Molecule
    for i in range(len(molecule)):
        if molecule[i].name == spot[1]:
            index = i
            break
    # Get the indexes of the perturbed Residues
    i1 = index - before
    if i1 < 0:
        i1 = 0
    i2 = index + after
    if i2 > len(molecule) - 1:
        i2 = len(molecule) - 1
    # Store those perturbed Residues
    perturbed = []
    for i in range(i1, i2+1):
        perturbed.append([molecule.name, molecule[i].name])
    return perturbed

def select_perturbation_region(experiment, group, spot):
    """Select the perturbation region"""
    # Do this based on the Experiment's type
    if experiment["Type"] in ['IPRO', "Mutator"]:
        return standard_perturbation_selection(group, spot)
    elif experiment["Type"] in ['OptMAVEn']:
        # Mave sure during optmaven maximum perturbed residues number are 3 instead of 5 for regular IPRO perturbation
        number = 1
        return standard_perturbation_selection(group, spot, number)

    else:
        text = "The select perturbation region function does not support the "
        text += experiment["Type"] + " IPRO Suite Experiment type"
        raise IPRO_Error(text)

def generate_random_angle():
    """Generate a random perturbation angle"""
    # Store the value here
    angle = 6
    while not -5 <= angle <= 5:
        angle = random.gauss(0, 1.5)
    return angle

def generate_angles(group, perturbed):
    """Generate the random dihedral angle perturbations"""
    # Store them here
    angles = {}
    # Go through the perturbed Residues
    for spot in perturbed:
        # Get the Molecule
        molecule = group[spot[0]]
        # Initialize angles
        if molecule.name not in angles:
            angles[molecule.name] = {}
        # Get the Residue
        residue = molecule[spot[1]]
        angles[molecule.name][residue.name] = {}
        # If appropriate, generate phi and psi angles
        if residue.name != molecule[0].name:
            angles[molecule.name][residue.name]["phi"] = generate_random_angle()
        if residue.name != molecule[-1].name:
            angles[molecule.name][residue.name]["psi"] = generate_random_angle()
    return angles

def standard_sequence_setting(experiment, group, perturbed):
    """Set freedoms and permissions based on sequence"""
    # Set everything to fixed in place
    for molecule in group:
        for residue in molecule:
            residue.freedom = "FIXED"
            residue.permission = "FIXED"
            # And store what the current kind of amino acid is for every Residue
            residue.currentKind = residue.kind
    # Go through the perturbed Residues
    for spot in perturbed:
        # Get the Molecule and the index of the Residue
        molecule = group[spot[0]]
        for i in range(len(molecule)):
            if molecule[i].name == spot[1]:
                index = i
                break
        # Set that Residue's values appropriately
        molecule[index].freedom = "RESTRAINED"
        # Only Design Positions can mutate, and only when a Standard IPRO
        # iteration is happening
        if molecule[index].design and experiment["Activity"] == "Standard":
            molecule[index].permission = "ROTAMER"
        else:
            molecule[index].permission = "ISOMER"
        # Go through Residues on either side of that perturbed Residue
        for i in range(index - 5, index + 6):
            # Skip indexes that are outside of the Molecule's range
            if i < 0 or i >= len(molecule):
                continue
            # If the Residue is fixed in place, let it move
            if molecule[i].freedom == "FIXED":
                molecule[i].freedom = "FREE"
                # ROTAMERS won't get a rotamer for a non-amino acid Residue, so
                # this is fine
                molecule[i].permission = "ISOMER"

def sequence_setting(experiment, group, perturbed):
    """Set the freedoms of the Residues by sequence"""
    # Use the standard method
    if experiment["Type"] in ['IPRO', 'Mutator', 'OptMAVEn']:
        standard_sequence_setting(experiment, group, perturbed)
    else:
        text = "The sequence setting function does not support the "
        text += experiment["Type"] + " IPRO Suite Experiment type"
        raise IPRO_Error(text)

def distance_setting(experiment, group):
    """Modify Residue freedoms and permissions based on distances"""
    # Only do this if the packing method is 'Distance'
    if experiment["Packing Method"] == "Distance":
        # Make a list of the Residue freedoms that are acceptable
        freedoms = ["RESTRAINED"]
        if experiment["Packing Selection"] == "FREE":
            freedoms.append("FREE")
        # Get all of the heavy Atoms in the docking selection
        atoms = []
        for molecule in group:
            for residue in molecule:
                if residue.freedom in freedoms:
                    for atom in residue:
                        if not ROTAMERS.is_hydrogen(atom):
                            atoms.append(atom)
        # Modify permissions in Design Molecules
        for molecule in group:
            if not molecule.design:
                continue
            for residue in molecule:
                # If the Residue already has settings, skip it
                if residue.freedom != "FIXED":
                    continue
                # Determine if it is sufficiently close
                close = False
                for atom in residue:
                    # Only consider heavy atoms
                    if ROTAMERS.is_hydrogen(atom):
                        continue
                    for atom2 in atoms:
                        dis = MOLECULES.calculate_distance(atom, atom2)
                        if dis <= experiment["Packing Cutoff"]:
                            close = True
                            break
                    if close:
                        break
                # If the Residue is close to the Packing selection, repack it
                if close:
                    residue.freedom = "FREE"
                    residue.permission = "ISOMER"

def dimer_match(mol1, mol2, angles):
    """When Dimers are present, make sure they are both modified."""
    # Match the freedoms and permissions of the Molecules that have been
    # identified as dimers. Also modify the angles dictionary
    for residue in mol1:
        # If this Residue is unavailable in the second molecule, skip it
        if residue.name not in mol2:
            continue
        # Match the angles dictionary for this Residue to the one in the second
        if mol1.name in angles and residue.name in angles[mol1.name]:
            if mol2.name not in angles:
                angles[mol2.name] = {}
            if residue.name not in angles[mol2.name]:
                angles[mol2.name][residue.name]=angles[mol1.name][residue.name]
        # Modify the freedom and permission
        if residue.freedom == "FIXED":
            pass
        elif residue.freedom == "FREE" and mol2[residue.name].freedom in \
        ["FIXED", "FREE"]:
            mol2[residue.name].freedom = "FREE"
        else:
            mol2[residue.name].freedom = "RESTRAINED"
        if residue.permission == "FIXED":
            pass
        elif residue.permission == "ISOMER" and mol2[residue.name].permission \
        in ["FIXED", "ISOMER"]:
            mol2[residue.name].permission = "ISOMER"
        else:
            mol2[residue.name].permission = "ROTAMER"
    # This function will be called twice, with the order of the Molecules
    # reversed, so there's no need to duplicate the information

def Dimers(experiment, group, angles):
    """Match freedom, permission, and perturbation information for Dimers"""
    # If there are Dimers
    if "Dimers" in experiment:
        # The dimers can ONLY be Design Molecules and every Design Molecule must
        # be in each Design Group, so everything can be matched up
        for pair in experiment["Dimers"]:
            dimer_match(group[pair[0]], group[pair[1]], angles)
            dimer_match(group[pair[1]], group[pair[0]], angles)

def refinement_check(experiment, gn):
    """Check to see if a refinement is ongoing"""
    # Only do this check under very specific conditions
    if experiment["Activity"] == "Standard" and experiment["Do Refinement"] \
    and gn == None:
        if "refinement" in os.listdir(experiment["Folder"]):
            return True
    return False

def Backbone_Perturbation(experiment, gn = None):
    """Perturb the Backbone of a Molecule"""
    # Make sure the experiment is an Experiment
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Backbone Perturbation function requires an Experiment class"
        text += " object as its input"
        raise IPRO_Error(text)
    # Find out if a refinement is happening
    refinement = refinement_check(experiment, gn)
    if refinement:
        return refinement
    # Assign the closest rotamers to all Residues
    #Closest_Rotamers(experiment, gn)
    # Start timing this
    startTime = time.time()
    # Determine which Group to make the perturbation selections for
    if gn == None:
        N = 1
    else:
        N = gn
    # Choose the perturbation location
    spot = select_perturbation_position(experiment)
    # Choose the perturbation region
    experiment["Summary"] += "Residue " + spot[1] + " in Molecule " + spot[0]
    experiment["Summary"] += " selected as the Perturbation Position\n"
    perturbed = select_perturbation_region(experiment, experiment[N], spot)
    if len(perturbed) > 1:
        items = []
        for spot in perturbed:
            items.append("Residue " + spot[1] + " in Molecule " + spot[0])
        text = screen_formatting("Perturbing" + list_items(items) + "\n")
        experiment["Summary"] += text
    # Generate the random perturbation angles
    angles = generate_angles(experiment[N], perturbed)
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # Assign freedoms and permissions
        if group.number == N:
            # Set based on sequence
            sequence_setting(experiment, group, perturbed)
            # Set based on distance (will only work if appropriate)
            distance_setting(experiment, group)
            # Match Dimer information
            Dimers(experiment, group, angles)
        # Otherwise, match things to the Number group
        else:
            for m in group:
                for r in m:
                    # If this is a Target Molecule, freeze it
                    if not m.design:
                        r.freedom = "FIXED"
                        r.permission = "FIXED"
                    else:
                        r.currentKind = r.kind
                        r.freedom = experiment[N][m.name][r.name].freedom
                        r.permission = experiment[N][m.name][r.name].permission
        # Mutate to Glycines
        ROTAMERS.Glycine(group, experiment)
        # Perturb the backbone
        if experiment["Force Field"] == "CHARMM":
            CHARMM.Perturbation(group, angles, experiment)
        else:
            text = "The Backbone Perturbation function does not support the "
            text += str(experiment["Force Field"]) + " force field."
            raise IPRO_Error(text)
        # Update the summary
        if gn == None:
            experiment["Summary"] += SHARING.summary_update(groupTime, \
                        "Perturbing Design Group " + str(group.number))
    experiment["Summary"] += SHARING.summary_update(startTime, \
                             "The Backbone Perturbation function")
    return refinement

def superimpose(current, original):
    """Superimpose the original copy of the Molecule on its current structure"""
    # Only do this if the two Molecules evaluate to equal
    if current.__eq__(original):
        # Center both Molecules
        move1 = current.center("all")
        move2 = original.center("all")
        # We know both Molecules have the same Atoms, so store them
        atoms1 = []
        atoms2 = []
        for r in current:
            for a in r:
                atoms1.append(a)
                atoms2.append(original[r.name][a.name])
        # Calculate the rotation matrix
        rmatrix = MOLECULES.calculate_rmatrix(atoms1, atoms2)
        MOLECULES.rotate(original, rmatrix)
        # Position everything properly
        MOLECULES.move(current, move1, "+")
        MOLECULES.move(original, move1, "+")

def Relaxation(experiment, gn = None, all = True):
    """Do an energy minimization"""
    # Make sure we have an Experiment
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Relaxation function requires an Experiment class object as "
        text += "its input."
        raise IPRO_Error(text)
    # Check to see if a refinement is happening
    refinement = refinement_check(experiment, gn)
    if refinement:
        return refinement
    # Start timing
    startTime = time.time()
    # Loop through the relevant Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # If appropriate, free all Residues
        if all:
            for molecule in group:
                for residue in molecule:
                    residue.freedom = "FREE"
        # If appropriate, superimpose the reference structures of Molecules
        for molecule in group:
            # Design Molecules never need to have structures superimposed
            if not molecule.design:
                if experiment["Superimpose"][group.number][molecule.name]:
                    superimpose(molecule, experiment[0][molecule.name])
        # Do the relaxation
        if experiment["Force Field"] == "CHARMM":
            CHARMM.Relaxation(group, experiment)
        else:
            text = "The Relaxation function does not support the "
            text += str(experiment["Force Field"]) + " force field."
            raise IPRO_Error(text)
        # update the summary
        if gn == None:
            experiment["Summary"] += SHARING.summary_update(groupTime, \
            "The Relaxation of Design Group " + str(group.number))
    experiment["Summary"] += SHARING.summary_update(startTime, \
                             "The Relaxation function")
    return refinement

def Calculate_Energy(experiment, gn = None):
    """Calculate the energies of an IPRO iteration"""
    # Do the standard introductory checks
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Calculate Energy function requires an Experiment class "
        text += "object as its input."
        raise IPRO_Error(text)
    refinement = refinement_check(experiment, gn)
    if refinement:
        return {}, refinement
    # Start timing
    startTime = time.time()
    # Store the energies here
    energies = {}
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # initialize the dictionary
        energies[group.number] = {}
        # If there are no Target Molecules, it is impossible to calculate an
        # interaction energy
        do = False
        for molecule in group:
            if not molecule.design:
                do = True
                break
        # Calculate the energies
        if experiment["Force Field"] == "CHARMM":
            complex = CHARMM.Energy(group, experiment, "all")
            energies[group.number]["Complex"] = complex
            if do:
                design = CHARMM.Energy(group, experiment, True)
                target = CHARMM.Energy(group, experiment, False)
                energies[group.number]["Interaction"] = complex-design-target
        else:
            text = "The Calculate Energy function does not support the "
            text += str(experiment["Force Field"]) + " force field."
            raise IPRO_Error(text)
        # update the summary
        if gn == None:
            experiment["Summary"] += SHARING.summary_update(groupTime, \
            "Energy Calculation for Design Group " + str(group.number))
        experiment["Summary"] +=SHARING.format_energies(energies[group.number],\
                                group.number, False)
    experiment["Summary"] += SHARING.summary_update(startTime, \
                             "The Calculate Energy function")
    return energies, refinement

def Optimal_Rotamers(experiment, gn = None):
    """Select an optimal combination of Rotamers for a system"""
    # Do the standard starting checks
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Optimal Rotamers function requires an Experiment class "
        text += "object as its input."
        raise IPRO_Error(text)
    refinement = refinement_check(experiment, gn)
    if refinement:
        return refinement
    # Start timing this
    startTime = time.time()
    # Determine the group that will have mutations done for it
    if gn == None:
        N = 1
    else:
        N = gn
    # Do this for each Design Group
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Time the group
        groupTime = time.time()
        # Make sure the amino acid sequences match across Design Groups
        if group.number != N:
            for molecule in group:
                for residue in molecule:
                    if residue.design and residue.permission == "ROTAMER":
                        residue.permittedKinds = \
                        [experiment[N][molecule.name][residue.name].kind]
        # Pick the Rotamers
        probability = True
        objective, solution = ROTAMERS.Optimal_Rotamers(group, experiment, probability)
        # Summarize the results
        message = ''
        if group.number == N:
            for molecule in group:
                for residue in molecule:
                    if residue.permission == "ROTAMER":
                        # Summarize the Residue
                        message += "Residue " + residue.name + " in Molecule "
                        message += molecule.name + ": " + residue.kind + "\n"
        message += "The objective value of the MILP was " + format(objective, \
                   '.3f') + "\n"
        experiment["Summary"] += SHARING.summary_update(groupTime, \
        "Selecting Rotamers for Design Group " + str(group.number), message)
    return refinement

def Docking(experiment, iteration, gn = None):
    """Carry out a local, rigid body docking of Target Molecules"""
    # Do the standard initial checks
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Docking function requires an Experiment class object as its"
        text += " input."
        raise IPRO_Error(text)
    refinement = refinement_check(experiment, gn)
    if refinement:
        return False, refinement
    # Make sure Docking should be run in this iteration
    if iteration % experiment["Docking Frequency"] != 0:
        return False, refinement
    # Time this
    startTime = time.time()
    # Collect the groups of Molecules that will move during docking
    dockingGroups = DOCKING_FUNCTIONS.collect_docking_groups(experiment, gn)
    # If there aren't any, say that
    if len(dockingGroups) == 0:
        experiment["Summary"] += "There were no Molecules that could move "
        experiment["Summary"] += "during docking.\n"
        return False, refinement
    # Store the names of all files generated by this function
    fileNames = []
    # Create all of the information and structures needed to run Docking
    information, movingMolecules = DOCKING_FUNCTIONS.create_moving_structures(\
                                   experiment, dockingGroups, fileNames)
    DOCKING_FUNCTIONS.create_static_structures(experiment, movingMolecules, \
                                               fileNames)
    DOCKING_FUNCTIONS.make_docking_information(experiment, dockingGroups, \
                                               information, fileNames)
    # Run docking
    i = os.system("./docking.out")
    if i != 0:
        text = "There was an error running the executable docking.out file."
        raise IPRO_Error(text)
    # Load the results from docking
    message = DOCKING_FUNCTIONS.load_docking_results(experiment, dockingGroups,\
                                                     fileNames)
    DOCKING_FUNCTIONS.load_docking_structures(experiment, dockingGroups)
    # Delete all of the files made by docking
    for fileName in fileNames:
        try:
            os.remove(fileName)
        except OSError:
            pass
    # Update the summary
    experiment["Summary"] += SHARING.summary_update(startTime, \
                             "The Docking function", message)
    # Return that docking ran and that a refinement is not currently ongoing
    return True, refinement

def max_iterations(experiment):
    """Determine the maximum current number of IPRO iterations"""
    if experiment["Activity"] == "Standard":
        return experiment["IPRO Iterations"]
    else:
        return experiment["Refinement Iterations"]

def Start_Iteration(experiment, gn = None):
    """Start an iteration of IPRO"""
    # Do the standard initial checks
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The Start Iteration function requires an Experiment class "
        text += "object as its input."
        raise IPRO_Error(text)
    refinement = refinement_check(experiment, gn)
    if refinement:
        return 0, refinement
    # Start sharing so this processor has exclusive access to the shared files
    SHARING.Start(experiment)
    # Determine what iteration this will be
    iteration = SHARING.iteration_counter(SHARING.get_current(), True) + 1
    # Determine what the expected maximum number of iterations is
    ITERATIONS = max_iterations(experiment)
    # Create a folder
    while iteration <= ITERATIONS:
        folder = "iteration" + str(iteration)
        do = SHARING.claim_calculations(folder)
        if do:
            break
        else:
            iteration += 1
    # If a folder was made
    if iteration <= ITERATIONS:
        # Start a summary
        experiment["Summary"] = "Started" + SHARING.time_stamp()
        # Get the last time structures were updated
        N = SHARING.identify_last_update("./Current/")
        # If necessary, update the structures
        if N > experiment["Last Update"]:
            SHARING.update_Current(experiment, "./Current/", gn)
            experiment["Last Update"] = N
        # Say what structures are being used
        if experiment["Last Update"] == 0:
            experiment["Summary"] += "Using the initial structures\n"
        else:
            experiment["Summary"] += "Using the structures from iteration "
            experiment["Summary"] += str(experiment["Last Update"]) + "\n"
        # Move into the appropriate folder
        os.chdir(folder)
        # Copy in the C++ and force field files
        SHARING.copy_standard_files(experiment)
        # Update the Experiment's structures
        for group in experiment:
            if gn not in [None, group.number]:
                continue
            for molecule in group:
                text =format(experiment["Current"][group.number][molecule.name])
                molecule.load(text)
    # End sharing
    SHARING.End(experiment)
    return iteration, refinement

def standard_decision(experiment, energies, iteration, gn):
    """Make a decision about retaining the results of an iteration"""
    # Keep track of whether the results are the best or kept by SA
    best = True
    keep = True
    # Determine what the maximum number of iteratons should be
    ITERATIONS = max_iterations(experiment)
    # Determine if simulated annealing should be used
    if iteration < ITERATIONS / 10:
        SA = False
    else:
        SA = True
        # Calculate the window for simulated annealing
        window = -(experiment["IPRO Annealing Temperature"] * GasConstant) * \
                   math.log(random.random())
    # Loop through all relevant design groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # If this is a standard comparison
        if experiment["Activity"] == "Standard":
            # Get the objective
            objective = group.objective
            # If the goal is to make binding better, don't do so at the expense
            # of significant penalties in the complex energy
            if objective in ['improve', 'maintain']:
                complex = experiment["Energies"][group.number]["Complex"]
                # too much worse is 40 kcal / mol or 5%, whichever is greater
                value = math.fabs(0.05 * complex)
                if value < 40:
                    value = 40.0
                # If the energy has worsened too much
                if energies[group.number]["Complex"] > complex + value:
                    best = False
                    keep = False
                    break
            # Compare the interaction energies
            current = energies[group.number]["Interaction"]
            reference = experiment["Energies"][group.number]["Interaction"]
        # Otherwise, just try to improve the complex energy
        else:
            objective = "improve"
            current = energies[group.number]["Complex"]
            reference = experiment["Energies"][group.number]["Complex"]
        # Make the decision differently depending on the objective
        if objective in ['improve', 'maintain']:
            if current <= reference:
                pass
            elif SA and current <= reference + window:
                best = False
            else:
                best = False
                keep = False
                break
        else:
            if current >= reference:
                pass
            elif SA and current >= reference - window:
                best = False
            else:
                best = False
                keep = False
                break
    return best, keep

def standard_decision_ompf(experiment, energies, iteration, gn):

    """Make a decision about retaining the results of an iteration"""
    print "ompf is in experiment Name!\n"
    # Keep track of whether the results are the best or kept by SA
    best = True
    keep = True
    # Determine what the maximum number of iteratons should be
    ITERATIONS = max_iterations(experiment)
    # Determine if simulated annealing should be used
    if iteration < ITERATIONS / 10:
        SA = False
    else:
        SA = True
        # Calculate the window for simulated annealing
        window = -(experiment["IPRO Annealing Temperature"] * GasConstant) * \
                   math.log(random.random())
    # Loop through all relevant design groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # If this is a standard comparison
        if experiment["Activity"] == "Standard":
            # Get the objective
            objective = group.objective
            # If the goal is to make binding better, don't do so at the expense
            # of significant penalties in the complex energy
            if objective in ['improve', 'maintain']:
                complex = experiment["Energies"][group.number]["Complex"]
                # too much worse is 40 kcal / mol or 5%, whichever is greater
                value = math.fabs(0.05 * complex)
                if value < 40:
                    value = 40.0
                # If the energy has worsened too much
                if energies[group.number]["Complex"] > complex + value:
                    best = False
                    keep = False
                    break
            # Compare the interaction energies
            current = energies[group.number]["Interaction"]
            reference = experiment["Energies"][group.number]["Interaction"]
        # Otherwise, just try to improve the complex energy
        else:
            objective = "improve"
            current = energies[group.number]["Complex"]
            reference = experiment["Energies"][group.number]["Complex"]
        # Make the decision differently depending on the objective
        if objective in ['improve', 'maintain']:
            if current <= reference:
                pass
            elif SA and current <= reference + window:
                best = False
            else:
                best = False
                keep = False
                break
        elif objective in ['eliminate']:

            pdbfile = os.path.join(experiment['folder'],'Current','Group1_MoleculeZ.pdb')
            
            experiment["Summary"] += "************\n"

            sumAB = []
            for i in range(-5,5):
                A, B = pore_area(pdbfile, VDW_radius=0.8, AngleSearchSize=10, RadiusSearchSize=20, Origin=[0,i,0])
                sumAB.append(A+B)

            min_sumAB = min(sumAB)
            experiment["Summary"] += "The min_sumAB for this iteration is %.3f\n"%min_sumAB


            if current >= reference:
                if min_sumAB < 7:
                    best = True
                    #experiment["Summary"] += "The min_sumAB for this iteration is %.3f\n"%min_sumAB
                    keep = False
                elif min_sumAB < 12:
                    best = False
                    keep = True
                    #experiment["Summary"] += "The min_sumAB for this iteration is %.3f\n"%min_sumAB
                else:
                    best = False
                    keep = False
            elif SA and current >= reference - window and min_sumAB < 12: 
                best = False
                #experiment["Summary"] += "The min_sumAB for this iteration is %.3f\n"%min_sumAB
            else:
                best = False
                keep = False
                break                            
    return best, keep


def optzyme_decision(experiment, energies, iteration, gn):
    """Make a decision about retaining the results from an OptZyme iteration"""
    # Keep track of whether the results are the best or if they are kept by
    # SA
    best = True
    keep = True
    # Determine the maximum number of iterations used in the OptZyme experiment
    ITERATIONS = max_iterations(experiment)
    # If OptZyme is over 10% finished, use simulated annealing
    if iteration < ITERATIONS / 10:
        SA = False
    else:
        SA = True
        # Calculate the window for simulated annealing
        window = -(experiment["IPRO Annealing Temperature"] * GasConstant) * \
                   math.log(random.random()) 
    # Loop through the OptZymeGroups
    for group in experiment["Optzyme Information"]:
        # Get the OptzymeGroup's objective
        objective = group.objective
        # Ensure that the BindingAssemblies are consistent with the OptzymeGroup
        property = group.property
        error = ''
        # A group optimizing km or kcat/km should have 1 BindingAssembly
        if property in ["km", "kcat/km"]:
            if len(group) != 1:
                error = "OptzymeGroup " + str(group.number) + " is "
                error += "optimizing " + property + " and should have 1 "
                error += "BindingAssembly but instead has " + str(len(group)) 
        # A group optimizing kcat needs 2 BindingAssemblies arranged properly
        elif property in ["kcat"]:
            if len(group) == 2 and group[0].objective == "improve" and \
            group[1].objective == "eliminate":
                pass
            else:
                error = "OptzymeGroup " + str(group.number) + " is "
                error += "optimizing " + property + " and should have 2 "
                error += "BindingAssemblies. The first Binding Assembly "
                error += "(TSA) should have 'improve' as the objective "
                error += "and the second Binding Assembly (substrate) should "
                error += "have 'eliminate' as the objective"
        # No other property should ever exist, but double check
        else:
            error = property + " is not an acceptable property for optimization"
            error += " in OptZyme"
        # If an error was found, terminate
        if error != '':
            raise IPRO_Error(text) 
        # Loop through all the relevant BindingAssemblies
        for assembly in group:             
            if gn not in [None, assembly.number]:
                continue
            # If this is a standard comparison
            if experiment["Activity"] == "Standard":
                if objective in ["improve", "maintain"]:
                    complex = experiment["Energies"][assembly.number]["Complex"]
                    # Too much worse is 40 kcal/mol or 5% of the complex energy,
                    # whichever is greater
                    value = math.fabs(0.05 * complex)
                    if value < 40:
                        value = 40.0
                    # If the complex energy has worsened too much, discard
                    # the results
                    if energies[assembly.number]["Complex"] > complex + value:
                        best = False
                        keep = False
                        break
                # Compare using the proper interaction energy term
                if property in ["km", "kcat/km"]:
                    current = energies[assembly.number]["Interaction"]
                    reference = \
                    experiment["Energies"][assembly.number]["Interaction"]
                # If the property is kcat, use the weighted difference of the
                # two interaction energies
                else:
                    current = energies[group[0]]["Interaction"] - \
                              group.weight * energies[group[1]]["Interaction"]
                    reference = experiment["Energies"][group[0]]["Interaction"]
                    reference = reference - group.weight * \
                                experiment["Energies"][group[1]]["Interaction"]
            # Otherwise, just try to improve the complex energy
            else:
                objective = "improve"
                current = energies[assembly.number]["Interaction"]
                reference = experiment["Energies"][assembly.number]["Complex"]
            # Make the decision differently depending on the objective    
            # Improve KM (small KM) -> lower IES
            # Improve kcat/KM (large kcat/KM) -> lower IETSA
            # Improve kcat (large kcat) -> lower IETSA - (RTTSA/RTS)*IES
            if objective in ['improve', 'maintain']:
                if current <= reference:
                    pass
                elif SA and current <= reference + window:
                    best = False
                else:
                    best = False
                    keep = False
                    break
            else:
                if current >= reference:
                    pass
                elif SA and current >= reference - window:
                    best = False
                else:
                    best = False
                    keep = False
                    break
    return best, keep                    



def iteration_decision(experiment, energies, iteration, gn = None):
    """Make a decision regarding whether or not to keep an iteration"""
    # If the normal decision function is needed

    if experiment["Type"] in ['IPRO', 'Mutator', 'OptMAVEn']:
        if 'ompf' in experiment["Name"]:
            print "I work!"
            #pdb.set_trace()
            best, keep = standard_decision_ompf(experiment, energies, iteration, gn)
        else:
            print "Standard..."
            #pdb.set_trace()
            best, keep = standard_decision(experiment, energies, iteration, gn)
    elif experiment["Type"] in ['OptZyme']:
        best, keep = optzyme_decision(experiment, energies, iteration, gn)

    else:
        text = "The iteration decision function does not support the "
        text += experiment["Type"] + " IPRO Suite Experiment type."
        raise IPRO_Error(text)
    return best, keep

def store_structures(experiment, gn = None):
    """Store the structures of the Experiment's Design Groups"""
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Go through the Moleculees
        for molecule in group:
            text = format(molecule)
            experiment["Current"][group.number][molecule.name].load(text)

def store_energies(experiment, energies, gn = None):
    """Store energies in the Experiment"""
    # Loop through the Design Groups
    for group in experiment:
        if gn not in [None, group.number]:
            continue
        # Reset the Energies dictionary for this group
        experiment["Energies"][group.number] = {}
        for energy in energies[group.number]:
            experiment["Energies"][group.number][energy] = \
                          energies[group.number][energy]

def End_Iteration(experiment, energies, iteration, gn):
    """End an iteration of IPRO."""
    # The standard initial checks
    if not isinstance(experiment, EXPERIMENT.Experiment):
        text = "The End Iteration function requires an Experiment class object "
        text += "as its input."
        raise IPRO_Error(text)
    refinement = refinement_check(experiment, gn)
    if refinement:
        return refinement
    # Move out of the iteration's folder
    os.chdir("../")
    # Start sharing
    SHARING.Start(experiment)
    # Load the appropriate reference energies
    if experiment["Activity"] == "Standard":
        SHARING.load_reference_Energies(experiment)
    else:
        SHARING.update_Energies(experiment, "./Best/", gn)
    # Determine the last completed iteration
    this = SHARING.iteration_counter(SHARING.get_current(), False) + 1
    # Determine whether or not keep the iteration results
    best, keep = iteration_decision(experiment, energies, this, gn)
    # If results should be kept, do so
    if best or keep:
        store_structures(experiment, gn)
        # If these are the best results so far, store the energies
        if best:
            store_energies(experiment, energies, gn)
    # If appropriate, share the results with other processors
    SHARING.Results(experiment, this, best, keep, gn)
    # Say what happened in the summary
    experiment["Summary"] = "\nIteration " + str(this) + "\n" + \
                            experiment["Summary"] + "The results were "
    if best:
        experiment["Summary"] += "the BEST so far\n"
    elif keep:
        experiment["Summary"] += "KEPT by simulated annealing\n"
    else:
        experiment["Summary"] += "DISCARDED\n"
    # Put a time stamp on the summary
    experiment["Summary"] += "Ended" + SHARING.time_stamp()
    # Get the name of the Summary file
    name = SHARING.summary_name(SHARING.get_current())
    f = open(name, "a")
    f.write(experiment["Summary"])
    f.close()
    # Delete the iteration
    os.system("rm -rf iteration" + str(iteration))
    # End sharing
    SHARING.End(experiment)
    return refinement

def IPRO_ITERATION(experiment, gn = None):
    """Do an iteration of IPRO."""
    # Start the iteration
    iteration, refinement = Start_Iteration(experiment, gn)
    if not refinement and iteration > max_iterations(experiment):
        return iteration
    # Do the steps of an iteration
    if not refinement:
        refinement = Backbone_Perturbation(experiment, gn)
    if not refinement:
        refinement = Optimal_Rotamers(experiment, gn)
    if not refinement:
        docked, refinement = Docking(experiment, iteration, gn)
    if not refinement:
        refinement = Relaxation(experiment, gn, True)
    if not refinement:
        energies, refinement = Calculate_Energy(experiment, gn)
    if not refinement:
        refinement = End_Iteration(experiment, energies, iteration, gn)
    # If a refinement has been started and an iteration folder was created,
    # delete the iteration folder
    if refinement and iteration > 0:
        os.chdir("../")
        os.system("rm -rf iteration" + str(iteration))
    return iteration

def initialize_group(experiment, gn):
    """Initialize a Design Group"""
    # Create the folder's name
    folder = "Group" + str(gn)
    # Try to claim it for calculations
    do = SHARING.claim_calculations(folder)
    # If this processor is doing the calculations
    if do:
        # Make a summary
        experiment["Summary"] = "Started" + SHARING.time_stamp()
        # Move into the folder and copy in files
        os.chdir(folder)
        SHARING.copy_standard_files(experiment)
        # Relax everything
        refinement = Relaxation(experiment, gn, True)
        # Assign closest rotamers
        Closest_Rotamers(experiment, gn)
        # Do another relaxation
        refinement = Relaxation(experiment, gn, True)
        # Calculate the initial energies
        energies, refinement = Calculate_Energy(experiment, gn)
        # Store the structures and energies
        store_structures(experiment, gn)
        store_energies(experiment, energies, gn)
        # Create a summary file
        experiment["Summary"] += "Ended" + SHARING.time_stamp()
        name = SHARING.summary_name(SHARING.get_current())
        f = open(name, "w")
        f.write(experiment["Summary"])
        f.close()
        # Move up a folder
        os.chdir("../")
        # Start sharing
        SHARING.Start(experiment)
        # output the structures and energies to the Current folder
        SHARING.output_Current(experiment, "./Current/", gn)
        SHARING.output_Energies(experiment, "./Current/", gn)
        # End sharing
        SHARING.End(experiment)
    return do

def initialize_optzyme_groups(experiment):
    """Create the OptZyme Groups for the experiment."""
    # Store the Optzyme Groups in a list
    groups = []
    # Go through the Optzyme Groups
    for i, group in enumerate(experiment["Optzyme Groups"]):
        # Find the corresponding Binding Assemblies
        assemblies = []
        # Go through each Assembly within the Optzyme Group
        for ba in group[3]:
            # State whether the matching Assembly has been found
            found = False
            # Go through the stored Binding Assemblies 
            for assembly in experiment:
                # Get the list of Target Molecules in the assembly
                target = []
                for mol in assembly:
                    if not mol.design:
                        target.append(mol.name)
                # See if the list matches the one provided from the Optzyme
                # Group
                if target.sort() == ba[1:].sort():
                    found = True
                    break
            # Raise an error if the Binding Assembly was not found        
            if not found:
                text = "A matching Binding Assembly for the inputs "
                text += "provided was not found"
                raise IPRO_Error(text)
            # Otherwise, store the information in the 'assemblies' list
            else:
                assemblies.append(assembly)
        # Generate the Optzyme Group
        optzymeGroup = MOLECULES.OptzymeGroup(i+1, assemblies, \
                        experiment["Force Field"], experiment["File Format"])
        # Store the attributes of the Optzyme Group
        optzymeGroup.property = group[0]
        optzymeGroup.objective = group[1]
        if group[0] != "kcat":
            optzymeGroup.weight = group[2]
        # Append this Optzyme Group to the list of Groups
        groups.append(optzymeGroup)
    experiment["Optzyme Information"] = groups 

def initialize_molecule(experiment, mn):
    """Initialize a Molecule"""
    # Create the folder's name
    folder = "Molecule" + mn
    do = SHARING.claim_calculations(folder)
    # If this processor is doing the calculations
    if do:
        # Make a summary
        experiment["Summary"] = "Started" + SHARING.time_stamp()
        # Move into the folder
        os.chdir(folder)
        SHARING.copy_standard_files(experiment)
        # Get a unique copy of the Molecule
        molecule = experiment[0][mn].duplicate()
        for residue in molecule:
            residue.freedom = "FREE"
        # Store the calculated energies here
        energies = {}
        # Do a relaxation and energy calculation
        if experiment["Force Field"] == "CHARMM":
            CHARMM.Relaxation(molecule, experiment)
            energies["Complex"] = CHARMM.Energy(molecule, experiment)
        else:
            text = "The initialize molecule function does not support the "
            text += str(experiment["Force Field"]) + " force field."
            raise IPRO_Error(text)
        # Format the energy
        text = SHARING.format_energies(energies)
        experiment["Summary"] += text
        # Create the summary file for the Molecule
        experiment["Summary"] += "Ended" + SHARING.time_stamp()
        name = SHARING.summary_name(SHARING.get_current())
        f = open(name, "w")
        f.write(experiment["Summary"])
        f.close()
        # Leave this folder, then start sharing
        os.chdir('../')
        SHARING.Start(experiment)
        # Output the Molecule's structure
        name = "./Current/" + molecule.generate_name()
        molecule.output(name, molecule.fileFormat, experiment["User"])
        # Output the Complex energy
        f = open("./Current/Molecule" + mn + "_Energy.txt", "w")
        f.write(text)
        f.close()
        # End sharing
        SHARING.End(experiment)

def start_initialization(experiment):
    """Make the folder to carry out an initialization in"""
    # This is the folder
    folder = "initialize"
    try:
        os.mkdir(folder)
    except OSError:
        pass
    try:
        os.mkdir(folder + "/Current")
    except OSError:
        pass
    # That's it. This is really very simple

def finish_initialization(experiment):
    """Finish an initialization"""
    # Sharing is already started and the structures and energies have already
    # been loaded into the Experiment
    # Make the initial experiment summary
    experiment["Summary"] = "Initial Calculations\n"
    # Determine when they started
    f = open(experiment["Folder"]+"initialize/Group1/Group1_Summary.txt","r")
    for line in f:
        if line.startswith("Started"):
            experiment["Summary"] += line
            break
    f.close()
    # Include the energies of each Design Group
    for group in experiment:
        experiment["Summary"] += SHARING.format_energies(\
                  experiment["Energies"][group.number], group.number, False)
    # Create the initial folder in results
    SHARING.output_best(experiment, 0, None)
    # Figure out the name of the output file
    if experiment["Type"] == "Mutator":
        folder = experiment["Folder"] + "results/wildtype/"
    else:
        folder = experiment["Folder"] + "results/initial/"
    # List the energies of each Target Molecule, too
    if experiment["Energy Calculation"] == "Binding":
        for molecule in experiment[0]:
            if not molecule.design:
                name = experiment["Folder"] + "initialize/Current/Molecule"
                name += molecule.name + "_Energy.txt"
                f = open(name, "r")
                energy = f.readline().split()[2]
                f.close()
                experiment["Summary"] += "The energy of Target Molecule "
                experiment["Summary"] += molecule.name + " is " + energy
                experiment["Summary"] += " kcal / mol\n"
                # Copy the file to the output folder
                os.system("cp " + name + " " + folder)
    # Move to the Experiment's folder
    os.chdir(experiment["Folder"])
    # Try to make the current folder
    try:
        os.mkdir("Current")
    except OSError:
        pass
    # Write the Structures to the Current folder
    SHARING.output_Current(experiment, "./Current/", None, 0)
    experiment["Last Update"] = 0
    # Start a refinment
    REFINEMENT.Start(experiment, 0, None)
    # Create the initial summary file
    experiment["Summary"] += "Ended" + SHARING.time_stamp()
    name = SHARING.summary_name(experiment["Folder"])
    f = open(name, "w")
    f.write(experiment["Summary"])
    f.close()
    # Remove the initialization folder
    os.system("rm -rf initialize")
    # End the sharing that was started elsewhere
    SHARING.End(experiment)

def load_initial_info(experiment, folder = None):
    """Load the initial information."""
    # Generate the folder where it should be stored
    if folder == None:
        if experiment["Type"] == "Mutator":
            folder = experiment["Folder"] + "results/wildtype/"
        else:
            folder = experiment["Folder"] + "results/initial/"
    # Try to load the information
    try:
        SHARING.update_Energies(experiment, folder)
        SHARING.update_Current(experiment, folder)
        SHARING.load_Target_Energies(experiment, folder)
        return True
    except IOError:
        return False

def INITIALIZE(experiment):
    """Initialize an Experiment"""
    # Set the Experiment's last update to -1
    experiment["Last Update"] = -1
    # Start sharing
    SHARING.Start(experiment)
    # If the initial information is avaible, be done
    have = load_initial_info(experiment)
    if have:
        SHARING.End(experiment)
        return have
    # If there isn't an initialize folder, make it
    if "initialize" not in os.listdir("./"):
        start_initialization(experiment)
    # Move into the folder, end sharing, and start trying to do the
    # initialization
    os.chdir("initialize")
    SHARING.End(experiment)
    # Go through the groups and Molecules
    for group in experiment:
        did = initialize_group(experiment, group.number)
    if experiment["Energy Calculation"] == "Binding":
        for molecule in experiment[0]:
            if not molecule.design:
                did = initialize_molecule(experiment, molecule.name)
    # Determine if the initial information is available or not
    SHARING.Start(experiment)
    have = load_initial_info(experiment, "./Current/")
    # If everything isn't available
    if not have:
        # End sharing and move to the Experiment's folder
        SHARING.End(experiment)
        os.chdir(experiment["Folder"])
        # Wait for the initialization folder to be gone
        SHARING.Wait("initialize")
        # Start sharing, get the initial info, and be done
        SHARING.Start(experiment)
        have = load_initial_info(experiment)
        SHARING.End(experiment)
        # If the information stil doesn't exist, throw an error
        if not have:
            text = "There is a problem in the INITIALIZE function"
            raise IPRO_Error(text)
    # Otherwise, just finish the initialization
    else:
        finish_initialization(experiment)
    # At this point, all structures are ready and all energies are known, so the
    # function can be done
    return have

def Wait(experiment):
    """Wait for an IPRO Experiment to be entirely completed"""
    # Keep track using this variable
    finished = False
    while not finished:
        # Periodically check the last completed iteration
        SHARING.Start(experiment)
        n = SHARING.iteration_counter(SHARING.get_current(), False)
        SHARING.End(experiment)
        # If all iterations are complete, the program doesn't need to wait any
        # more
        if n == experiment["IPRO Iterations"]:
            finished = True
        # If another processor has called for a refinement, do it
        REFINEMENT.DO(experiment)
        # If the experiment isn't done yet, take a 15 second break before
        # continuing the while loop
        if not finished:
            time.sleep(15)
