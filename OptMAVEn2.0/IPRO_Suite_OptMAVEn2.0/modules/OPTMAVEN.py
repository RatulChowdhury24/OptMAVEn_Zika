#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite OptMAVEn Functions"
# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions that are needed for running the Optimal Method of
Antibody variable Region Engineering program."""

# Import standard PYTHON modules
import os
import sys
import math
import random
# And include IPRO Suite modules
from STANDARDS import *
import SHARING
import MOLECULES
import EXPERIMENT
import ROTAMERS
import IPRO_FUNCTIONS
import CHARMM
from CPLEX import OptMAVEn_selector
from IPRO_FUNCTIONS import generate_random_angle
from DEIMMUNIZATION import load_human_sequences 
#import DOCKING_FUNCTION
#import CPLEX
#import DEIMMUNIZATION.py

class OptMAVEnError(IPRO_Error):
    """An Error for problems running OptMAVEn"""
    def __init__(self, error = ''):
        """The initialization of the OptMAVEnError class"""
        IPRO_Error.__init__(self, error)

def modify_parameters(molecule):
    """Modify force field parameters specifically for Molecule positioning"""
    # For now, just divide each VDW radius by two
    for residue in molecule:
        for atom in residue:
            if atom.vdwRadius != None:
                atom.vdwRadius /= 2

def prepare_MAPs_parts(experiment = {}):
    """Get the MAPs parts ready for use in OptMAVEn"""
    # Determine which domains to use
    domains = []
    if "OptMAVEn Domains" not in experiment:
        domains.extend(["H", "K", "L"])
    else:
        if "Heavy" in experiment["OptMAVEn Domains"]:
            domains.append("H")
        if "Light" in experiment["OptMAVEn Domains"]:
            domains.extend(["K", "L"])
    # Store the names of the files made by this function
    fileNames = []
    # Write out the types of domains being designed
    text = ''
    for d in domains:
        text += d + " "
    f = open("domains.txt", "w")
    f.write(text[:-1])
    f.close()
    fileNames.append("domains.txt")
    # Get file format and force field information
    if "Force Field" in experiment:
        ff1 = experiment["Force Field"]
    else:
        ff1 = defaultField
    if "File Format" in experiment:
        ff2 = experiment["File Format"]
    else:
        ff2 = defaultFormat
    # Go through those domains and the regions of structure
    for D in domains:
        for R in ["V", "CDR3", "J"]:
            # Store the text of the formatted structures here
            text = ''
            # Access the MAPs parts
            path = InstallFolder + "databases/MAPs/" + D + R + "/"
            # Loop through possible number options
            for i in range(1, 501):
                name = path + D + R + "_" + str(i) + ".pdb"
                try:
                    f = open(name, "r")
                except IOError:
                    # If the file couldn't be opened, we're out of parts
                    break
                # Read in the part
                lines = f.readlines()
                f.close()
                # Make a Molecule out of it
                part = MOLECULES.Molecule(lines, None, None, D, True, ff1, ff2)
                # parameterize the part
                ROTAMERS.parameterize(part, experiment)
                # Modify its parameters
                modify_parameters(part)
                # Get the formatted text, listing the rotamer number as i
                for residue in part:
                    text += format(residue, "all - " + str(i))
            # Write the formatted text of this file to a folder
            f = open(D + R + ".txt", "w")
            f.write(text)
            f.close()
            # Store the file name
            fileNames.append(D + R + ".txt")
    # Return the names of the files
    return fileNames

def load_integer_cuts():
    """Load the data for MAPs integer cuts"""
    # The information is stored here
    name = InstallFolder + "databases/OptMAVEn/MAPs_Integer_Cuts.txt"
    try:
        f = open(name, "r")
    except IOError:
        text = "The MAPs Integer Cut information could not be found"
        raise OptMAVEnError(text)
    # Store the cuts here
    cuts = []
    for line in f:
        items = line.split()
        if len(items) == 4:
            cuts.append(items)
    f.close()
    return cuts

def load_epitope_from_groups(experiment, dockingGroups):
    """ Return epitope residues from the docking groups"""
    epitope = []
    if "epitope" not in experiment.keys():
        text = "The epitope residues are not defined"
        raise OptMAVEnError(text)    
    else:
        # The epitope defined in the experiment should be a dict of dict.
        # The molecule id is the key of the inner dict
        # Get the molecule ID
        for epimol in experiment["epitope"]:
            # Get the epitope residue number for each molecule 
            for residueID in experiment["epitope"][epimol]:
                #Get the residue from the molecule
                for dockingGroup in dockingGroups:
                    for molID in dockingGroup:
                        if molID == epimol:
                            molecule = experiment[molID[0]][molID[1]]
                            epitope.append(molecule[residueID])
    
    return epitope

def load_epitope_from_molecules(experiment, antigens):
    """ Returen epitope residues from a molecule list"""
    epitope = []
    
    if "epitope" not in experiment.keys():
        text = "The epitope residues are not defined"
        raise OptMAVEnError(text)    
    else:
        # The epitope defined in the experiment should be a dict of dict.
        # The molecule id is the key of the inner dict
        # Get the molecule ID
        for epimol in experiment["epitope"]:
            # Get the epitope residue number for each molecule 
            for residueID in experiment["epitope"][epimol]:
                #Get the residue from the molecule
                for mol in antigens:
                        if mol.name == epimol:
                            epitope.append(mol[residueID])

    return epitope

def calculate_groups_center(experiment, groups):
    """ Calculate the coordinate center of the groups"""
    # Define a three elements list to store the average coordinates
    average = [0, 0, 0]
    # Count the number of atoms in the groups
    count = 0
    # Summary the coordinate for x, y , z respectively
    for group in groups:
        for mol in group:
            molecule = experiment[mol[0]][mol[1]]
            for residue in molecule:
                for atom in residue:
                    count += 1
                    for i in range(3):
                        average[i] += atom[i]
    # Get the average x, y, z coordinates
    for j in range(3):
        if count == 0:
            text = "No atoms in the groups were found"
            raise OptMAVEnError(text)
        else:     
            average[j] /= count

    return average

def calculate_molecule_center(molecule):
    """ Calculate the coordinate center of the selected molecules """
    average = [0.0, 0.0, 0.0]
    count = 0
    for residue in molecule:
        for atom in residue:
            count += 1
            for i in range(3):
                average[i] += atom[i]
    for j in range(3):
        if count == 0:
            text = "No atoms in the molecules were found"
            raise OptMAVEnError(text)
        else:     
            average[j] /= count

    return average

def move_groups_orgin(experiment, groups, average = []):
    """Move the center of the groups to the orgin"""
    #Calculate the coordinate center of the groups
    average = calculate_groups_center(experiment, groups)
    #Move all the atoms in the groups to make the center of the groups in the rogin
    for group in groups:
        for mol in group:
            molecule = experiment[mol[0]][mol[1]]
            for residue in molecule:
                for atom in residue:
                    for i in range(3):
                        atom[i] -= average[i]

def move_molecule_orgin(molecule, average = []):
    """Move the center of the molecules to the orgin"""
    #Calculate the coordinate center of the groups
    if not average:
        average = calculate_molecule_center(molecule)
    #Move all the atoms in the molecules to make the center of the molecules in the rogin
    for residue in molecule:
        for atom in residue:
            for i in range(3):
                atom[i] -= average[i]
       
def move_residues_orgin(experiment, residues, average = []):
    """ Move the center of the selected residues to the orgin """
    # Calculate the centert of the residues
    if not average:
        average = calculate_residues_center(residues) 
   
    # Store the molecules names which the residues belong to
    molName = []
    for residue in residues:
        if residues.moleculeName not in molName:
            molName.append(residues.moleculeName)
    #Move all the atoms in the groups to make the center of the groups in the rogin
    for group in groups:
        for mol in group:
            molecule = experiment[mol[0]][mol[1]]
            if molecule.name in molName:
                for residue in molecule:
                    for atom in residue:
                        for i in range(3):
                            atom[i] -= average[i]

def move_groups_original(experiment, groups, change):
    """ Move the coordinate center of the groups from the orgin to their original position"""
    #Move all the atoms in the groups to make the center of the groups in the rogin
    for group in groups:
        for mol in group:
            molecule = experiment[mol[0]][mol[1]]
            for residue in molecule:
                for atom in residue:
                    for i in range(3):
                        atom[i] += change[i]
       
#def move_residues_orginal(experiment, groups, residues, change):
#    """ Move the groups in order to center the selected residues back to their original position  """

def generate_random_number(mean, std, decimal = 3):
    """ Generate random float number with mean and std """
    rnd = round(random.gauss(mean, std), decimal)
    return rnd

def move_groups_epitope_znegative(experiment, dockingGroups):
    """ Move the groups so that the epitopes residues have the most negative Z coordinates in the antigen
    Do random roations for several iteration and make the  
    """
    # Define iteration for rotation tries
    iterations = 100
    #t_start = 3640
    #t_end = 2190
    #angle_step = 2
    #pi = 3.14
    #delta = (t_end - t_start)/(iterations - 1)
    # Define antigens list for storing the antigen molecules
    antigens = []
    # Load the epitopes 
    epitopes = load_epitope_from_groups(experiment, dockingGroups)
    # Store the molecule names that epitopes belong to
    molNames = []
    for residue in epitopes:
        if residue.moleculeName not in molNames:
            molNames.append(residue.moleculeName)
    # Store the antigen molecules in antigens
    for name in molNames:
        for group in dockingGroups:
            for molID in group:
                if name == molID:
                    antigens.append(experiment[molID[0]][name])
    # Define bestZ as the summary of all the z coordinates of the atoms of epitopes
    bestZ = 0.0
    for residue in epitopes:
        for atom in residue:
            bestZ += atom[2]
    # Start the iterations
    for i in range(iterations):
        # Copy molecues in the antigens to a new antigens list.
        antigens_new = []
        for mol in antigens:
            mol_new = mol.duplicate()
            antigens_new.append(mol_new)
        # Generate random rotation angles for x, y and z axies
        # The maximum allowed perturbation angle is 5 degree
        anglex = generate_random_number(0, 5) * pi / 180.0
        angley = generate_random_number(0, 5) * pi / 180.0
        anglez = generate_random_number(0, 5) * pi / 180.0
        # Define the x, y, z unit vector
        xvector = [1, 0, 0]
        yvector = [0, 1, 0]
        zvector = [0, 0, 1]
        # Generate the rotation matix based on the radian degree and unit vector
        xrmatrix = MOLECULES.calculate_rmatrix(anglex, xvector)
        yrmatrix = MOLECULES.calculate_rmatrix(angley, yvector)
        zrmatrix = MOLECULES.calculate_rmatrix(anglez, zvector)
        # Rotate the molecules in the new antigens 
        for mol in antigens_new:
            MOLECULES.rotate(mol, xrmatrix)
            MOLECULES.rotate(mol, yrmatrix)
            MOLECULES.rotate(mol, zrmatrix)
        # Calculate the summary of the z coordinates of the epitopes
        epitopes = load_epitope_from_molecules(experiment, antigens_new)
        sumZ = 0.0
        for residue in epitopes:
            for atom in residue:
                sumZ += atom[2]
        # If the new summary is lower than the best summary so far, keep it
        if sumZ < bestZ:
            bestZ = sumZ
            antigens = []
            for mol in antigens_new:
                mol_new = mol.duplicate()
                antigens.append(mol_new)
        else:
            continue
            #t = t_start - i * delta
            #diff = sumZ - bestZ 
            #rnd = generate_random_number(0, 1.0) 
            # Simulated annealing to determinze whether keep the roataion
            #if math.exp(-diff/t) > rnd:
            #    bestZ = sumZ
            #    antigens = []
            #    for mol in antigens_new:
            #        antigens.append(mol)

    return antigens

def move_antigen_epitope_znegative(antigen, epitopeResiduesNames):
    """ Move the antigen so that the epitope residues have the most negative Z coordinates in the antigen
    Do random roations for several iteration and make the summary of the Z coordinates smallest  
    """
    # Define iteration for rotation tries
    iterations = 100
    pi = 3.142
    # Extract the epitope residues from the antigen
    epitope = []
    for residue in antigen:
        if residue.name in epitopeResiduesNames:
            epitope.append(residue) 
    # Define bestZ as the summary of all the z coordinates of the atoms of epitopes
    bestZ = 0.0
    for residue in epitope:
        for atom in residue:
            bestZ += atom[2]
        
    # Initialize the best antigen from antigen
    antigen_best = antigen.duplicate()
    
    # Start the iterations
    for i in range(iterations):
        # Generate random rotation angles for x, y and z axies
        # The maximum allowed perturbation angle is 5 degree
        anglex = generate_random_number(0, 5) * pi / 180.0
        angley = generate_random_number(0, 5) * pi / 180.0
        anglez = generate_random_number(0, 5) * pi / 180.0
        # Define the x, y, z unit vector
        xvector = [1, 0, 0]
        yvector = [0, 1, 0]
        zvector = [0, 0, 1]
        # Generate the rotation matix based on the radian degree and unit vector
        xrmatrix = MOLECULES.calculate_rmatrix(anglex, xvector)
        yrmatrix = MOLECULES.calculate_rmatrix(angley, yvector)
        zrmatrix = MOLECULES.calculate_rmatrix(anglez, zvector)
        # Rotate the molecules in the new antigens 
        antigen_new = antigen_best.duplicate()
        MOLECULES.rotate(antigen_new, xrmatrix)
        MOLECULES.rotate(antigen_new, yrmatrix)
        MOLECULES.rotate(antigen_new, zrmatrix)
        # Calculate the summary of the z coordinates of the epitopes from the new rotated antigen
        sumZ = 0.0
        for residue in antigen_new:
            if residue.name in epitopeResiduesNames:
                for atom in residue:
                    sumZ += atom[2]
        #print "bestZ: ", bestZ, "sumZ: ", sumZ
        # If the new summary is lower than the best summary so far, keep it
        if sumZ < bestZ:
            bestZ = sumZ
            antigen_best = antigen_new.duplicate()
        else:
            continue

    return antigen_best

def move_antigen_epitope_zpositive(antigen, epitopeResiduesNames):
    """ Move the antigen so that the epitope residues have the most positive Z coordinates in the antigen
    Do random roations for several iteration and make the summary of the Z coordinates largest  
    """
    # Define iteration for rotation tries
    iterations = 100
    pi = 3.142
    # Extract the epitope residues from the antigen
    epitope = []
    for residue in antigen:
        if residue.name in epitopeResiduesNames:
            epitope.append(residue) 
    # Define bestZ as the summary of all the z coordinates of the atoms of epitopes
    bestZ = 0.0
    for residue in epitope:
        for atom in residue:
            bestZ += atom[2]
        
    # Initialize the best antigen from antigen
    antigen_best = antigen.duplicate()
    
    # Start the iterations
    for i in range(iterations):
        # Generate random rotation angles for x, y and z axies
        # The maximum allowed perturbation angle is 5 degree
        anglex = generate_random_number(0, 5) * pi / 180.0
        angley = generate_random_number(0, 5) * pi / 180.0
        anglez = generate_random_number(0, 5) * pi / 180.0
        # Define the x, y, z unit vector
        xvector = [1, 0, 0]
        yvector = [0, 1, 0]
        zvector = [0, 0, 1]
        # Generate the rotation matix based on the radian degree and unit vector
        xrmatrix = MOLECULES.calculate_rmatrix(anglex, xvector)
        yrmatrix = MOLECULES.calculate_rmatrix(angley, yvector)
        zrmatrix = MOLECULES.calculate_rmatrix(anglez, zvector)
        # Rotate the molecules in the new antigens 
        antigen_new = antigen_best.duplicate()
        MOLECULES.rotate(antigen_new, xrmatrix)
        MOLECULES.rotate(antigen_new, yrmatrix)
        MOLECULES.rotate(antigen_new, zrmatrix)
        # Calculate the summary of the z coordinates of the epitopes from the new rotated antigen
        sumZ = 0.0
        for residue in antigen_new:
            if residue.name in epitopeResiduesNames:
                for atom in residue:
                    sumZ += atom[2]
        # If the new summary is lower than the best summary so far, keep it
        if sumZ > bestZ:
            bestZ = sumZ
            antigen_best = antigen_new.duplicate()
        else:
            continue

    return antigen_best

def parameter_antigen_experiment(experiment, dockingGroups, fileNames):
    """ Prepare the antigens for use in OptMAVEn
        Modify the residue permisson and parameter them 
        Output all the information to a "Antigen.txt"
    """
    # Loop through the docking groups
    for dockingGroup in dockingGroups:
        # Go through each molecule
        for molID in dockingGroup:
            # Get that Molecule
            molecule = experiment[molID[0]][molID[1]]
            # Modify its Residues' permissions so that all Atoms will be output
            for residue in molecule:
                residue.permission = "FIXED"
            #Parameterize the Molecule
            ROTAMERS.parameterize(molecule, experiment)
            # Get the energy calculation - ready text
            contents = format(molecule, "energy")
            # Put those contents in a file
            #fileName = "Antigen_" + str(molID[0] + str(molID[1]) + ".txt"
            fileName = "Antigen.txt"
            f = open(fileName, "a")
            f.write(contents)
            f.close()
            # Store that file name so it can be deleted later
            fileNames.append(fileName)

def parameter_output_antigen(antigen, outFile = "Antigen.txt"):
    """ Parameter the antigen molecule and output the antigen to the Antigen.txt     """
    #Missing_Atoms(antigen)
    for residue in antigen:
        #residue.permission = "FIXED"
        ROTAMERS.parameterize_Residue(residue)
    modify_parameters(antigen)
    # Get the energy calculation - ready text
    contents = format(antigen, "energy")
    # Put those contents in a file
    #outFile = "Antigen.txt"
    f = open(outFile, "w")
    f.write(contents)
    f.close()
    

def parameter_output_antigen_normal(antigen, outFile = "Antigen.txt"):
    """ Parameter the antigen molecule and output the antigen to the Antigen.txt     """
    #Missing_Atoms(antigen)
    for residue in antigen:
        #residue.permission = "FIXED"
        ROTAMERS.parameterize_Residue(residue)
    # Get the energy calculation - ready text
    contents = format(antigen, "energy")
    # Put those contents in a file
    #outFile = "Antigen.txt"
    f = open(outFile, "w")
    f.write(contents)
    f.close()

def count_molecule_atom(molecule):
    """ Count the total number of atoms in a molecule """
    count = 0
    for residue in molecule:
        for atom in residue:
            count += 1
   
    return count

def measure_molecule_overlap(receptor, ligand, cutoff):
    """ Measure the overlap between two molecules, return the overlap percentage 
        Cutoff is used to determinze whether two atoms are overlapped based on the their distance 
    """
    # Count the total number of atoms in the receptor and ligand
    count_receptor = count_molecule_atom(receptor)
    count_ligand = count_molecule_atom(ligand)
    # Store the lowest number as length
    #length = min(count_receptor, count_ligand)
    length = count_ligand
    # The number of overlap    
    overlap = 0
    for receptor_residue in receptor:
        for receptor_atom in receptor_residue:
            for ligand_residue in ligand:
                for ligand_atom in ligand_residue:
                    # Calculate the distance between the two atoms
                    distance = MOLECULES.calculate_distance(receptor_atom, ligand_atom)
                    # If the distance < cutoff, indicates these two atoms are overlapped 
                    if distance < cutoff:
                        overlap += 1
    if length == 0:
        text = " The receptor or ligand has no any atoms "
        raise OptMAVEnError(text)
    else:
        return overlap/length

def initial_antigen_positions(antigen, epitopeResiduesNames):
    """ Generate several antigen initial positions relative to the antibody """
    antigens = []
    # Move the coordinate center of the antigen to the orgin
    move_molecule_orgin(antigen)
    # Move the antigen to make the epitope having the most negative Z coordinates
    antigen_new = move_antigen_epitope_znegative(antigen, epitopeResiduesNames)
    # Move the new antigen to make the center of the epitope in the orgin
    epitope = []
    for residue in antigen_new:
        if residue.name in epitopeResiduesNames:
            epitope.append(residue)
    average = calculate_molecule_center(epitope)
    move_molecule_orgin(antigen_new, average)
    # Define the translations (angstrom) along X, Y, and Z axies
    # Z axis
    zmin = -14
    zmax = -6
    # Y axis
    ymin = -4
    ymax = 4
    # x axis
    xmin = -4
    xmax = 4
    # Move the antigen along the X, Y and Z axis
    for zmove in range(zmin, zmax + 1):
        for ymove in range(ymin, ymax +1):
            for xmove in range(xmin, xmax + 1):
                coor = [xmove, ymove, zmove] 
                antigen_current = antigen_new.duplicate()
                MOLECULES.move(antigen_current, coor)
                antigens.append(antigen_current)

    return antigens

def molecule_overlap_filter(antibodyParts, antigens, maxOverlap):
    """ Filter the antigen list and remove those with high overlaped with antibody parts """ 
    antigens_filtered = []
    for antigen in antigens:
       # Measue the overlap percentage between the antibody and antigen
       percentage = 0
       for part in antibodyParts:
           percentage += measure_molcule_overlap(part, antigen)
           # If there are two many overlaps between antibody and antigen, don't keep that antigen position
           if percentage > maxOverlap:
               pass
           # Store the moved antigen to antigen list
           else:
               antigens_filtered.append(antigen)

    return antigens_filtered 

def prepare_original_antigen(antigen, epitopeResiduesNames):
    """ Generate several antigen initial positions relative to the antibody """
    # Move the coordinate center of the antigen to the orgin
    move_molecule_orgin(antigen)
    # Move the antigen to make the epitope having the most negative Z coordinates
    antigen_new = move_antigen_epitope_znegative(antigen, epitopeResiduesNames)
    # Move the new antigen to make the center of the epitope in the orgin
    epitope = []
    for residue in antigen_new:
        if residue.name in epitopeResiduesNames:
            epitope.append(residue)
    average = calculate_molecule_center(epitope)
    move_molecule_orgin(antigen_new, average)
    parameter_output_antigen(antigen_new)
    
def prepare_antigen_epitope():
    """Return the antigen epitope residues name list"""
    epitopeResiduesNames = ["44", "45", "46", "47", "48"]
    f = open("epitope.txt")
    for residueName in epitopeResiduesNames:
        f.write(residueName)
    f.close()
    return epitopeResiduesNames

def prepare_antigen_translation():
    """ Prepare the antigen x, y, z translations and output to translation_xyz.txt
    """
    # Translate the antigen along z, y, x axis
    # The center of antigen epitope should be already moved to the orgin
    # z from 3.75 to 16.25 (Angstrom)  1.25 interval
    # y from -5 to 10   2.50 interval
    # x from -10 to 5   2.50 interval
    xmin = -10
    xmax = 5
    xinterval = 2.50
    xstep = int((xmax - xmin) / xinterval)
    
    ymin = -5
    ymax = 10
    yinterval = 2.50
    ystep = int((ymax - ymin) / yinterval)
    
    zmin = -16.25
    zmax = -3.75
    zinterval = 1.25
    zstep = int((zmax - zmin) / zinterval)
    
    zmove = 0.0
    ymove = 0.0
    xmove = 0.0
    
    f = open("translation_xyz.txt", "w")

    for x in range(xstep + 1):
        xmove = xmin +  x * xinterval
        for y in range(ystep + 1):
            ymove = ymin + y * yinterval
            for z in range(zstep + 1):
                zmove = zmin + z * zinterval
                text = str(xmove) + ":" + str(ymove) + ":" + str(zmove) + "\n"
                f.write(text)
                
    f.close()    

def prepare_antigen_rotation(rotation = 60.0):
    """ Prepare the antigen rotations angle (radian) along the z axis and output to rotation_angle.txt
    """
    pi = 3.142
    rotationRadian = round(pi * rotation / 180.0, 3)
    print rotationRadian
    radian = 0.0
    if rotation == 0:
        text = "Rotation is 0 and a meaningful rotation needs to be given"
        raise OptMAVEnError(text)

    f = open("rotation_angle.txt", "w")
    step = int(360 / rotation)
    for i in range(step):
        radian += rotationRadian
        f.write(str(radian) + "\n")

    f.close()
         
def prepare_antibody_antigen(antigen, experiment = {}): 
    """ Prepare the all the files for generating the initial positions of antigens using cpp code
    """
    # Prepare the antibody parts and output the domains.txt and all the parts information 
   # prepare_MAPs_parts()
    # Output the epitope residue names to epitope.txt
    #epitopeResiduesNames = prepare_antigen_epitope()
    # Output the antigen information to Antigen.txt
    # The center of the original antigen is first moved to the orgin and rotated to have the most z negative for the epitope residues, which brings the epitope toward negative z
    ####prepare_original_antigen(antigen, epitopeResiduesNames)
    # Output the antigen translation coordinates along x, y and z axies to translation_xyz.txt
    #prepare_antigen_translation()
    # Output the antigen rotation angles along the z axis to rotation_angle.txt
    #prepare_antigen_rotation()

def load_maps_energy(energyFile = "MAPs_Energies.txt"):
    """ Load the energy from the MAPS Energies file"""
    # Define a list to store part name and corresponding interaction energy between the part and the antigen.
    partsEnergies = []
    # Open the MAPS_Energies.txt whose format is:  partName  partNumber  Energy
    parts = []
    f = open(energyFile, "r")
    for line in f:
        # Read each line in the file to a list
        line = line.strip(' \t\n')
        items = line.split()
        #print len(items)
        if len(items) == 0:
            continue
        if len(items) != 2 or items[0] in parts:
            #print items[0], len(items)
            continue
        parts.append(items[0])
        part = items[0].split('_')[0]
        partNum = items[0].split('_')[1]
        data = []
        data.append(part)
        data.append(partNum)
        #print items[1]
        #if items[1][0].isalpha():
            #print items[0]
        try:
            float(items[1])
        except ValueError:
            continue
        
        data.append(float(items[1]))
        
        # The partName (HV, HCDR3, HJ, KV, KCDR3, KJ, LV, LCDR3, LV)
        # The partNumber (1, 2, 3...)
        # Store the part name and energy
        partsEnergies.append(data)
        # Store the part name and energy
    f.close()
    return partsEnergies

def select_parts_cplex(energyFile, solutionCuts = []):
    """Select an optimal combination of parts for a particular system based on interaction energy between parts and antigen"""
    # Load the integer cuts
    struCuts = load_integer_cuts()
    # Load the MAPS energies
    energies = load_maps_energy(energyFile)
    # Using cplex to get the optimal solutions. 
    # The solution is a dictionary and the key and values are partname and number respectively
    # The energy is the total interaction energy between the selected parts and antigen
    #solution, energy = OptMAVEn_selector(energies, struCuts, solutionCuts)
    solution, energy = OptMAVEn_selector(energies, struCuts, solutionCuts)

    # Store the selected parts in a dictionary
    selected_parts = []
    for part in solution:
        #print part, solution[part]
        #selected_parts[part] = solution[part]
        selected_parts.append(part)
        selected_parts.append(solution[part])

    return selected_parts, energy

def calculate_molecule_rmsd(moleculeA, moleculeB):
    """ Calculate the rmsd between two conformations of a molecule
        No superimposing before the calculation
    """
    rmsd = 0
    # All atoms rmsd calculation
    count = 0
    caA = []
    caB = []
    for residueA in moleculeA:
        for atomA in residueA:
            if atomA.name == "CA":
                caA.append(atomA)
    for residueB in moleculeB:
        for atomB in residueB:
            if atomB.name == "CA":
                caB.append(atomB)
    if len(caA) == len(caB):
        for atomA, atomB in zip(caA, caB):
                #print atomA.residueKind, atomB.residueKind
                count += 1
                rmsd += MOLECULES.calculate_distance(atomA, atomB)
    if count == 0:
        raise OptMAVEnError()
    rmsd = math.sqrt(rmsd/count)
   
    return rmsd

def merge_molecule_parts(parts):
    """ Merge all the provided antibody parts to one molecule and rename the atom sequentially
        The order of antibody part in the parts list shoulbe be as follows (HV, HCDR3, HJ) or (KV, KCDR3, KJ) or (LV, LCDR3, LJ)
        Three parts should be included to make a antibody H or L chain molecule
    """     
    if not len(parts) == 3:
        text = " The parts list should include 3 antibody part "
        raise OptMAVEnError(text)
    atomLines = []
    for part in parts: 
        for residue in part:
            for atom in residue:
                atomLines.append(format(atom, "PDB"))
    molecule = MOLECULES.Molecule(atomLines)
    molecule.renumber()

    return molecule
    
def generate_residues_list(molecule):
    """ Extract each residue name in an antibody and save them to a residue list
        which can be used to random selecting the perturbed residues for OPTMAVEN refinement
        If a residue belongs to the CDRs, this residue name is saved three times
        which will make sure CDRs residues can be bias randomly selected with three times higher probablity than
        the non-CDRs residues.
    """
    residues = []
    for residue in molecule:
            try:
                r = int(residue.name)
            except ValueError:
                r = int(residue.name[:-1])
            if r in range(27, 39) or r in range(56, 66) or r in range(105, 118):
                    N = 3
            else:
                    N = 1
            for i in range(N):
                    residues.append(residue.name)
    return residues

def select_perturbed_molecule(gn=None):
    # Random select which molecule in the group to be perturbed
    if gn == None:
        N = 1
    else:
        N = gn
    return gn

def select_perturbed_residues(molecule):
    # Random select a residue in the residues list
    spots = generate_residues_list(molecule)
    number = random.randint(0, len(spots))
    residueName = spots[number]
    # Get the index of the Residue in the Molecule
    for i in range(len(molecule)):
        if molecule[i].name == residueName:
            index = i
            break
    # Choose how many extra spots there will be
    more = random.randint(0, 2) 
    # And how many come before
    before = random.randint(0, more)
    after = more - before
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
        perturbed.append([molecule[i].name, i])
    return perturbed

def generate_perturbed_angles(molecule, perturbed):
    """Generate the random dihedral angle perturbations"""
    # Store them here
    angles = {}
    if molecule.name not in angles:
        angles[molecule.name] = {}
        # Go through the perturbed Residues
        for item in perturbed:
            residueName = item[0]
            residueIndex = item[1]
        # Get the angles
            residue = molecule[residueIndex]
            angles[molecule.name][residue.name] = {}
            if residueName != molecule[0].name:
                angles[molecule.name][residueName]["phi"] = generate_random_angle()
            if residue.name != molecule[-1].name:
                angles[molecule.name][residueName]["psi"] = generate_random_angle()
    return angles

def optmaven_backbone_perturbation(experiment, gn = None):
    
    molecule = select_perturbed_molecule(group)
    perturbed = select_perturbed_residues(molecule) 
    # Generate the random perturbed angles 
    angles = generate_perturbed_angles(molecule, perturbed)
    positions = []
    for item in perturbed:
        residueIndex = item[1]
        positions.append(residueIndex)
    # Load human 9 mers sequence
    humanSeqs = load_human_sequences()
    # Calculate the permitted kinds for each position based on comparing the mutations between the current sequence and the human 9 mers sequence database to maximize the humanization of the designed sequence
    kinds = generate_permitted_kinds(molecule, positions, database);
    
    CHARMM.Perturbation(group, angles)

def optmaven_parameter_antigen(experiment):
    # Parameter the antigen and output it to the Antigen.txt
    
    fileFolder = "/gpfs/group/cdm/IPRO_Suite/input_files/"
    topFile = "top_all27_prot_na.rtf"
    parFile = "par_all27_prot_na.prm"
    solFile = "solvation.dat"
    shutil.copyfile(fileFolder + topFile, topFile)
    shutil.copyfile(fileFolder + parFile, parFile)
    shutil.copyfile(fileFolder + solFile, solFile)
    
    antigen = experiment["molecules"][0][2]
    parameter_output_antigen(antigen)
    
    #os.remove(topFile)
    #os.remove(parFile)
    #os.remove(solFile)
       
def optmaven_output_epitope(experiment):
    # Output the epitope residues into Epitope.txt
    f = open("Epitope.txt", "w")
    for k in experiment["Epitope Positions"]:
        for name in experiment["Epitope Positions"][k]:
            f.write(name, "\n")

    f.close()

def optmaven_antigen_positions(experiment):
    cmd = "/gpfs/group/cdm/IPRO_Suite/modules/CPP/initialization/initialantigen.out " + prefix + " " + "epitope.txt"
    os.system(cmd)

def optmaven_calculate_energy(experiment):

    domainFolder = "/gpfs/home/tul12/work/OPTMAVEN/position_antigen/native_antigen/"
    gridFolder = os.getcwd()
    domainFile = "domains.txt"
    HCDR3File = "HCDR3.txt"
    HVFile = "HV.txt"
    HJFile = "HJ.txt"
    KCDR3File = "KCDR3.txt"
    KVFile = "KV.txt"
    KJFile = "KJ.txt"
    LCDR3File = "LCDR3.txt"
    LVFile = "LV.txt"
    LJFile = "LJ.txt"
    conformationFile = "conformations.txt"

    shutil.copyfile(domainFolder + domainFile, domainFile)
    shutil.copyfile(domainFolder + HCDR3File, HCDR3File)
    shutil.copyfile(domainFolder + HJFile, HJFile)
    shutil.copyfile(domainFolder + HVFile, HVFile)
    shutil.copyfile(domainFolder + KCDR3File, KCDR3File)
    shutil.copyfile(domainFolder + KJFile, KJFile)
    shutil.copyfile(domainFolder + KVFile, KVFile)
    shutil.copyfile(domainFolder + LCDR3File, LCDR3File)
    shutil.copyfile(domainFolder + LJFile, LJFile)
    shutil.copyfile(domainFolder + LVFile, LVFile)
    shutil.copyfile(domainFolder + conformationFile, conformationFile)

    f = open('conformations.txt')
    for prefix in f:
        prefix = prefix.rstrip('\n')
        prefix = prefix.rstrip()
        pdbName = prefix + ".pdb"
        current_folder_path, current_folder_name = os.path.split(os.getcwd())
        folder = gridFolder + current_folder_name + "/"
        
        inputFile = MoleculeFile(pdbName, folder)
        antigen = inputFile['A']
        outFile = prefix + ".txt" 
        #Missing_Atoms(antigen)
        parameter_output_antigen(antigen, outFile)

    cmd = "/gpfs/group/cdm/IPRO_Suite/modules/CPP/initialization/energy_grid.out " + "conformations.txt"
    i = os.system(cmd)
    os.remove(domainFile)
    os.remove(HCDR3File)
    os.remove(HVFile)
    os.remove(HJFile)
    os.remove(KCDR3File)
    os.remove(KJFile)
    os.remove(KVFile)
    os.remove(LCDR3File)
    os.remove(LJFile)
    os.remove(LVFile)
    os.remove(conformationFile)

def optmaven_optimal_parts(optimizedPartsFile):
    f = open("conformations.txt", "r")
    #fout = open("optimized_parts.txt", "w")
    for prefix in f:
        prefix = prefix.rstrip(' \t\n')
        print prefix
        energyFile = prefix + "_maps_energies.txt"  
        if os.stat(energyFile).st_size != 0:  
            parts, energy = select_parts_cplex(energyFile)
            print len(parts), energy
            fout = open(optimizedPartsFile, "a")
            for part in parts:
                fout.write(part + " ")
            fout.write(str(round(energy,2)))
            fout.write(" ")
            fout.write(prefix)
            fout.write("\n")
            fout.close()    
        else:
            text = energyFile + " does not exist or it is empty"
            continue
        print prefix
    f.close()

def sort_optimized_parts(optimizedPartsFile):
    f = open(optimizedPartsFile, "r")
    fout = open(optimizedPartsFile + "_sorted", "w")
    lines = []
    for line in f:
        line = line.strip(' \t\n\r')
        content = line.split()
        lines.append(content)
    f.close()
    print len(lines)
    lines.sort(key=lambda x: float(x[-2]))
    print len(lines)
    for line in lines:
        for item in line:
            fout.write(item + " ")
        fout.write("\n")
    fout.close()

def optmaven_select_parts(optimizedPartsFile, howMany):
    f = open(optimizedPartsFile, "r")
    lines = []
    for line in f:
        line = line.strip(' \t\n\r')
        content = line.split()
        lines.append(content)
    lines = lines.sort(key=lambda x: float(x[-2]))
    parts = []
    antigenNames = []
    for n in howMany:
        #del lines[n][-1]
        antigens.append(lines[n].pop())
        lines[n].pop()
        parts.append(lines[n])
   
    return antigenNames, parts

def optmaven_antigen_refinement(experiment, optimizedFile, howMany):
                    
    parts, antigenNames = optmaven_select_parts(optimizedFile, howMany)
    for name in antigenNames:
        cmd = "/gpfs/group/cdm/IPRO_Suite/modules/CPP/initialization/refinement.out " + antigenName
        os.system(cmd)

def read_optimized_file(optimizedFile):
    f = open(optimizedFile, 'r')
    conformations = [] 
    for line in f:
        line = line.strip(' \t\n\r')
        items = line.split()
        conformations.append(items.pop())
    f.close()

    return conformations 

def optmaven_construct_antibodys(optmizedFile, howMany):
       
    rootFolder = "/gpfs/group/cdm/IPRO_Suite/databases/MAPs/"
    selected, antigenNames = optmaven_select_parts(optimizedFile, howMany)
    antibodys = [] 
    for line in selected:
        heavy = []
        light = []
        for item in line:
            names = item.split("_")
            folder = rootFolder + names[0] + "/"
            pdbName = item + ".pdb"
            shutil.copyfile(folder + ".pdb", pdbName)
            inputFile = MoleculeFile(pdbName)
            antigen = inputFile(item[0]) 
            parameter_output_antigen(antigen, item + ".txt")
            if item.startswith("H"):
                heavy.append(item)
            elif item.startswith(("K", "L")):
                light.append(item)
        if len(heavy) != 3 or len(light) != 3:
            text = "Needs three parts for constructing the H or L chain of a antibody molecule"
            raise OptMAVEnError(text)
        # Sort the heavy and light lists and make sure the antibody parts looks like V, CDR3, J for both H and L chain
        newHeavy = ["", "", ""]
        for partName in heavy:
            if "HV" in partName:
                newHeavy[0] = partName
            elif "HCDR3" in partName:
                newHeavy[1] = partName
            elif "HJ" in partName:
                newHeavy[2] = partName
        newLight = ["", "", ""]
        for partName in light:
            if "HV" in partName:
                newLight[0] = partName
            elif "HCDR3" in partName:
                newLight[1] = partName
            elif "HJ" in partName:
                newLight[2] = partName
        parts = []
        for partName in newHeavy:
            inputFile = MoleculeFile(partName + ".pdb")
            part = inputFile("H")
            parts.append(part)
        #H = merge_molecule_parts(parts):    
        parts = []
        for partName in newLight:
            inputFile = MoleculeFile(partName + ".pdb")
            part = inputFile("L")
            parts.append(part)
        #L = merge_molecule_parts(parts):    
        antibodys.append([H, L])
        return antibodys    
                 
def optmaven_ipro_iteration(experiment, gn = None):
    """Do an iteration of IPRO."""
    # Start the iteration
    iteration, refinement = IPRO_FUNCTIONS.Start_Iteration(experiment, gn)
    print iteration, refinement
    print "finish start iteration"
    if not refinement and iteration > IPRO_FUNCTIONS.max_iterations(experiment):
        return iteration
    # Do the steps of an iteration
    if not refinement:
        refinement = optmaven_backbone_perturbation(gn)
        print "finish refinement"
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

def calculate_interaction_energy(group):
    complex1 = CHARMM.Energy(group, None, "all")
    design  = CHARMM.Energy(group, None, True)
    target = CHARMM.Energy(group, None, False)
    print complex1, design, target
    energy = complex1 - design -target

    return energy
