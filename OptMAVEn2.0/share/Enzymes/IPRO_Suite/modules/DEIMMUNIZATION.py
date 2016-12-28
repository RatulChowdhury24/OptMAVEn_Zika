#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite DEIMMINIZATION Functions"
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
import copy
# And include IPRO Suite modules
from STANDARDS import *
import SHARING
import MOLECULES
import EXPERIMENT
import ROTAMERS
import DOCKING_FUNCTIONS
import CPLEX
import shutil
class DeimmunizationError(IPRO_Error):
    """An Error for problems in DEIMMUNIZATION mocule"""
    def __init__(self, error = ''):
        """The initialization of the DeimmunizationError class"""
        IPRO_Error.__init__(self, error)

def validate_molecule_length(molecule):
    """ Validate the moleucle to make sure the molecule containing at least 9 amino acids
    """
    # Define a variable to store the allowd minimum amino acid numbers
    minAA = 9
    # Get the total numbers of residues
    length = len(molecule)
    # If the amino acids number is less than 9, raise error
    if length < minAA:
        text = "The allowed minimum amino acid number in the molecule to humanize is 9"
        raise DeimmunizationError
         
def select_perturbation_residues(molecule):
    """ Random select 1-3 residues for perturbation"""
    # Validte the length of the molecule
    validate_molecule_length(molecule)
    # Set the maximum allowed perturbation numbers
    maxNum = 3
    # Select a random number between 1 and maxNum
    ranNum = random.randint(1, maxNum)
    # Define the perturbed list to hold the perturbed residue name
    perturbed = []
    count = 0
    while count < ranNum:
         # Random select a residue number in the molecule 
         index = random.randint(1, len(molecule))
         perturbed.append([molecule.name, index, molecule[index].name])
         count += 1

    return perturbed

def validate_perturbed_residues(molecule, perturbed):
    """ Validate perturbation list and make sure the indexes of the perturbed residues are not equal and in the range of molecule length.
    """
    length = len(molecule)
    numPerturbed = len(perturbed)
    #if 

def extract_sequence(molecule, begin, gap = 8):
    """ Given the begin point, get the 9mer sequence from a molecule
    """
    # Validte the length of the molecule
    validate_molecule_length(molecule)
    # Define the sequence variable to store the amino acid sequence
    sequence = ''
    # Get the total numbers of residues
    length = len(molecule)
    # The begin point should be equal or larger than 1 and equal or less than len(molecule) - gap 
    if isinstance(begin, int) and begin >= 0 and begin <= (length - gap) : 
        for i in range(gap + 1):
            # Get the amino acid kind
            aa = molecule[begin + i].kind
            # Convert the three-letter amino acid name to one-letter character 
            # Store it in the sequence
            sequence += convertAA["PDB"][aa]  
    if len(sequence) == gap + 1 :
         return sequence
    else:
        text += "The amino acid number in the sequence is not correct"
        raise DeimmunizationError(text)

def determine_begin_end(molecule, index):
    """ Given a amino acid number index, determine the begin and end points in the molecule for extracting the 9mers sequence """
    # Validte the length of the molecule
    validate_molecule_length(molecule)
    # Validate_the index1, index2, index3, 
    # Define the constant gap for accessing the 9mers sequence from  a certain residue
    gap = 8
    # Get the total numbers of residues
    length = len(molecule)
    # Set the begin point 
    begin = index - gap
    # Set the end point
    end = index + gap
    if begin  < 0:
        begin = 0 
    if end > length:
        end = length
    
    return begin, end    

def extract_all_sequences(molecule, begin, end, gap = 8):
    """ Given the begin, end points, extract all the possible 9mers sequences from a molecule
        If all is set to True, all the 20 amino acids for the index are considered. If flase, only
        the amino acid in the original molecule for the index are considered
    """
    # Validte the length of the molecule
    validate_molecule_length(molecule)
    # Get the total numbers of residues
    length = len(molecule)
    # Get the gap between the end and begin points
    diff = end - begin
    # Define a list to store all the extracted sequences
    sequences = []
    #if isinstance(begin, int) and isinstance(end, int) and diff >= gap and begin > 0 and end < length:
    for i in range(diff - gap):
        sequence = extract_sequence(molecule, begin + i)
        sequences.append(sequence)
            
    return sequences

def load_human_sequences():
    """ Load the human 9mer sequences to a list from "Human_9mer_Sequences.txt" """
    # Define sequences list variable to store the sequences
    sequences = []
    # Open the human 9mer sequences file
    f = open("Human_9mer_Sequences.txt", "r")
    # Store each sequence to the list
    for line in f:
          sequences.append(line.strip())

    return sequences

def count_mutations_sequence(seq1, seq2):
    """ Count the numbers of mutations between two sequences"""
    # If the lengths of these two sequences are not equal, raise error 
    if len(seq1) != len(seq2) and len(seq1) == 0 and len(seq2) == 0:
        text = "The two sequences have different amino acids"
        raise DeimmunizationError(text)
    # Define the mutations
    mutCount = 0
    # Calculate the mutations
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            mutCount += 1

    return mutCount

def count_minimum_mutations(seq, database):
    """ Count the minimum mutations between a sequence and sequences from human sequence  database"""
    # Define minMutations variable
    minMutCount = 10000
    # Get the minimum mutations by comparing to everty sequence in the human sequence database
    for humanseq in database:
        mutCount = count_mutations_sequence(seq, humanseq)
        if mutCount < minMutCount:
            minMutCount = mutCount

    return minMutCount

def count_total_mutations(seqs, database):
    """ Count the total minimum mutations for a list of sequences from human sequence database """
    total = 0
    for seq in seqs:
        total += count_minimum_mutations(seq, database)
    return total

def count_total_mutations_cpp(seqs):
    """ Count the minimum mutations between a sequence and sequences from human sequence  database using cpp code"""
    folder = "/gpfs/group/cdm/IPRO_Suite/modules/CPP/humanization/"
    name = "humanization.out"
    shutil.copyfile(folder + name, name)
    cmd = "chmod a+x " + name
    os.system(cmd)
    seqFile = "sequences.txt"
    f = open(seqFile, 'w')
    for s in seqs:
        f.write(s + "\n")
    f.close()
    cmd = "./humanization.out " + seqFile
    os.system(cmd)
    countFile = "counts.txt"
    if os.path.exists(countFile):
        f = open(countFile, 'r')
        firstline = f.readline().strip(' \t\n')
        return int(firstline)
    else:
        text = "humanization.out cpp code do not give the right counts of the mutations, please check"
        raise DeimmunizationError(text) 

def extract_sequence_molecule(molecule):
    """ Extract the one letter amino acid sequence from a molecule """
    sequence = ""
    for residue in molecule:
        sequence += convertAA["PDB"][residue.kind]
    return sequence

def make_humanization_cuts(molecule, spots):
    """ Using cpp code to make the humaniation cuts """
    positions = []
    indexNameDict = {}
    for i, residue in enumerate(molecule):
        if residue.name in spots:
            positions.append(i)
            indexNameDict[i] = residue.name
    length = len(molecule)
    positionNum = len(positions)
    #print length, positionNum 
    allowedNum = 3
    
    if not positionNum-1 in range(allowedNum):
        text = "The humanization of the sequences only allows no more than 3 residues mutated at one time"
        raise DeimmunizationError(text)
    for p in positions:
        if not p in range(length): 
            text = "The position index are not correct, please check"
            raise DeimmunizationError(text)
            
    sequence = extract_sequence_molecule(molecule)
    f = open("sequence.txt", "w")
    f.write(sequence)
    f.close()

    f = open("spots.txt", "w")
    positions.sort()
    for p in positions:
        f.write(str(p) + "\n")
    f.close()
    
    # Using the huminzation.out cpp code to get the cuts
    folder = "/gpfs/group/cdm/IPRO_Suite/modules/CPP/humanization/"
    name = "humanization.out"
    shutil.copyfile(folder + name, name)
    cmd = "chmod a+x " + name
    os.system(cmd)
    
    cmd = "./humanization.out sequence.txt spots.txt"
    os.system(cmd)
    
    cutFile = "cuts.txt"
    cuts = []
    if os.path.exists(cutFile):
        f = open(cutFile, 'r')
        for line in f:
            line = line.strip(' \t\n')
            items = line.split(';')
            solutions = {}
            for item in items:
                if item:
                    item = item.strip(' \t\n')
                    final = item.split(':')
                    solutions[indexNameDict[int(final[0])]] = final[1]
            cuts.append(solutions)
    else:
        text = "humanization.out cpp code do not give the right cuts, please check"
        raise DeimmunizationError(text) 
    
    return cuts

def generate_humanization_cuts(molecule, spots, database = []):
    """ Generate the permitted amino acid kinds given specific positions in a molecule
        to make this designed molecule sequence  more humanized.
        positions list must store the residue number instead of residue name 
    """ 
    # Convert amino acids names to number
    positions = []
    indexNameDict = {}
    for i, residue in enumerate(molecule):
        if residue.name in spots:
            positions.append(i)
            indexNameDict[i] = residue.name
    length = len(molecule)
    positionNum = len(positions)
    #print length, positionNum 
    allowedNum = 3
    
    if not positionNum-1 in range(allowedNum):
        text = "The humanization of the sequences only allows no more than 3 residues mutated at one time"
        raise DeimmunizationError(text)
     
    cuts= []
    
    if positionNum == 1:
        position = positions[0]
        #print "Position in generate_humanization_cuts function", position
        if position < 0 or position > length:
            text = " The provided residue position is out of the molecule residues range"
            raise DeimmunizationError(text)
        begin, end = determine_begin_end(molecule, position)
        seqs = extract_all_sequences(molecule, begin, end)
        #print "seqs: ", len(seqs)
        #count = count_total_mutations(seqs, database)
        count = count_total_mutations_cpp(seqs)
        #print "count: ", count
        iteraction = 0
        for aa in aminoAcids["PDB"]:
            molecule_mut = molecule.duplicate()
            molecule_mut[position].kind = aa
            seqs_mut = extract_all_sequences(molecule_mut, begin, end)
            #mutCount = count_total_mutations(seqs_mut, database)
            mutCount = count_total_mutations_cpp(seqs_mut)
            #print "seqs_mut: ", len(seqs_mut)
            #print "mutCount: ", mutCount
            #print iteration
            iteration += 1
            if mutCount > count:
                solution = {}
                solution[indexNameDict[position]] = aa
                cuts.append(solution) 

    elif positionNum == 2:
        position1 = positions[0]
        position2 = positions[1]
        if  position1 < 0 or position1 > length or position2 < 0 or position2 > length:
            text = " The provided residue position is out of the molecule residues range"
            raise DeimmunizationError(text)
        if position1 > position2 :
            position1, position2 = position2, position1
        begin_all, end1 = determine_begin_end(molecule, position1)
        begin2, end_all = determine_begin_end(molecule, position2)
        seqs = extract_all_sequences(molecule, begin_all, end_all)
        #count = count_total_mutations(seqs, database)
        count = count_total_mutations_cpp(seqs)
        #print "seqs: ", len(seqs)
        #print "count: ", count
        iteration = 0
        for aa1 in aminoAcids["PDB"]:
           molecule1 = molecule.duplicate()
           molecule1[position1].kind = aa1
           for aa2 in aminoAcids["PDB"]:
                molecule2 = molecule1.duplicate()
                molecule2[position2].kind = aa2
                begin_mut_all, end_mut1 = determine_begin_end(molecule, position1)
                begin_mut2, end_mut_all = determine_begin_end(molecule, position2)
                seqs_mut = extract_all_sequences(molecule2, begin_mut_all, end_mut_all)
                #mutCount = count_total_mutations(seqs_mut, database)
                mutCount = count_total_mutations_cpp(seqs_mut)
                #print "seqs_mut: ", len(seqs_mut)
                #print "mutCount: ", mutCount
                #print iteration
                iteration += 1
                if mutCount > count:
                    solution = {}
                    solution[indexNameDict[position1]] = aa1
                    solution[indexNameDict[position2]] = aa2
                    cuts.append(solution) 
           
    elif positionNum == 3:
        position1 = positions[0]
        position2 = positions[1]
        position3 = positions[2]
        if  position1 < 0 or position1 > length or position2 < 0 or position2 > length or position3 < 0 or position3 > length:
            text = " The provided residue positions are out of the molecule residues range"
            raise DeimmunizationError(text)
        # Make sure position1 < position2 < position3
        indexs = [position1, position2, position3]
        indexs.sort()
        position1 = indexs[0]
        position2 = indexs[1]
        position3 = indexs[2]
        begin_all, end1 = determine_begin_end(molecule, position1)
        begin3, end_all = determine_begin_end(molecule, position3)
        seqs = extract_all_sequences(molecule, begin_all, end_all)
        #count = count_total_mutations(seqs, database)
        count = count_total_mutations_cpp(seqs)
        #print "seqs: ", len(seqs)
        #print "count: ", count
        iteration = 0
        for aa1 in aminoAcids["PDB"]:
           molecule1 = molecule.duplicate()
           molecule1[position1].kind = aa1
           for aa2 in aminoAcids["PDB"]:
                molecule2 = molecule1.duplicate()
                molecule2[position2].kind = aa2
                for aa3 in aminoAcids["PDB"]:
                    molecule3 = molecule2.duplicate()
                    molecule3 = molecule2.duplicate()
                    molecule3[position3].kind = aa3
                    begin_mut_all, end_mut1 = determine_begin_end(molecule, position1)
                    begin_mut3, end_mut_all = determine_begin_end(molecule, position3)
                    seqs_mut = extract_all_sequences(molecule2, begin_mut_all, end_mut_all)
                    #mutCount = count_total_mutations(seqs_mut, database)
                    mutCount = count_total_mutations_cpp(seqs_mut)
                    #print "seqs_mut: ", len(seqs_mut)
                    #print "mutCount: ", mutCount
                    #print iteration
                    iteration += 1 
                    if mutCount > count:
                        solution = {}
                        solution[indexNameDict[position1]] = aa1
                        solution[indexNameDict[position2]] = aa2
                        solution[indexNameDict[position3]] = aa3
                        cuts.append(solution) 
        
    return cuts
