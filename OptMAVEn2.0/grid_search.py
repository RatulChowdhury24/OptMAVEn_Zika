#!/usr/bin/env python


import sys
import os
import shutil

sys.path.append('/storage/home/mfa5147/work/IPRO_Suite/modules')
#mfa# sys.path.append("/gpfs/home/jgp5130/work/zika")
import MOLECULES
from MOLECULES import MoleculeFile
from OPTMAVEN import parameter_output_antigen 
from CHARMM import Missing_Atoms

fileFolder = "/gpfs/group/cdm/IPRO_Suite/input_files/"
topFile = "top_all27_prot_na.rtf"
parFile = "par_all27_prot_na.prm"
solFile = "solvation.dat"

antibody_clash_dir = '/storage/home/mfa5147/scratch/GCN4/position/structures'

# Get starting directory.
start_dir = os.getcwd() #mfa#

# Get name of structure file from arguments.
structure_in = sys.argv[1] #mfa#

# Get name of epitope file.
epitope_file = sys.argv[2] #mfa#

f = open(epitope_file) #mfa#
for line in f:
    newline = line.replace(":", " ")
    newline.rstrip()
    strs = newline.split()
    prefix = strs[0]
    print prefix
    if not os.path.isdir(prefix):
        os.mkdir(prefix)
    os.chdir(prefix)
    print(prefix)
    for angle in range(0, 360, 20):
	print(angle)
        if not os.path.isdir(str(angle)):
            os.mkdir(str(angle))
        os.chdir(str(angle))
        if not os.path.exists(topFile):
            shutil.copyfile(fileFolder + topFile, topFile)
        if not os.path.exists(parFile):
            shutil.copyfile(fileFolder + parFile, parFile)
        if not os.path.exists(solFile):
            shutil.copyfile(fileFolder + solFile, solFile) 
	# Copy the structure of the antigen into the current directory.
	shutil.copyfile(os.path.join(start_dir, structure_in), structure_in) #mfa#
	shutil.copyfile(os.path.join(antibody_clash_dir, "MoleculeH.txt"), "MoleculeH.txt") #mfa#
	shutil.copyfile(os.path.join(antibody_clash_dir, "MoleculeK.txt"), "MoleculeK.txt") #mfa#
        """#mfa# shutil.copyfile("/gpfs/group/cdm/tong/OPTMAVEN/position_antigen/GCN4/position/structures/moleculep_postrelaxation.pdb", "moleculep_postrelaxation.pdb")
        shutil.copyfile("/gpfs/group/cdm/tong/OPTMAVEN/position_antigen/GCN4/position/structures/MoleculeH.txt", "MoleculeH.txt")
        shutil.copyfile("/gpfs/group/cdm/tong/OPTMAVEN/position_antigen/GCN4/position/structures/MoleculeK.txt", "MoleculeK.txt")
	"""
	inputFile = MoleculeFile(structure_in) #mfa#
        #mfa# inputFile = MoleculeFile("moleculep_postrelaxation.pdb")
	G = inputFile[0] #mfa#
        #mfa# G = inputFile["P"]
	print('Missing Atoms')
        Missing_Atoms(G)
        parameter_output_antigen(G, "gcn4.txt")
        #inputFile = MoleculeFile("moleculea_postrelaxation.pdb")
        #H = inputFile["A"]
        #Missing_Atoms(H)
        #parameter_output_antigen(H, "MoleculeH.txt")
        #inputFile = MoleculeFile("moleculeb_postrelaxation.pdb")
        #L = inputFile["B"]
        #Missing_Atoms(L)
        #parameter_output_antigen(L, "MoleculeL.txt")
        fe = open("epitope.txt", "w")
	print(' '.join([epi for epi in strs[1:]])) #mfa#
        for epi in strs[1:]:
            #mfa# print epi
            fe.write(epi + "\n")
        fe.close()
        current = os.getcwd()
        cmd = "/gpfs/group/cdm/tong/jordan/IPRO_Suite_all_probability/modules/CPP/initialization/initialantigen_clash.out " + "gcn4.txt " + "epitope.txt " + "MoleculeH.txt " + "MoleculeK.txt " + str(angle)
        runScriptFile = "gcn4_" + str(angle)
	#mfa#
        text = """#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf
set -u
cd %s

%s """ %(current, cmd)

        fs = open(runScriptFile, 'w')
        fs.write(text)
        fs.close()
        cmd = "chmod +x " + runScriptFile
        os.system(cmd)
        os.chdir("../")
f.close()
