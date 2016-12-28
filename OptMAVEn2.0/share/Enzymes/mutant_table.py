import copy
import os
import sys
sys.path.append("/gpfs/work/mjg5185/IPRO_Suite/modules/")
from MOLECULES import MoleculeFile
from STANDARDS import convertAA
from STANDARDS import aminoAcids

# If you want the computer to ask you the questions, set ASK to True. If you
# want to specify them in the code, set 'ASK' to False and set the 5 variables
# in lines 14-18 
ASK = False

if not ASK:
    path1 = "./"
    path2 = path1
    file1 = "Experiment_Details.txt"
    folder = "results"

if ASK:
    question = "What is the path to the folder that contains the "
    question += "'Experiment_Details.txt' file?"
    path1 = raw_input(question + "\n")
    if path1[-1] != "/":
        path1 += "/"     

    question = "What is the path to the folder that contains the "
    question += "'results' folder? You may answer 'same' if it is "
    question += "identical to the previous folder."
    path2 = raw_input(question + "\n")
    if path2 == "same":
        path2 = path1
    if path2[-1] != "/":
        path2 += "/"

    stay = True
    while stay: 
        question = "Is the 'Experiment_Details.txt' file actually named "
        question += "'Experiment_Details.txt'?"
        file1 = raw_input(question + "\n")
        if file1 in ["yes", "Yes", "YES", "y", "Y"]:
            file1 = "Experiment_Details.txt"
            stay = False
            break
        elif file1 in ["no", "No", "NO", "n", "N"]:
            question = "What is the name of the file?"
            file1 = raw_input(question + "\n")
            stay = False
            break 
        else:
            print "This answer is not permitted. Try again."

try:
    file = open(path1 + file1)
    lines = file.readlines()
    file.close()
except IOError:
    message = "This file could not be found at the specified path."
    print message 
    raise IOError

if ASK:
    stay = True
    while stay:
        question = "Is the 'results' folder actually named "
        question += "'results'?"
        folder = raw_input(question + "\n")
        if folder in ["yes", "Yes", "YES", "y", "Y"]:
            folder = "results"
            stay = False
            break
        elif folder in ["no", "No", "NO", "n", "N"]:      
            question = "What is the name of the folder?"
            folder = raw_input(question + "\n")
            stay = False
            break
        else:
            print "This answer is not permitted. Try again."

try:
    current = os.getcwd()
    os.chdir(path2 + folder)
    os.chdir(current)
except OSError:
    message = "This folder could not be found at the specified path."
    print message
    raise OSError    


# Actual code

def gather_mutant_data(folder, design_positions, group):
    sequence = []
    ch = "NA"
    for id in design_positions:
        if id[0] != ch:
            mol = MoleculeFile("Group" + str(group) + "_Molecule" + id[0] + \
            ".pdb", "./" + folder + "/")[0]
            ch = id[0]
        AA = convertAA["PDB"][mol[id[1:]].kind]
        sequence.append(AA)
    return sequence

def gather_energy_data(folder, group):
    with open("./" + folder + "/Group" + str(group) + "_Energies.txt") as file:
        for line in file:
            items = line.split()
            if len(items) > 2:
                if items[:2] == ["Interaction", "Energy:"]:
                    IE = float(items[2])
                elif items[:2] == ["Complex", "Energy:"]:
                    CE = float(items[2])
    return IE, CE                    


def mutant_info(mutants, folder, groups, sequences, design_positions):
    sequence = gather_mutant_data(folder, design_positions, 1)
    keys = mutants.keys()
    keys.sort()
    group_keys = groups.keys()
    group_keys.sort()
    if sequence in sequences:
        for key in keys:
            if mutants[key]["sequence"] == sequence:
                break
        mutants[key]["instances"] += 1
        for group in group_keys:
            IE, CE = gather_energy_data(folder, group)
            mutants[key]["average"][group].append(IE)
            mutants[key]["complex average"][group].append(CE)
            if groups[group].upper() in ["MAINTAIN", "IMPROVE"]:
                if IE < mutants[key]["best"][group]:
                    mutants[key]["best"][group] = IE
            elif groups[group].upper() in ["REDUCE", "ELIMINATE"]:
                if IE > mutants[key]["best"][group]:
                    mutants[key]["best"][group] = IE
            else:
                print "Unexpected objective type"
                raise KeyError

    else:
        sequences.append(sequence)
        key = len(keys)
        mutants[key] = {"sequence": sequence, "best": {}, "average": {}, \
        "instances": 1, "complex average": {}}
        for group in group_keys:
            IE, CE = gather_energy_data(folder, group)
            mutants[key]["best"][group] = IE
            mutants[key]["average"][group] = [IE]
            mutants[key]["complex average"][group] = [CE]


def get_average(data, entries = None):
    if entries != None:
        if len(data) != entries:
            message = "The number of data entries does not match the " + \
            "number of mutant instances"
            print message
            raise IndexError
    divideBy = len(data)
    sum = 0.0
    for i in data:
        sum += i
    sum /= float(divideBy)
    info = "%.4f" % sum
    return info


design_positions = []
dg = {}
count = 0
for line in lines:
    items = line.split()
    if len(items) >= 7:
        if items[:2] == ["Design", "Position:"]:
            ID = items[-1] + items[3]    
            if ID in design_positions:
                print "This design position is specified twice"
                raise KeyError
            else:
                design_positions.append(ID)
        if items[:2] == ["Design", "Group:"]:
            count += 1
            dg[count] = items[2]


mutants = {}
sequences = []
os.chdir(path2 + folder)
mutant_info(mutants, "initial", dg, sequences, design_positions)
list = os.listdir(".")
List = []
for file in list:
    if file[:9] == "iteration":
        List.append(int(file[9:]))
List.sort()
for i in List:
    mutant_info(mutants, "iteration" + str(i), dg, sequences, design_positions)
os.chdir(current) 

text = "Mutant,"
for i in design_positions:
    text += i + "," 
group_keys = dg.keys()
group_keys.sort()
for group in group_keys:
    text += "Best IE (" + str(group) +"),Average IE (" + str(group) + \
    "),Average CE (" + str(group) + "),"
text += "Instances\n"
mutant_keys = mutants.keys()
mutant_keys.sort()
for key in mutant_keys:
    if key == 0:
        text += "WT,"
    else:    
        text += str(key) + ","
    for i in mutants[key]["sequence"]:
        text += i + ","
    number = mutants[key]["instances"]
    for group in group_keys:
        text += "%.4f" % mutants[key]["best"][group] + ","
        text += get_average(mutants[key]["average"][group], number) + ","
        text += get_average(mutants[key]["complex average"][group], \
        number) + ","
    text += str(number) + "\n"       


with open("Mutant_Information.csv", "w") as file:
    file.write(text)
