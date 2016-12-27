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
ASK = True 

if not ASK:
    path1 = "./"
    path2 = path1
    file1 = "Experiment_Details.txt"
    folder = "results"
    CS = "StarryNightColors"

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


schemes = ["AlpineColors", "Aquamarine", "ArmyColors", "AtlanticColors", \
"AuroraColors", "AvocadoColors", "BeachColors", "BlueGreenYellow", \
"BrassTones", "BrightBands", "BrownCyanTones", "CandyColors", \
"CherryTones", "CMYKColors", "CoffeeTones", "DarkBands", \
"DarkRainbow", "DarkTerrain", "DeepSeaColors", "FallColors", \
"FruitPunchColors", "FuchsiaTones", "GrayTones", "GrayYellowTones", \
"GreenBrownTerrain", "GreenPinkTones", "IslandColors", "LakeColors", \
"LightTemperatureMap", "LightTerrain", "MintColors", "NeonColors", \
"Pastel", "PearlColors", "PigeonTones", "PlumColors", "Rainbow", \
"RedBlueTones", "RedGreenSplit", "RoseColors", "RustTones", \
"SandyTerrain", "SiennaTones", "SolarColors", "SouthwestColors", \
"StarryNightColors", "SunsetColors", "TemperatureMap", \
"ThermometerColors", "ValentineTones", "WatermelonColors"]

if ASK:
    Flag = True
    while Flag:
        os.system("clear")
        question = "What color scheme would you like to use for plotting? "
        question += "Your choices include "
        for scheme in schemes:
            question += scheme
            if scheme == schemes[-2]:
                question += ", and "
            elif scheme == schemes[-1]:
                question += ".\n"
            else:
                question += ", "
        print question
        answer = raw_input()

        found = False
        for scheme in schemes:
            if scheme.upper() == answer.upper():
                found = True
                CS = scheme
                break

        if not found:
            print "The color scheme provided is not allowable. Try again."
        else:
            Flag = False
            break


if not ASK:
    found = False
    for scheme in schemes:
        if scheme.upper() == CS.upper():
            found = True
            CS = scheme
            break
if not found:
    print "The color scheme provided is not allowable."
    raise IOError


order = []
mutants = {}
mols = []
for line in lines:
    items = line.split()
    if len(items) >= 7:
        if items[:2] == ["Design", "Position:"]:
            ID = items[-1] + items[3]    
            if items[-1] not in mols:
                mols.append(items[-1])
            if ID in mutants.keys():
                print "This design position is specified twice"
                raise KeyError
            for aa in aminoAcids['PDB']:
                order.append(ID + "_" + convertAA["PDB"][aa])
            mutants[ID] = {}      


os.chdir(path2 + folder)
list = os.listdir(".")
count = 0
WT = {}
for i in list:
    if i[:9] == "iteration":
        count += 1
        for ch in mols:
            mol = MoleculeFile("Group1_Molecule" + ch + ".pdb", "./" + i +
            "/")[0]
            for ID in mutants.keys():
                if ID[0] == ch:
                    if count in mutants[ID]:
                        message = "The mutant was specified twice"
                        print message
                        raise KeyError
                    mutants[ID][count] = convertAA["PDB"][mol[ID[1:]].kind]
    elif i == "initial":
        for ch in mols:
            mol = MoleculeFile("Group1_Molecule" + ch + ".pdb", "./" + i +
            "/")[0]
        for ID in mutants.keys():
            if ID[0] == ch:
                WT[ID] = convertAA["PDB"][mol[ID[1:]].kind]
os.chdir(current)

MUTANTS = copy.deepcopy(mutants)
m_keys = MUTANTS.keys()

# Check against WT for copies
for i in MUTANTS[m_keys[0]].keys():
    COPY = True
    for id in m_keys:
        if MUTANTS[id][i] != WT[id]:
            COPY = False
            break
    if COPY:
        if i in mutants[m_keys[0]].keys():
            for id in m_keys:
                del mutants[id][i]

# Check that other mutants are not copies
for i in MUTANTS[m_keys[0]].keys():
    for j in MUTANTS[m_keys[0]].keys():
        if i >= j:
            continue
        COPY = True
        for id in m_keys:
            if MUTANTS[id][i] != MUTANTS[id][j]:   
                COPY = False
                break
        if COPY:
            if j in mutants[m_keys[0]].keys():
                for id in m_keys:
                    del mutants[id][j]

data = {}
COUNT = {}
for id in mutants.keys():
    for i in mutants[id].keys():
        label = id + "_" + mutants[id][i]
        if label not in data.keys():
            data[label] = {}
            COUNT[label] = 0.0

for label in data.keys():
    for label2 in data.keys():
        data[label][label2] = 0.0

for id1 in mutants.keys():
    for i in mutants[id1].keys():
        label1 = id1 + "_" + mutants[id1][i]
        COUNT[label1] += 1.0
        for id2 in mutants.keys():
            label2 = id2 + "_" + mutants[id2][i]
            data[label1][label2] += 1.0

for label1 in data.keys():
    for label2 in data.keys():
        data[label1][label2] /= COUNT[label1]
            
ORDER = []  
for label in order:
    if label in data.keys():
        ORDER.append(label)


text = """Format:
Label1,Label2,Probability

Labels are Formatted as M###_A, where M is the molecule name, ### is the
residue's position within the molecule, and A is the one-letter amino acid
abbreviation

*** The probability is calculated by dividing by the number of instances for the
mutant specified by label 1. If looking at the covariance matrix, go through the
row and each column is the probability***


"""

for label1 in data.keys():
    for label2 in data.keys():
        val = '%0.2f' % data[label1][label2] 
        text += label1 + "," + label2 + "," + val + "\n"

with open("Mutant_Covariance.txt", "w") as file:
    file.write(text)


# Write the Mathematica script
script = '<< PlotLegends`\n'
script += 'CS = "' + CS + '";\n'
script += 'matrix = {'
for i in range(len(ORDER)):
    k = len(ORDER) - i - 1
    label = ORDER[k]
    script += '{'
    for j in range(len(ORDER)):
        Label = ORDER[j]
        script += str(data[label][Label])
        if j != len(ORDER) - 1:
            script += ', '
    script += '}'
    if i == len(ORDER) - 1:
        script += '};\n'
    else:
        script += ',\n'
script += 'figure = ShowLegend[MatrixPlot[matrix, ColorFunction -> CS,'
script += 'FrameTicks -> {{'    
for i in range(len(ORDER)):
    k = len(ORDER) - i - 1
    label = ORDER[k]
    script += '{' + str(i+1) + ', "' + label + ' (' + str(int(COUNT[label])) + \
    ')"}'
    if i != len(ORDER) - 1:
        script += ', '
script += '}, {'
for i in range(len(ORDER)):
    label = ORDER[i]
    script += '{' + str(i+1) + ', Rotate["' + label + '", Pi/2]}'
    if i != len(ORDER) - 1:
        script += ', '
script += '}}, FrameLabel -> {Style["Mutant", FontSize -> 20],' 
script += ' Style["Probability", FontSize -> 20]},'
script += ' ColorFunction -> (ColorData[CS]['
script += 'Rescale[#, {0, 1}]] &), ColorFunctionScaling -> False'
script += '],{ColorData[CS][1-#] &, 10, "1", "0"'
script += ', LegendPosition -> {1.05, 0.2}}];\n'
script += 'Export["Mutant_Covariation.png", figure, ImageSize -> 2000,'
script += ' ImageResolution -> 1000];\nExit[]\n'

with open("Covariance.m","w") as file:
    file.write(script)

message = os.system("math -noprompt -script Covariance.m")
if message == 32512:
    print "Load the Mathematica module before executing"
    raise IOError

os.system("display Mutant_Covariation.png")    
