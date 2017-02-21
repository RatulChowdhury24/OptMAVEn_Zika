import copy, os
import sys
import math
import numpy
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
sys.path.append("/gpfs/work/rzc158/IPRO/IPRO_Suite/modules/")
import IPRO_FUNCTIONS
from MOLECULES import MoleculeFile
from MOLECULES import calculate_distance
from ROTAMERS import parameterize


# This code assumes that the porin is aligned to the y-axis!


# Generic formula for ellipse
def myFunction(a, b, c, d, e, x, y):
    z = a + b*x + c*x*x + d*y + e*x*y - y*y
    return z

def pore_area(PDBfilepath, VDW_radius=0.8, AngleSearchSize=10, RadiusSearchSize=20, outputName=None,
PDB_chain = "", Origin = [0.0, 0.0, 0.0], User = "Ratul"):
    # Inputs
    # PDB: PDB File Name
    # PDB_chain: Only necessary if multiple chains in "PDB"
    # Origin: Approximate center of ellipse
    # outputName: All output file names minus file extention
    # VDW_radius: Value from I-CAVER
    # AngleSearchSize: 360 degrees divided by this number number of searches
    # RadiusSearchSize: % increase of VDW radius after each step
    # User: User Name
    # exportScatter: True if export scatter plot
    # exportPYMOL: True if export PYMOL figure

    # Code
    # Create the experiment dictionary
    experiment = {}
    experiment["User"] = User

    # Create the Molecule class object
    currentdir = os.getcwd()

    pdbfpath, pdbfname = os.path.split(PDBfilepath)

    os.chdir(pdbfpath)
    porin_mols = MoleculeFile(pdbfname)
    os.chdir(currentdir)
    
    if PDB_chain != "":
        porin = porin_mols[PDB_chain]
    else:
        porin = porin_mols[0]

    # Parameterize the Molecule
    VDW = {"C": 1.7, "H": 1.2, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.8, \
    "S": 1.8}
    max_Radius = 1.8
    # Store the results in these lists   
    Edge_Atoms = []
    Edge_Coordinates = []
    # Generate the search angle
    for angle in range(AngleSearchSize):
        theta = (angle*math.pi*2/AngleSearchSize)
        # Generate the search radius
        r = float(VDW_radius)
        Found = False
        # Continuously search until a clash is found
        while not Found:
            # Generate the new coordinates
            x = r*math.sin(theta)
            y = Origin[1]
            z = r*math.cos(theta)
            new = [x, y, z]
            # Go through all of the atoms
            Min = 1000.0
            min_Atom = "NONE"
            for res in porin:
                for atom in res:
                    diff = calculate_distance(atom, new) - VDW[atom.name[0]] - \
                    VDW["H"]
                    if diff < Min:
                        Min = diff
                        min_Atom = atom.residueName + "_" + atom.name
            if Min < 0.0:
                Edge_Atoms.append(min_Atom)
                Edge_Coordinates.append([round(r*math.sin(theta),2), round(Origin[1],2),
                round(r*math.cos(theta),2)])
                Found = True
                break
            r += (float(RadiusSearchSize)/100.0)*float(VDW_radius)
    
    twoD = []
    x = []
    y = []
    for coor in Edge_Coordinates:
        twoD.append([coor[0], coor[2]])
        x.append(coor[0])
        y.append(coor[2])
    # Remove any outliers
    mean =  numpy.mean(x)
    std = numpy.std(x)
    tot_rem = []
    for i in range(len(Edge_Coordinates)):
        if x[i] < (mean-2*std):
            tot_rem.append(i)
        if x[i] > (mean+2*std):
            tot_rem.append(i)
    new2D = []
    newX = []
    newY = []
    for i in range(len(Edge_Coordinates)):
        if i not in tot_rem:
            new2D.append(twoD[i])
            newX.append(x[i])
            newY.append(y[i])
    twoD = new2D
    x = newX
    y = newY
    A, B, area, r2, model, coeffs = regression(twoD)   
    

    text = "Area:   " + format(area, '.3f') + "\n"
    text += "R^2:    " + format(r2, '.3f') + "\n"
    text += "Removed " + str(len(tot_rem)) + " Atoms from Model\n"
    text += "Model:  " + model + "\n\nEdge Atoms:\n"
    for atom in Edge_Atoms:
        items = atom.split("_")
        res = "Position: " + items[0]
        ATOM = "Atom: " + items[1]
        text += res.ljust(17) + ATOM + "\n"
    text += "\n\nEdge Coordinates:\n"
    for coor in Edge_Coordinates:
        text += "(" + format(coor[0], '.3f') + ", " + format(coor[1], '.3f')
        text += ", " + format(coor[2], '.3f') + ")\n"
    if outputName is not None:
        fileName = outputName + ".txt"
        with open(fileName, "w") as file:
            file.write(text)
    return A, B

def regression(points):
    xs = []
    ys = []
    combo = []
    for point in points:
        xs.append(float(point[0]))
        ys.append(float(point[1]))
        combo.append(2.0*float(point[1])*float(point[0]))

    C2 = numpy.array(xs)
    C1 = numpy.power(C2, 2)
    C3 = numpy.array(combo)
    C4 = numpy.array(ys)
    C5 = numpy.power(C4, 2)
    C6 = numpy.ones(C2.size)
    F = numpy.zeros((C2.size,5))
    F[:,0] = C1
    F[:,1] = C2
    F[:,2] = C3
    F[:,3] = C4
    F[:,4] = C6
    Y = numpy.ones((C2.size,1))
    Y[:,0] = C5
    M = numpy.mat(F.T) * numpy.mat(F)
    B = numpy.mat(F.T) * numpy.mat(Y)
    sol = numpy.linalg.solve(M, B)
    sol0 = -1.0*float(sol.item(0))
    sol1 = -1.0*float(sol.item(1))
    sol2 = -1.0*float(sol.item(2))
    sol3 = -1.0*float(sol.item(3))
    sol4 = -1.0*float(sol.item(4))

    theta = 0.5*math.atan(2.0*sol2/(sol0-1.0))
    A = sol0*math.pow(math.cos(theta),2)
    A += 2.0*sol2*math.sin(theta)*math.cos(theta)
    A += math.pow(math.sin(theta),2)
    B = 0
    C = sol0*math.pow(math.sin(theta),2)
    C -= 2.0*sol2*math.sin(theta)*math.cos(theta)
    C += math.pow(math.cos(theta),2)
    D = sol1*math.cos(theta) + sol3*math.sin(theta)
    E = -sol1*math.sin(theta) + sol3*math.cos(theta)
    F = sol4
    a = math.sqrt((-4*F*A*C+C*D*D+A*E*E)/(4*A*C*C))
    b = math.sqrt((-4*F*A*C+C*D*D+A*E*E)/(4*A*A*C))
    
    x0= -D/(2*A)
    y0= -E/(2*C)

    area = math.pi * a * b

    # Report Error of Equation
    # Total number of values
    n = float(len(points))
    # Mean
    sum = 0.0
    for point in points:
        sum += math.pow(float(point[1]),2)
    mean = sum/n
    # SSTot
    sstot = 0.0
    for point in points:
        sstot += math.pow(math.pow(float(point[1]),2) - mean, 2)
    # SSRes
    ssres = 0.0
    for point in points:
        X = float(point[0])
        Y = float(point[1])
        yModel = 1.0*sol0*math.pow(X, 2) + 2.0*sol2*X*Y + 1.0*sol1*X
        yModel = yModel + 1.0*sol3*Y + 1.0*sol4 + Y*Y
        ssres += math.pow(yModel,2)
    # Calculate the R^2 value for the model
    r2 = 1.0 - ssres/sstot
    # Store the equation for the model
    model = "y^2 = " + format(-sol0, '.3f') + "x^2"
    if sol1 > 0:
        model += " - " + format(sol1, '.3f') + "x"
    else:
        model += " + " + format(-sol1, '.3f') + "x"
    if sol2 > 0:
        model += " - " + format(2*sol2, '.3f') + "xy"
    else:
        model += " + " + format(-2*sol2, '.3f') + "xy"
    if sol3 > 0:
        model += " - " + format(sol3, '.3f') + "y"
    else:
        model += " + " + format(-sol3, '.3f') + "y"
    if sol4 > 0:
        model += " - " + format(sol4, '.3f')
    else:
        model += " + " + format(-sol4, '.3f')
    coeffs = [-sol0, -sol1, -sol2, -sol3, -sol4]
    
    return 2*a, 2*b, area, r2, model, coeffs

if __name__ == "__main__":

    for i in range(-5,5):
        A, B = pore_area("../mutant33.pdb", VDW_radius=0.8, AngleSearchSize=10, RadiusSearchSize=20, Origin=[0,i,0])
        print "(%.3f, %.3f)\n" %(A,B)
