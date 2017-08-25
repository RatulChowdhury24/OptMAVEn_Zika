#!/usr/bin/env python

# The name of this file
__name__ = "GAMS IPRO Suite Module"
# The documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions related to the use of GAMS for optimization
purposes in the IPRO suite of programs."""

# Include PYTHON modules
import os
import sys
# Include the STANDARDS module, which contains all "default" variables, the
# "supported" lists, the backboneAtoms and aminoAcids dictionaries, and the
# IPRO_Error class
from STANDARDS import *

# Create an error class for problems in the module
class GAMS_Error(IPRO_Error):
    """An error for problems in the GAMS module of the IPRO suite."""
    def __init__(self, error = ''):
        """The initialization of the GAMS_Error class."""
        IPRO_Error.__init__(self, error)

# Have a function that runs a GAMS script
def execute_GAMS_script(script, inputFile, outputFile):
    """Run a GAMS script."""
    # Wrap everything in a try statement to catch any errors
    try:
        # Make the input file
        f = open(inputFile, "w")
        f.write(script)
        f.close()
        # Run GAMS without printed output
        os.system(GAMSCommand + " " + inputFile + " lo=0 o /dev/null")
        # Open the output file
        f = open(outputFile, "r")
        lines = f.readlines()
        f.close()
        return lines
    # If there was an error, complain about it
    except (IOError, TypeError, OSError):
        text = "Running a GAMS script has failed. Here are the expected input "
        text += "and output file names:\nInput: " + str(inputFile)+"\nOutput: "
        text += str(outputFile)
        raise GAMS_Error(text)

def cut_previous(rotamers, previous, names, equations):
    """Use integer cuts to remove previous solutions from being possible."""
    # If the previous entry isn't a list, just return the names and equations
    # without modification
    if not isinstance(previous, list):
        return names, equations
    # Loop through the entries
    for I, group in enumerate(previous):
        # If it is not a dictionary, skip it
        if not isinstance(group, dict):
            continue
        # Make sure each entry in the dictionary corresponds to a residue
        # and a rotamer
        problem = False
        for spot in group:
            if spot not in rotamers:
                problem = True
                break
            if not isinstance(group[spot], int) or not 0 <= group[spot] < \
            len(rotamers[spot]):
                problem = True
                break
        if problem:
            continue
        # It has been confirmed that each Residue has a specific rotamer
        # designated in this solution, so use it to create integer cuts.
        # Create the name of the cut
        name = "CUT" + str(I+1)
        names += "\n\t" + name
        equation = name + " ..\n"
        # Loop through the different spots
        for spot in group:
            # identify the kind of Residue that was selected
            kind = rotamers[spot][group[spot]].kind
            # Loop through the rotamers for this Residue and include all of
            # them that are of this kind
            for i, rot in enumerate(rotamers[spot]):
                if rot.kind == kind:
                    equation += 'x("' + str(spot) + '","' + str(i+1)
                    equation += '") + '
        # Cut the last two characters off of the equation string to get rid
        # of the last +_
        equations += equation[:-2] + '=l= ' + str(len(group) - 1) + ';\n\n'
    # Return the names and equations
    return names, equations

def match_sequences(rotamers, experiment, names, equations):
    """Make sure the sequences of Dimers stay the same."""
    # Use a try statement to check the experiment input for whether or not it
    # contains a Dimers specification. The try statement is to make sure the
    # experiment is something that the 'in' function works on
    try:
        if "Dimers" not in experiment:
            return names, equations
    except (KeyError, TypeError, AttributeError):
        return names, equations
    # We know that the experiment contains information about pairs of molecules
    # that should have matching sequences, so this function may continue
    # The number of matching equations made so far
    I = 0
    # Separate the Rotamers by Molecule name, Residue name, and amino acid kind
    sorted = {}
    for i in rotamers:
        for j, r in enumerate(rotamers[i]):
            if r.moleculeName not in sorted:
                sorted[r.moleculeName] = {}
            if r.name not in sorted[r.moleculeName]:
                sorted[r.moleculeName][r.name] = {}
            if r.kind not in sorted[r.moleculeName][r.name]:
                sorted[r.moleculeName][r.name][r.kind] = []
            # Store the rotamer's position number and rotamer number
            sorted[r.moleculeName][r.name][r.kind].append([r.number, j+1])
    # Loop through the pairs of dimers
    for pair in experiment["Dimers"]:
        # If either doesn't have rotamers, skip this pair
        if pair[0] not in sorted or pair[1] not in sorted:
            continue
        # Go through all of the positions in the first Molecule
        for spot in sorted[pair[0]]:
            # Only do equations when both Residues in both Molecules allow
            # more than one kind of amino acid. Refer to the dimer_match
            # function in IPRO_FUNCTIONS for more information.
            if spot not in sorted[pair[1]] or len(sorted[pair[0]][spot]) <= 1 \
            or len(sorted[pair[1]][spot]) <= 1:
                continue
            # Go through each kind of amino acid. They should be guaranteed to
            # match
            for kind in aminoAcids[experiment["File Format"]]:
                # If this residue kind does not exist within the rotamers at
                # either position, skip the residue kind
                if kind not in sorted[pair[0]][spot].keys() and kind not in \
                sorted[pair[1]][spot].keys():
                    continue
                # Increment the equation counter
                I += 1
                # Generate a name for this equation
                name = "MATCH" + str(I)
                # Include this name in the list of equation names
                names += "\n\t" + name
                # Start creating the equation
                equation = name + " ..\n"
                # If this residue kind is absent in one set of rotamers, then it
                # should not be selected
                if kind not in sorted[pair[0]][spot].keys():
                    equation += '0 =e= '
                else:    
                    # Include the rotamers of this kind for the first Molecule
                    # on the LHS of the equation
                    for data in sorted[pair[0]][spot][kind]:
                        equation += 'x("' + str(data[0]) + '","' + \
                        str(data[1]) + '") + '
                    # Cut off the last two characters and put an equals
                    # statement
                    equation = equation[:-2] + "=e= "
                # If this residue kind is absent in the other set of rotamers,
                # then it should also not be selected
                if kind not in sorted[pair[1]][spot].keys():
                    equations += equation + '0 ;\n\n'
                else:     
                    # Put the RHS of the equation
                    for data in sorted[pair[1]][spot][kind]:
                        equation += 'x("' + str(data[0]) + '","' + \
                        str(data[1]) + '") + '
                    # Store the equation, cutting off the last two characters
                    equations += equation[:-2] + ";\n\n"
    # Return the names of the equations as well as the equations themselves
    return names, equations

# Have a function that creates special restraints for use in GAMS based on the
# type of IPRO suite program that is being run.
def special_rotamer_restraints(residues, rotamers, experiment = None, \
                               previous = None):
    """Create special restraints for use in selecting optimal rotamers."""
    # Store the names of the equations in this string
    names = ''
    # And the equation declarations themselves in this one
    equations = ''
    # Call the different special restraint functions
    names, equations = cut_previous(rotamers, previous, names, equations)
    names, equations = match_sequences(rotamers, experiment, names, equations)
    return names, equations

def make_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, \
                          previous = None):
    """Make a GAMS script to select an optimal combination of rotamers."""
    # Wrap everything in a try statement so that if there is any problem a
    # particular error message can be given
    try:
        # Retrieve any special equations
        eq_names, eqs = special_rotamer_restraints(residues, rotamers, \
                                                   experiment, previous)
        # Declare names for the input and output files for gams
        inputFile = "select_rotamers.gms"
        outputFile = "selected_rotamers.txt"
        # Identify the different positions that are receiving rotamers
        spots = rotamers.keys()
        spots.sort()
        # Identify the maximum number of rotamers at any single position
        max = 0
        for spot in spots:
            if len(rotamers[spot]) > max:
                max = len(rotamers[spot])
        # Start creating the script
        script = """* GAMS script for selecting an optimal rotamer arrangement

$INLINECOM /* */

option subsystems;

OPTIONS
\tsysout = off
\tsolprint = off
\treslim = 900
\titerlim = 1000000
\tdomlim = 0
\tlimcol = 0
\tlimrow = 0
\toptca = 0.0
\toptcr = 0.0
\twork = 50000000;

SETS
\ti\tthe positions with rotamers\t/
"""
        for spot in spots:
            script += "\t\t\t" + str(spot) + "\n"
        script += "/\n\tr\trotamers per residue /1*" + str(max) + """/;

ALIAS(i,j);
ALIAS(r,s);

PARAMETERS
\tE_rc(i,r)
\tE_rr(i,r,j,s)
\texist(i,r);

"""
        # Now its time to actually extract the Rotamer-Constant and
        # Rotamer-Rotamer energies from the RCE and RRE dictionaries,
        # respectively. Have a large maximum energy value to include, which is
        # necessary to keep GAMS from treating the problem as infeasible
        cutoff = 10000.0
        # First declare the E_rc (rotamer-constant portions) values
        for spot in spots:
            # Loop through the possible rotamers
            for i in range(max):
                # If this rotamer exists
                if i < len(RCE[spot]):
                    # If the value is greater than a cutoff (large, but needed
                    # to keep GAMS in the feasible range)
                    if RCE[spot][i] > cutoff:
                        RCE[spot][i] = cutoff
                    # Include the E_rc value
                    script += 'E_rc("' + str(spot) + '","' + str(i+1) + '") = '
                    script += format(RCE[spot][i], '.3f') + ";\n"
                    # And declare that this rotamer exists
                    script += 'exist("' + str(spot) +'","'+str(i+1)+'") = 1;\n'
                # If the rotamer doesn't exist, declare that
                else:
                    script += 'exist("'+str(spot)+'","'+str(i+1)+'") = 0;\n'
        # Now include the rotamer-rotamer energies
        for group in RRE:
            # If the energy is too large, use the cutoff value
            if group[4] > cutoff:
                group[4] = cutoff
            script += 'E_rr("' + str(group[0]) + '","' + str(group[1]) + '","'
            script += str(group[2]) + '","' + str(group[3]) + '") = '
            script += format(group[4], '.3f') + ';\n'
        # Keep declaring things in the script. I apologize that some of these
        # lines go longer than 80 characters, but that's allowed in GAMS and its
        # the easiest way to do this.
        script += """
VARIABLES
\tobj
\tz(i,r,j,s);

BINARY VARIABLES
\tx(i,r);

EQUATIONS
\tOBJECTIVE
\tROTAMER_CHOICE
\tZ_RESTRAINT1
\tZ_RESTRAINT2""" + eq_names + """;

OBJECTIVE ..
obj =e= sum((i,r)$(exist(i,r)), x(i,r) * E_rc(i,r)) + sum((i,r,j,s)$(ord(i) lt ord(j) and (exist(i,r) and exist(j,s))), z(i,r,j,s) * E_rr(i,r,j,s));

ROTAMER_CHOICE(i) ..
sum(r$(exist(i,r)), x(i,r)) =e= 1;

Z_RESTRAINT1(i,r,j)$(ord(i) lt ord(j) and exist(i,r)) ..
x(i,r) =e= sum(s$(exist(j,s)), z(i,r,j,s));

Z_RESTRAINT2(i,j,s)$(ord(i) lt ord(j) and exist(j,s)) ..
x(j,s) =e= sum(r$(exist(i,r)), z(i,r,j,s));

""" + eqs + """MODEL best_rotamers /
\tOBJECTIVE
\tROTAMER_CHOICE
\tZ_RESTRAINT1
\tZ_RESTRAINT2""" + eq_names + """/;
best_rotamers.workspace = 1500;

file results /""" + outputFile + """/;
put results;

z.lo(i,r,j,s) = 0;
z.up(i,r,j,s) = 1;

SOLVE best_rotamers USING mip MINIMIZING obj;

put best_rotamers.modelstat /;
put obj.l /;

LOOP(i,
        LOOP(r$(x.l(i,r)),
                put i.tl, r.tl /;
        );
);
"""
        # Return the script, input file, and output file
        return script, inputFile, outputFile
    # If there was any error
    except (IOError, OSError, KeyError, TypeError, ValueError, AttributeError, \
            IPRO_Error):
        text = "There was an error in the make_rotamer_selector function in "
        text += "the " + __name__ + ". Good luck determining why."
        raise GAMS_Error(text)

def optimal_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, \
                             previous = None):
    """Use GAMS to select an optimal arraingment of rotamers."""
    # Create the GAMS script
    script, inputFile, outputFile = make_rotamer_selector(residues, rotamers, \
                                    RCE, RRE, experiment, previous)
    # Run GAMS
    lines = execute_GAMS_script(script, inputFile, outputFile)
    # Store the solution here
    solution = {}
    for i, line in enumerate(lines):
        if i == 0:
            status = int(float(line) + 0.01)
        elif i == 1:
            obj = float(line)
        else:
            solution[int(line.split()[0])] = int(line.split()[1]) - 1
    # If the status is not acceptable (refer to pages 88-9 in "GAMS-A User's
    # Guide" by Richard E. Rosenthal for meanings), raise an error. 1 is an
    # optimal solution and 8 is an integer solution
    if status not in [1, 8]:
        text = "The optimal_rotamer_selection function in GAMS did not converge"
        text += " to an appropriate solution. Good luck determining why."
        raise GAMS_Error(text)
    # With the solution collected, delete the input and output files
    os.remove(inputFile)
    os.remove(outputFile)
    # Return the status, objective function value, and identified solution
    return obj, solution
