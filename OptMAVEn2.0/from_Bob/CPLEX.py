#!/usr/bin/env python

# Name of the file
__name__ = "IPRO Suite CPLEX Interface"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes and Dr. Tong Li of the Costas Maranas Lab in
the Chemical Engineering Department of the Pennsylvania State University.

This file contains functions for accessing the CPLEX Modeling environment to
identify optimal selections. This includes rotamers in IPRO and canonical
structures in OptCDR."""

# Include normal PYTHON modules
import os
import sys
# Include all contents from the STANDARDS module
from STANDARDS import *
# Allow access to the CPLEX solver
sys.path.append(CPLEXFolder)
import cplex

class CPLEX_Error(IPRO_Error):
    """An Error for problems in the IPRO Suite CPLEX Module"""
    def __init__(self, error = ''):
        """The initialization of the CPLEX_Error"""
        IPRO_Error.__init__(self, error)

def cut_previous_selection(model, rotamers, previous):
    """Create an integer cut to eliminate previous solutions"""
    # Use a lot of error checking to make sure the generated equations are valid
    if isinstance(previous, list) and len(previous) > 0:
        # Loop through the options in the list
        for I, solution in enumerate(previous):
            # If it isn't a dictionary of answers, skip it
            if not isinstance(solution, dict) or len(solution) != len(rotamers):
                continue
            # Make sure that each entry in the solution refers to a specific
            # rotamer in the rotamers dictionary of lists of rotamers
            problem = False
            for rn in solution:
                if not isinstance(solution[rn], int) or rn not in rotamers or \
                not 0 <= solution[rn] < len(rotamers[rn]):
                    problem = True
                    break
            if problem:
                continue
            # We've validated that there is a choice for every Residue and that
            # each choice refers to a specific rotamer, so it is possible to
            # make the restraint. Create a list of variables and coefficients
            vars = []
            coefs = []
            # Go through the Residues in the solution
            for rn in solution:
                # Get the kind of the relevant rotamer
                kind = rotamers[rn][solution[rn]].kind
                # Go through the rotamers that match this kind and include them
                # in the sum by storing variable names in vars and coefficients
                # of 1 in coefs
                for i, rot in enumerate(rotamers[rn]):
                    if rot.kind == kind:
                        vars.append("X_" + str(rn) + "_" + str(i+1))
                        coefs.append(1)
            # Add an integer cut to the model to eliminate the possibility of
            # this solution being selected again
            model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, \
                                coefs)], senses = ["L"], rhs = [len(solution)-1])

def match_sequences(model, rotamers, experiment):
    """Make sure Dimer sequences match"""
    # Use a try statement to make sure the experiment input contains information
    # about pairs of molecules that are required to have the same sequence. The
    # try statement is necessary in case the experiment object cannot be acted
    # upon by the 'in' operator
    flag = False
    try:
        if "Dimers" in experiment:
            flag = True
    except (KeyError, TypeError, AttributeError):
        pass
    # If there is information in the experiment input about pairs of molecules
    # that must have the same sequence
    if flag:
        # Separate the rotamers by Molecule name, Residue name, and amino acid
        # kind
        sorted = {}
        for N in rotamers:
            for i, rot in enumerate(rotamers[N]):
                if rot.moleculeName not in sorted:
                    sorted[rot.moleculeName] = {}
                if rot.name not in sorted[rot.moleculeName]:
                    sorted[rot.moleculeName][rot.name] = {}
                if rot.kind not in sorted[rot.moleculeName][rot.name]:
                    sorted[rot.moleculeName][rot.name][rot.kind] = []
                # Store the rotamer's position number and rotamer number
                sorted[rot.moleculeName][rot.name][rot.kind].append([N, i+1])
        # Now go through the Dimers
        for pair in experiment["Dimers"]:
            # If there are no rotamers for one of the Molecules, it can be
            # skipped
            if pair[0] not in sorted or pair[1] not in sorted:
                continue
            # Only create equations when there are more than one kind of amino
            # acid allowed for each of the Residues
            for rn in sorted[pair[0]]:
                if rn not in sorted[pair[1]] or len(sorted[pair[0]][rn]) < 2:
                    continue
                # go through each kind of amino acid. We know they match up, but
                # check that anyway
                for kind in sorted[pair[0]][rn]:
                    if kind not in sorted[pair[1]][rn]:
                        continue
                    # Store the created variables and coefficients
                    vars = []
                    coefs = []
                    # Include all of the rotamers from the first Molecule
                    for data in sorted[pair[0]][rn][kind]:
                        vars.append("X_" + str(data[0]) + "_" + str(data[1]))
                        coefs.append(1)
                    # Include all of the rotamers from the second Molecule, but
                    # subtract them instead
                    for data in sorted[pair[1]][rn][kind]:
                        vars.append("X_" + str(data[0]) + "_" + str(data[1]))
                        coefs.append(-1)
                    # Store the linear constraint for this kind of amino acid at
                    # this pair of dimer positions in the model
                    model.linear_constraints.add(lin_expr = \
                    [cplex.SparsePair(vars, coefs)], senses = ["E"], rhs = [0])

def special_rotamer_restraints(model, rotamers, experiment = None, previous = \
                               None):
    """Create special integer cuts for use in the CPLEX model"""
    # Add any previous solutions to the model
    cut_previous_selection(model, rotamers, previous)
    # And match the sequences between dimers
    match_sequences(model, rotamers, experiment)

def make_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, previous = \
                          None):
    """Make a CPLEX model to select an optimal combination of rotamers"""
    # First, make a CPLEX variable
    model = cplex.Cplex()
    # Set it up as an MILP that minimizes a value
    model.set_problem_type(cplex.Cplex.problem_type.MILP)
    model.objective.set_sense(model.objective.sense.minimize)
    # Create the objective function, which minimizes the sum of the selected
    # rotamers with the constant portions of the system and with each other
    objV = []
    objC = []
    # Go through the positions receiving rotamers
    spots = rotamers.keys()
    spots.sort()
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            # Store the binary variable for this position
            objV.append("X_" + str(i) + "_" + str(r))
            # And its energy with the constant portions of the system. Have a
            # large maximum energy value to help make sure CPLEX treats the
            # problem as feasible
            e = RCE[i][r - 1]
            if e > 10000:
                e = 10000.0
            objC.append(e)
    # Now include the information about rotamer - rotamer energies
    # Go through the list of rotamer rotamer energies
    for data in RRE:
        # Create a 'Z' variable encoding the position / rotamer combination
        objV.append("Z_"+str(data[0]) + "_" + str(data[1]) + "_" + str(data[2])\
                    + "_" + str(data[3]))
        # Get the rotamer-rotamer energy, and limit it to a maximum value to
        # help keep the problem feasible to solve
        e = data[4]
        if e > 10000:
            e = 10000.0
        objC.append(e)
    # Put the objective function in the CPLEX model
    model.variables.add(names = objV, obj = objC, lb = [0] * len(objV), ub = \
                        [1] * len(objV), types = ['B'] * len(objV))
    # Include the restraint to select exactly one rotamer at each position
    for i in spots:
        vars = []
        for r in range(1, len(rotamers[i]) + 1):
            vars.append("X_" + str(i) + "_" + str(r))
        model.linear_constraints.add(lin_expr = [cplex.SparsePair(vars, [1] * \
                                     len(vars))], senses = ["E"], rhs = [1])
    # Create the two linearization restraints that relate Z to X
    # For each i, r, and j, the sum over s of Z(i,r,j,s) = x(i,r) ->
    # (sum over s of Z(i,r,j,s)) - x(i,r) = 0
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            for j in spots:
                # Don't double count combinations
                if j <= i:
                    continue
                # Store the X(i,r) variable with a -1 coefficient
                vars = ["X_" + str(i) + "_" + str(r)]
                coefs = [-1]
                # Loop through s (which is the rotamers at position j) and store
                # the Z values with a coefficient of 1
                for s in range(1, len(rotamers[j]) + 1):
                    vars.append("Z_" + str(i) + "_" + str(r) + "_" + str(j) + \
                                "_" + str(s))
                    coefs.append(1)
                # Store this linear constraint
                model.linear_constraints.add(lin_expr=[cplex.SparsePair(vars, \
                                             coefs)], senses = ["E"], rhs = [0])
    # for each j, s, and i, the sum over r of Z(i,r,j,s) = x(j,s). Do the same
    # things that were just done for i, r, and j
    for i in spots:
        for j in spots:
            if j <= i:
                continue
            for s in range(1, len(rotamers[j]) + 1):
                vars = ["X_" + str(j) + "_" + str(s)]
                coefs = [-1]
                for r in range(1, len(rotamers[i]) + 1):
                    vars.append("Z_" + str(i) + "_" + str(r) + "_" + str(j) + \
                                "_" + str(s))
                    coefs.append(1)
                model.linear_constraints.add(lin_expr=[cplex.SparsePair(vars, \
                                             coefs)], senses = ["E"], rhs = [0])
    # Add in any special restraints
    special_rotamer_restraints(model, rotamers, experiment, previous)
    # Return the model
    return model

def optimal_rotamer_selector(residues, rotamers, RCE, RRE, experiment = None, previous = \
                             None):
    """Use CPLEX to select an optimal combination of rotamers"""
    # Make the CPLEX model
    model = make_rotamer_selector(residues, rotamers, RCE, RRE, experiment, previous)
    # Suppress the printed output to the screen
    model.set_results_stream(None)
    # Set the time limit (seconds) for cplex optimizer
    model.parameters.timelimit.set(900)
    # Solve the model
    model.solve()
    # Get the status of the solution
    status = model.solution.get_status()
    # 101 should refer to an optimal solution and 108 to an integer solution if
    # the solver ran out of time. Other status values are not acceptable, and
    # tracing the cause of any errors may be very difficult. Good luck!
    if status not in [101, 108]:
        text = "The optimal rotamer selector function did not find an optimal "
        text += "solution. Good luck determining why."
        raise CPLEX_Error(text)
    # Get the objective value of the selector
    objective = model.solution.get_objective_value()
    # Store the rotamers that were selected
    solution = {}
    spots = rotamers.keys()
    spots.sort()
    # Go through the positions and rotamers for that position
    for i in spots:
        for r in range(1, len(rotamers[i]) + 1):
            # Get the value of the binary variable
            x = model.solution.get_values("X_" + str(i) + "_" + str(r))
            # It is a binary variable, so it should have a value of 0 or 1. So
            # if the value is > 0.01, it means the rotamer was selected
            if x > 0.01:
                # If this is the first rotamer selected at this position, store
                # it
                if i not in solution:
                    solution[i] = r - 1
                # If there are multiple rotamers for this position, raise an
                # error as that is a significant problem.
                else:
                    text = "The CPLEX solution has multiple selections for "
                    text += "Residue Number " + str(i)
                    raise CPLEX_Error(text)
    # Return the objective value and the solution
    return objective, solution
