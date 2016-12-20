#!/usr/bin/env python

# The name of the file
__name__ = "IPRO Suite Restraint Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for the creation and collection of restraints used
to maintain structures during IPRO Suite energy minimizations."""

# Include standard installation PYTHON modules
import os
import sys
# Include IPRO Suite modules, too
from STANDARDS import *
import CHARMM
import MOLECULES
import EXPERIMENT

class RestraintError(IPRO_Error):
    """An error for problems in the RESTRAINTS module"""
    def __init__(self, error = ''):
        """The initialization function of the Restraint Error class."""
        IPRO_Error.__init__(self, error)

def CHARMM_PDB_label_maker(atom, format = 1):
    """Create a label describing a PDB formatted Atom in CHARMM."""
    # Store the label here
    label = ''
    # If desired, include the Molecule's name in the specification. Technically
    # this isn't necessary (since each Residue has a unique number), but it is
    # clearer
    if format == 1:
        label += "atom ml"
        if atom.moleculeName != ' ':
            label += atom.moleculeName.lower()
        label += " "
    # Always include the Atom's Residue's number and the atom's name 
    label += str(atom.residueNumber) + " " + atom.name.lower()
    return label

def CHARMM_PDB_fix_atoms(molecules):
    """Fix Atoms in place in CHARMM when the PDB file format is used"""
    # Store the text to fix the Atoms in place here
    fix = ''
    # Loop through the molecules
    for molecule in molecules:
        # Loop through the Residues
        for residue in molecule:
            # Skip Residues where nothing is fixed in place
            if residue.freedom != "FIXED" and len(residue.fixedAtoms) == 0:
                continue
            # Since something should be fixed in place, make sure the text is
            # prepared
            if fix == '':
                fix = "defi fixed sele none end\n"
            # If the entire Residue should be fixed in place
            if residue.freedom == "FIXED":
                fix += "defi fixed sele fixed .or. resi " + str(residue.number)
                fix += " end\n"
            # Or if there are Atoms that should be fixed in place
            elif len(residue.fixedAtoms) > 0:
                # Loop through the Atom names in the list
                for an in residue.fixedAtoms:
                    # This should be impossible, but if the Atom isn't in the
                    # Residue, skip it
                    if an not in residue:
                        continue
                    # Include the Atom in the fixing
                    fix += "defi fixed sele fixed .or. "
                    fix += CHARMM_PDB_label_maker(residue[an]) + " end\n"
    # If anything is fixed in place, finish the definition of the fixed list in
    # CHARMM
    if fix != '':
        fix += "cons fix sele fixed end\n\n"
    return fix

def fix_atoms(molecules):
    """Create text to fix Atoms in place during energy minimizations"""
    # Do this differently based on the file format and force field
    if molecules[0].forceField == "CHARMM" and molecules[0].fileFormat == "PDB":
        return CHARMM_PDB_fix_atoms(molecules)
    # If the file format and force field combination aren't supported, raise an
    # error
    else:
        text = "The fix atoms function does not support the " + \
        molecules[0].fileFormat + " file format and " + molecules[0].forceField\
        + " force field combination."
        raise RestraintError(text)

def CHARMM_PDB_position(group, experiment):
    """Make restraints to keep PDB formatted Atoms near a position in CHARMM"""
    # It has already been confirmed by the position restraints function that the
    # experiment is an Experiment with position restraints and that the group is
    # a Design Group.
    # In CHARMM, these position restraints are referred to as harmonic
    # restraints
    # Store the generated restraints here
    harmonic = ''
    # Keep track of how many harmonic restraints have been made
    count = '0'
    # Loop through the harmonic restraints in the experiment
    for info in experiment["Restraints"]["Position"]:
        # Get the specifications of the restraint - the group(s), molecule(s),
        # residue(s) and atom(s) it applies to, as well as the force constant
        gn = info[0]
        mn = info[1]
        rn = info[2]
        an = info[3]
        f = format(info[4], '.3f')
        # The default is to maintain the position close to the original
        # structure. However, it IS possible to maintain the position close to
        # an Atom in another Design Group
        if len(info) == 6:
            D = info[5]
        else:
            D = 0
        # If the restraint is for a particular group that is not this one, skip
        # it
        if gn not in ['all', group.number]:
            continue
        # Also skip it if it applies to a specific Molecule that is not in this
        # Design Group
        elif mn != 'all' and mn not in group:
            continue
        # Since some Residues and Atoms may or may not be fixed in place, try to
        # make the restraint. Store the positions here
        COORS = ''
        # And the atom specifications here
        ATOMS = ''
        # Get the relevant Molecules
        if mn == 'all':
            molecules = group
        else:
            molecules = [group[mn]]
        # Loop through those relevant Molecules
        for molecule in molecules:
            # Get the relevant Residues
            if rn == 'all':
                residues = molecule
            else:
                residues = [molecule[rn]]
            # Loop through those Residues
            for residue in residues:
                # If the Residue is fixed in place, skip it
                if residue.freedom == "FIXED":
                    continue
                # Get the relevant Atoms
                atoms = []
                if an == 'all':
                    # If the Residue is an amino acid, only use the N, CA, and C
                    if residue.kind in aminoAcids[residue.fileFormat]:
                        for atomName in ['N', 'CA', 'C']:
                            atoms.append(residue[atomName])
                    # Otherwise, use all non-hydrogen atoms
                    else:
                        for atom in residue:
                            if not MOLECULES.is_hydrogen(atom):
                                atoms.append(atom)
                else:
                    atoms.append(residue[an])
                # Loop through those Atoms
                for atom in atoms:
                    # If the Atom is fixed in place, skip it
                    if atom in residue.fixedAtoms:
                        continue
                    # Otherwise the Atom should be in the restraint. If this is
                    # the first Atom, start the restraint
                    if ATOMS == '':
                        count = str(int(count) + 1)
                        ATOMS = "\ndefi rain" + count + " sele none end\n"
                    # Include the atom in the restraint
                    ATOMS += "defi rain" + count + " sele rain" + count + \
                          " .or. " + CHARMM_PDB_label_maker(atom) + " end\n"
                    # Get the Atom that this Atom is being kept close to
                    atom2 =experiment[D][molecule.name][residue.name][atom.name]
                    # Include the Atom's coordinates in COORS
                    COORS += "coor set comp xdir " + format(atom2[0], '.3f') + \
                             " ydir " + format(atom2[1], '.3f') + " zdir " + \
                             format(atom2[2], '.3f') + " sele " + \
                             CHARMM_PDB_label_maker(atom) + " end\n"
        # If a restraint was generated, include it
        if ATOMS != '':
            # Create the comparison set of Atom structures in CHARMM
            if harmonic == '':
                harmonic += "coor copy comp\n\n"
            # Constrain this set of Atoms
            harmonic += COORS + "\n" + ATOMS + "\ncons harm force " + f + \
                        "mass sele rain" + count + " end comp\n"
    # Return the harmonic restraints
    return harmonic

def position_restraints(molecules, experiment):
    """Create restraints that keep Atoms near specified positions."""
    # Make sure the input is an Experiment class object with the proper values
    # and a Design Group. Otherwise, don't do anything
    if not isinstance(experiment, (EXPERIMENT.Experiment,dict)) or "Restraints"\
    not in experiment or "Position" not in experiment["Restraints"] or not \
    isinstance(molecules, MOLECULES.DesignGroup):
        return ''
    # Create the restraints differently based on file format and force field
    if experiment["Force Field"] == "CHARMM" and experiment["File Format"] == \
    "PDB":
        return CHARMM_PDB_position(molecules, experiment)
    # Otherwise raise an error
    else:
        text = "The position_restraints function does not support the "
        text += str(experiment["File Format"]) + " file format and "
        text += str(experiment["Force Field"]) + " force field combination"
        raise RestraintError(text)

def CHARMM_PDB_distance(group, experiment):
    """Create restraints that maintain the distance between two Atoms"""
    # The distance restraints function has already made sure there are
    # restraints and that the group is a Design Group
    # In CHARMM, distance restraints are imposed using Nuclear-Overhouser
    # Effects (NOE restraints)
    # Store the restraints here
    noe = ''
    # Go through the NOE restraints in the experiment
    for info in experiment["Restraints"]["Distance"]:
        # Make sure this particular NOE restraint should be applied
        if info[0] not in [group.number, 'all']:
            continue
        # For each Atom, make sure the Molecule is actually present in this
        # Design Group AND that both Atoms aren't fixed in place
        count = 0
        for i in range(1, 3):
            # if the Molecule isn't in the group, we can stop checking
            if info[i][0] not in group:
                count = 5000
                break
            # Get the Residue
            residue = group[info[i][0]][info[i][1]]
            # If the Residue can't move or the Atom is always fixed in place,
            # increment the counter
            if residue.freedom == "FIXED" or info[i][2] in residue.fixedAtoms:
                count += 1
        # If the counter is greater than 1, the restraint shouldn't be made
        if count > 1:
            continue
        # Get the two atoms
        atom1 = group[info[1][0]][info[1][1]][info[1][2]]
        atom2 = group[info[2][0]][info[2][1]][info[2][2]]
        # Create the NOE restraint
        if noe == '':
            noe = "noe\n"
        noe += "\nassign sele " + CHARMM_PDB_label_maker(atom1) + " end sele "+\
            CHARMM_PDB_label_maker(atom2) + " end -\nkmin " + format(info[3], \
            '.3f') + " rmin " + format(info[4], '.3f') + " kmax " + \
            format(info[5], '.3f') + " rmax " + format(info[6], '.3f')+" fmax "\
            + format(info[7], '.3f') + "\n"
    # If NOE restraints were declared, end their declaration
    if noe != '':
        noe += "\nend\n\n"
    return noe

def distance_restraints(molecules, experiment):
    """Create restraints that keep two Atoms a certain distance apart"""
    # Make sure the inputs are acceptable. Otherwise, return a blank string
    if not isinstance(experiment, (EXPERIMENT.Experiment,dict)) or "Restraints"\
    not in experiment or "Distance" not in experiment["Restraints"] or not \
    isinstance(molecules, MOLECULES.DesignGroup):
        return ''
    # Create them based on file format and force field
    if experiment["Force Field"] == "CHARMM" and experiment["File Format"] == \
    "PDB":
        return CHARMM_PDB_distance(molecules, experiment)
    else:
        text = "The distance_restraints function does not support the "
        text += str(experiment["File Format"]) + " file format and "
        text += str(experiment["Force Field"]) + " force field combination"
        raise RestraintError(text)

def CHARMM_PDB_dihedral(group, experiment):
    """Create dihedral restraints for PDB formatted Atoms in CHARMM"""
    # Store the restraints here
    dihedral = ''
    # Loop through the Dihedral restraints
    for info in experiment["Restraints"]["Dihedral"]:
        # Make sure the restraint applies to this Design Group
        if info[0] not in ['all', group.number]:
            continue
        # Determine if the Restraint should actually be made
        count = 0
        for i in range(1, 5):
            # If a Molecule is missing, don't try to make the restraint
            if info[i][0] not in group:
                count = 5000
                break
            # Get the specified Residue
            residue = group[info[i][0]][info[i][1]]
            if residue.freedom == "FIXED" or info[i][2] in residue.fixedAtoms:
                count += 1
        # If the restraint should not be made, skip it
        if count > 3:
            continue
        # Create the restraint
        dihedral += "cons dihe"
        for i in range(1, 5):
            atom = group[info[i][0]][info[i][1]][info[i][2]]
            dihedral += " " + CHARMM_PDB_label_maker(atom, 2)
        dihedral += " force " + format(info[5], '.3f') + " min "
        dihedral += format(info[6], '.3f') + "\n"
    # Add an extra end line character if the restraints were actually created
    if dihedral != '':
        dihedral += "\n"
    return dihedral

def dihedral_restraints(molecules, experiment):
    """Restrain dihedral angles between Atoms."""
    # Make sure the inputs are acceptable
    if not isinstance(experiment, (EXPERIMENT.Experiment,dict)) or "Restraints"\
    not in experiment or "Dihedral" not in experiment["Restraints"] or not \
    isinstance(molecules, MOLECULES.DesignGroup):
        return ''
    # Do this differently based on the file format and force field
    if experiment["Force Field"] == "CHARMM" and experiment["File Format"] == \
    "PDB":
        return CHARMM_PDB_dihedral(molecules, experiment)
    else:
        text = "The dihedral_restraints function does not support the "
        text += str(experiment["File Format"]) + " file format and "
        text += str(experiment["Force Field"]) + " force field combination."
        raise RestraintError(text)
