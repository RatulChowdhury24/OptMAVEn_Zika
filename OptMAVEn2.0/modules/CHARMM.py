#!/usr/bin/env python

# The name of this module
__name__ = "IPRO Suite CHARMM Module"
# The documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This module contains functions for creating and executing CHARMM scripts to
carry out various functions in the IPRO suite of programs."""

# Include PYTHON modules
import os
import sys
# Include the contents of the STANDARDS module, which contains all 'default'
# variables, the 'supported' lists, the aminoAcids and backboneAtoms
# dictionaries, and the IPRO_Error class
from STANDARDS import *
# Include the classes from the MOLECULES module
import MOLECULES
# Include the Experiment, Solvation, Restraints, and I/O modules
import EXPERIMENT
import SOLVATION
import RESTRAINTS

# Create an Error for problems in this module
class CHARMM_Error(IPRO_Error):
    """An error for problems in the CHARMM module."""
    def __init__(self, error = ''):
        """The initialization routine of the CHARMM_Error class."""
        IPRO_Error.__init__(self, error)

# Have a function to run a CHARMM script. Please note that this function will
# have to be edited by each new group to use the IPRO suite so that CHARMM works
# properly
def execute_CHARMM_script(script, procedure = None, gn = None):
    """Create and execute a CHARMM script for the IPRO suite of programs."""
    # Validate that the script is a string so it can be written to a file
    if not isinstance(script, str):
        text = "The execute_CHARMM_script requires a string as the 'script' "
        text += "input to function, not:\n" + str(script)
        raise CHARMM_Error(text)
    # Validate the procedure
    procedure = validate_procedure(procedure)
    # Make the group number input useable
    if isinstance(gn, int) and gn >= 0:
        num = str(gn)
    elif isinstance(gn, str) and ' ' not in gn and '\t' not in gn and '\n' not\
    in gn:
        num = gn
    else:
        num = ''
    # Make names for the input and output files for CHARMM
    if len(num) > 0:
        input = procedure.lower() + "_" + num + ".inp"
    else:
        input = procedure.lower() + ".inp"
    output = input[:-3] + "out"
    # Make the input file
    file = open(input, "w")
    file.write(script)
    file.close()
    # Run CHARMM. The CHARMM Command is located in the STANDARDS module
    success = os.system(CHARMMCommand + " < " + input + " > " + output)
    # Make sure CHARMM ran successfully
    if success != 0:
        text = "Running a CHARMM script has failed. Good luck determining why."
        raise CHARMM_Error(text)
    # Return the names of the input and output files
    return input, output

# This module contains 4 major functions: Missing_Atoms, Relaxation,
# Perturbation, and Energy. Each of those functions can work with the same
# inputs, which are validated here.
def validate_molecules(molecules):
    """Make sure the inputs to a CHARMM function are acceptable."""
    # The input should be a Design Group, a Molecule, or a list of Molecules. If
    # it is a Design Group, get the group number. Make sure the final result is
    # iterable and each entry is a Molecule
    error = False
    gn = None
    # molecules can be one of three things: a DesignGroup
    if isinstance(molecules, MOLECULES.DesignGroup):
        gn = molecules.number
    # A Molecule
    elif isinstance(molecules, MOLECULES.Molecule):
        molecules = [molecules]
    # Or a list of Molecules
    elif isinstance(molecules, list) and len(molecules) > 0:
        for obj in molecules:
            if not isinstance(obj, MOLECULES.Molecule):
                error = True
                break
    else:
        error = True
    # If molecules is not one of those three things
    if error:
        text = "The functions in the CHARMM module will work with either a "
        text +="Design Group, Molecule, or list of Molecules as input, but not:"
        text += "\n" + str(molecules)
        raise CHARMM_Error(text)
    # Make sure that all Molecules use the CHARMM force field and the identified
    # file format, and that each Molecule has a unique name
    names = []
    for molecule in molecules:
        # Check that it uses CHARMM
        if molecule.forceField != "CHARMM":
            text = "All Molecules used by a function in the CHARMM module must "
            text += "use the CHARMM force field."
            raise CHARMM_Error(text)
        # And they have a consistent file format
        elif molecule.fileFormat != molecules[0].fileFormat:
            text = "All Molecules used by a function in the CHARMM module must "
            text += "use the same file format."
            raise CHARMM_Error(text)
        # And that it has a unique name
        elif molecule.name in names:
            text = "All Molecules used by a function in the CHARMM module must "
            text += "have a unique name."
            raise CHARMM_Error(text)
        # Otherwise store the Molecule's name
        else:
            names.append(molecule.name)
    # Return the pointer to the molecules, the identified group number, and the
    # file format
    return molecules, gn

def validate_procedure(procedure):
    """Make sure the specified procedure can be used in naming things."""
    # If it is not a string, use "charmm"
    if not isinstance(procedure, str):
        return "charmm"
    else:
        # Split on white space and replace it with underscores
        items = procedure.split()
        procedure = ''
        for i, item in enumerate(items):
            if i != 0:
                procedure += "_"
            procedure += item
        return procedure

def introduction(purpose = None):
    """Create a header to start a CHARMM script."""
    # If a purpose has been specified, use that as the start of the file
    if isinstance(purpose, str):
        text = "* " + purpose.strip()
    # Otherwise use a generic phrase
    else:
        text = "* IPRO Suite CHARMM script"
    # Set options: warning and bomb levels to -2, and print level left at 0 but
    # optioned out (its always been that way, so we leave it for nostaligic
    # reasons)
    text += "\n*\n\nwrnl -2\n!prnl -2\nbomb -2\n\n"
    return text

def load_input_files(experiment = None):
    """Load the Topology and Parameter input files in a CHARMM script."""
    # Loop through topology and parameter data
    sets = [["Topology", defaultCHARMMTopologies, "rtf"], \
            ["Parameter", defaultCHARMMParameters, "para"]]
    # Store the output in this string
    output = ''
    for set in sets:
        # Get the list of files that will be used
        try:
            files = experiment["CHARMM " + set[0] + " Files"]
        except (KeyError, TypeError, AttributeError, IPRO_Error):
            files = set[1]
        # Make sure those file names are useable. They must be a non-empty list
        error = False
        if not isinstance(files, list) or len(files) == 0:
            error = True
        else:
            for F in files:
                # Each entry in the files must be a string that is lower case
                # and doesn't contain any white space
                if not isinstance(F, str) or ' ' in F or '\n' in F or '\t' in F\
                or not F.islower():
                    error = True
        # If the files aren't all useable
        if error:
            text = "The CHARMM force field requires that the " + set[0]
            text += " Files be a list of lowercase strings with no whitespace "
            text += "in them."
            raise CHARMM_Error(text)
        # If an error wasn't raised, create a message saying what is being done
        output += "! Load the " + set[0] + " file(s)\n"
        # Loop through the files, creating the proper content to load them in
        # CHARMM
        for file in files:
            # Open the file
            output += "open read unit 10 form name " + file + " card\n"
            output += "read " + set[2] + " card unit 10"
            # If this isn't the first file, append it to what has already been
            # read in
            if file != files[0]:
                output += " append"
            output += "\nclose unit 10\n\n"
    return output

def load_molecules(molecules, procedure, who = defaultUser, which = "all"):
    """Create text to load Molecules in a CHARMM script."""
    # This function assumes that the molecules and procedure have already been
    # validated.
    # Validate the which input
    if which not in ['all', 'ALL', True, False]:
        text = "The load_molecules function does not recognize the following "
        text += "as a valid 'which' input:\n" + str(which)
        raise CHARMM_Error(text)
    # Create the string
    text = ''
    # Make sure the Residues and Atoms are sequentially numbered
    an = 1
    rn = 1
    # Loop through the Molecules
    for molecule in molecules:
        # Only include the appropriate molecules
        if which in ['all', 'ALL', molecule.design]:
            # Renumber the Molecule
            molecule.renumber(an, rn)
            an = molecule.lastAtom + 1
            rn = molecule.lastResidue + 1
            # Generate the name of the Molecule's file
            name = molecule.generate_name("pre" + procedure, "CHARMM")
            # output the Molecule to that file
            molecule.output(name, "CHARMM", who)
            # Tell CHARMM what it is doing with this Molecule
            if molecule.design:
                text += "! Load DESIGN Molecule "
            else:
                text += "! Load TARGET Molecule "
            # Include the Molecule's name
            if not molecule.name.isspace():
                text += molecule.name + "\n"
            # If the name is all white space, put an underscore so the user
            # knows the Molecule
            else:
                text += "_\n"
            # Have different methods based on the file format of the Molecule
            # (the validate_molecules function made sure these Molecules all
            # have the same file format)
            if molecule.fileFormat == "PDB":
                # Open the Molecule
                text += "open read unit 10 form name " + name + "\n"
                # Read in its sequence and generate a setup
                text += "read sequ pdb offi unit 10\nclose unit 10\ngene ml"
                if molecule.name != ' ':
                    text += molecule.name.lower()
                # Read in the Molecule's coordinates from the PDB file
                text += " setup\nopen read unit 10 form name " + name + "\n"
                text += "read coor pdb unit 10\nclose unit 10\n\n"
            # If the file format isn't supported
            else:
                error = "The load_molecules function does not support the "
                error += str(molecule.fileFormat) + " file format."
                raise CHARMM_Error(error)
    return text

def ic_fill():
    """Tell CHARMM to include missing Atoms."""
    text = "! Add missing Atoms and assign them coordinates\n"
    text += "ic fill preserve\nic param\nic build\nhbuild\n\n"
    return text

def energy_line(experiment, procedure, h = '', n = '', c = '', s = ''):
    """Generate a string saying how energy should be calculated in CHARMM."""
    # It is assumed that the procedure has been validated
    # Try to get the list of energy terms to use from the experiment
    try:
        terms = experiment["CHARMM Energy Terms"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        terms = defaultCHARMMEnergies
    # Generate the line
    text = "skip all excl"
    for term in terms:
        text += " " + term
    # Include additional terms, as indicated by the provided strings
    if h != '' and "harm" not in terms:
        text += " harm"
    if n != '' and "noe" not in terms:
        text += " noe"
    if (procedure == "perturbation" or c != '') and "cdih" not in terms:
        text += " cdih"
    if s != '' and "gbener" not in terms:
        text += " gbener"
    text += "\n"
    return text

def minimize(molecules, experiment, procedure, gn):
    """Generate the text to run an energy minimization in CHARMM."""
    # It is assumed that the molecules and procedure have already been validated
    # Get information about how and whether solvation should be used
    solvation = SOLVATION.get_string(experiment, procedure)
    # Don't include harmonic, NOE, or CDIH restraints during backbone
    # perturbations. They'll be included in the follow up backbone relaxation,
    # so they will still be obeyed
    if procedure != "perturbation":
        harmonic = RESTRAINTS.position_restraints(molecules, experiment)
        noe = RESTRAINTS.distance_restraints(molecules, experiment)
        cdih = RESTRAINTS.dihedral_restraints(molecules, experiment)
    else:
        harmonic = ''
        noe = ''
        cdih = ''
    # Get the information about fixing atoms in place
    fix = RESTRAINTS.fix_atoms(molecules)
    # Create the text
    text = "! Carry out an energy minimization\n\n"
    text += solvation + "\nnbon nbxm 5\n" + harmonic + noe + cdih + fix
    text += energy_line(experiment, procedure, harmonic, noe, cdih, solvation)
    text += "mini abnr nstep "
    # Try to access the number of iterations, but use the default if necessary
    try:
        text += str(experiment["CHARMM Iterations"])
    except (KeyError, TypeError, AttributeError, EXPERIMENT.ExperimentError):
        text += str(defaultCHARMMIterations)
    text += " nprint 50 -\n tolgrd 0.01 tolenr 0.0001 tolstp 0.00\n\n"
    return text

def output_molecules(molecules, procedure, which = "all"):
    """Tell CHARMM to output the structures of Molecules."""
    # This function assumes the molecules and procedure have already been
    # validated
    # Validate the which input
    if which not in ['all', 'ALL', True, False]:
        text = "The output_molecules function does not recognize the following "
        text += "as a valid 'which' input:\n" + str(which)
        raise CHARMM_Error(text)
    # Store the created text in this string
    text = ''
    # Loop through the molecules
    for molecule in molecules:
        # Only consider the Molecules that are needed
        if which in ['all', 'ALL', molecule.design]:
            # Generate the name for the molecule to be output to
            name = molecule.generate_name("post" + procedure, "CHARMM")
            # Explain what is happening in CHARMM
            if molecule.design:
                text += "! Output DESIGN Molecule "
            else:
                text += "! Output TARGET Molecule "
            # Include the Molecule's name
            if not molecule.name.isspace():
                text += molecule.name + "\n"
            else:
                text += "_\n"
            # Have a different method for each supported file format
            if molecule.fileFormat == "PDB":
                # Open the file to write the structure to
                text += "open write unit 10 name " + name + " card\n"
                # Write the proper structure
                text += "write coor sele segi ml"
                if not molecule.name.isspace():
                    text += molecule.name.lower()
                text += " end pdb unit 10 card\nclose unit 10\n\n"
            # If the file format isn't supported raise an error
            else:
                text = "The output_molecules function does not support the "
                text += str(molecule.fileFormat) + " file format."
                raise CHARMM_Error(text)
    # Return the text
    return text

def load_structures(molecules, procedure, which = "all"):
    """Load the structures of Molecules after a CHARMM script."""
    # It is assumed that the molecules and procedure have been validated
    # Check the which input
    if which not in ['all', 'ALL', True, False]:
        text = "The load_structures function does not recognize " + str(which)
        text += " as a valid which input."
        raise CHARMM_Error(text)
    # Loop through the Molecules
    for mol in molecules:
        # Make sure this Molecule should be searched for
        if which not in ['all', 'ALL', mol.design]:
            continue
        # Get the name of the Molecule's structure
        name = mol.generate_name("post" + procedure, "CHARMM")
        # Try to load the structure
        try:
            file = open(name, "r")
            lines = file.readlines()
            file.close()
            mol.load(lines)
        # If there is a problem raise an error
        except IOError:
            text = "The " + name + " file could not be opened after a "
            text += procedure + " CHARMM script. Please review the CHARMM "
            text += "output file to determine why this occurred."
            raise CHARMM_Error(text)

def clean_up(molecules, input, output, procedure, energy = None):
    """Clean up files after a CHARMM procedure."""
    # It is assumed that the procedure and molecules are validated
    # Remove the Molecules' input and output structures
    for mol in molecules:
        try:
            os.remove(mol.generate_name("pre" + procedure, "CHARMM"))
        except OSError:
            pass
        # Energy functions don't make structure outputs, so don't try to remove
        # them if they were never made
        if energy != None:
            continue
        try:
            os.remove(mol.generate_name("post" + procedure, "CHARMM"))
        except OSError:
            pass
    # Move the input, output, and possibly energy files
    os.rename(input, "charmm.inp")
    os.rename(output, "charmm.out")
    if energy != None:
        os.rename(energy, "charmm_energy.txt")

# Create the actual functions that are run from outside of this file
def Missing_Atoms(molecules, experiment = None):
    """Add missing Atoms to Molecules."""
    # Validate the Molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure, making sure it is OK (it is, but whatever)
    procedure = validate_procedure("add_missing_atoms")
    # Determine who is running this experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    script = introduction("Add missing Atoms to Molecules.")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill() 
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the new structures of the Molecules
    load_structures(molecules, procedure, "all")
    # Clean up after the procedure
    clean_up(molecules, input, output, procedure)

def Relaxation(molecules, experiment = None):
    """Use CHARMM to run an energy minimization."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure
    procedure = validate_procedure("relaxation")
    # Determine who is running the experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    if gn:
        script = introduction("The Relaxation of Design Group " + str(gn))
    else:
        script = introduction("A Structure Relaxation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill() + minimize(molecules, experiment, procedure, gn)
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the new structures of the Molecules
    load_structures(molecules, procedure, "all")
    # Clean everything up
    clean_up(molecules, input, output, procedure)

def Perturbation(molecules, angles, experiment = None):
    """Use CHARMM to do a backbone perturbation."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Declare the procedure
    procedure = validate_procedure("perturbation")
    # Determine who is running the experiment
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Create the CHARMM script
    if isinstance(gn, int):
        script = introduction("The Perturbation of Design Group " + str(gn))
    else:
        script = introduction("A Backbone Perturbation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, "all")
    script += ic_fill()
    # Specify the names of the Atoms to use in the restraints based on the file
    # format
    if molecules[0].fileFormat == "PDB":
        A1 = ' C '
        A2 = ' N '
        A3 = ' CA '
        A4 = ' C '
        A5 = ' N '
    # If the file format isn't supported
    else:
        text = "The CHARMM Perturbation function does not support the "
        text += str(molecules[0].fileFormat) + " file format."
        raise CHARMM_Error(text)
    # Loop through the Molecules
    for molecule in molecules:
        # Calculate the Molecule's dihedral angles
        molecule.calculate_dihedrals()
        # Loop through the Residues by index
        for i in range(len(molecule)):
            # Get the Residue, the Residue's Molecule's name, the Residue's
            # name, and the Residue's number
            res2 = molecule[i]
            mn = res2.moleculeName
            rn = res2.name
            n = res2.number
            # If this Residue has perturbed angles
            if mn in angles and rn in angles[mn]:
                # Do the phi angle - iff there is a phi value listed for this
                # Residue, it is not the first Residue in the Molecule, the
                # listed value is a floating point number, and the Residue has a
                # phi dihedral angle
                if "phi" in angles[mn][rn] and i != 0 and res2.phi != None and \
                isinstance(angles[mn][rn]["phi"], float):
                    # Try to access the appropriate Atoms to make sure they 
                    # exist and that at least one of them can move
                    error = False
                    try:
                        # Get the previous Residue
                        res1 = molecule[i-1]
                        # Make a list of the Residues / Atoms to validate the
                        # existance of and make sure they're not all fixed in
                        # place
                        sets = [[res1, A1], [res2, A2], [res2, A3], [res2, A4]]
                        count = 0
                        for set in sets:
                            atom = set[0][set[1].strip()]
                            if set[0].freedom == "FIXED" or atom.name in \
                            set[0].fixedAtoms:
                                count += 1
                        # If all 4 Atoms are fixed in place, skip this restraint
                        if count == 4:
                            error = True
                    # If there was an error accessing any of the Atoms, skip the
                    # restraint
                    except MOLECULES.MoleculeError:
                        error = True
                    # Get the new phi angle
                    if not error:
                        phi = res2.phi + angles[mn][rn]["phi"]
                        if phi > 180.0:
                            phi -= 360.0
                        if phi < -180.0:
                            phi += 360.0
                        phi = format(phi, '.3f') + "\n"
                        # Create the restraint
                        script += "cons dihe " + str(n-1) + A1+str(n)+A2+str(n)
                        script += A3 + str(n) + A4 + "force 32800.0 min "+phi
                # Do the same for psi
                if "psi" in angles[mn][rn] and i != len(molecules) - 1 and \
                res2.psi != None and isinstance(angles[mn][rn]["psi"], float):
                    error = False
                    try:
                        res3 = molecule[i+1]
                        sets = [[res2, A2], [res2, A3], [res2, A4], [res3, A5]]
                        count = 0
                        for set in sets:
                            atom = set[0][set[1].strip()]
                            if set[0].freedom == "FIXED" or atom.name in \
                            set[0].fixedAtoms:
                                count += 1
                        if count == 4:
                            error = True
                    except MOLECULES.MoleculeError:
                        error = True
                    if not error:
                        psi = res2.psi + angles[mn][rn]["psi"]
                        if psi > 180.0:
                            psi -= 360.0
                        if psi < -180.0:
                            psi += 360.0
                        script += "cons dihe " + str(n) + A2 + str(n) + A3
                        script += str(n) + A4 + str(n+1) + A5 + " force 32800.0"
                        script += " min " + format(psi, '.3f') + "\n"
    # Now finish the script by minimizing the structure and outputting the
    # molecules
    script += minimize(molecules, experiment, procedure, gn)
    script += output_molecules(molecules, procedure, "all") + "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Load the Molecules' new structures
    load_structures(molecules, procedure, "all")
    # Clean up the files
    clean_up(molecules, input, output, procedure)

def Energy(molecules, experiment = None, which = "all"):
    """Use CHARMM to calculate the complex energy of a group of Molecules."""
    # Validate the molecules
    molecules, gn = validate_molecules(molecules)
    # Create the procedure
    procedure = validate_procedure("energy")
    # Determine who is doing the calculation
    try:
        user = experiment["User"]
    except (KeyError, TypeError, AttributeError, IPRO_Error):
        user = defaultUser
    # Determine if solvation is or should be used
    solvation = SOLVATION.get_string(experiment, procedure)
    # Generate the CHARMM script to do the energy calculation
    if isinstance(gn, int):
        script = introduction("Energy Calculation for Design Group " + str(gn))
    else:
        script = introduction("A CHARMM energy calculation")
    script += load_input_files(experiment)
    script += load_molecules(molecules, procedure, user, which)
    script += ic_fill() + solvation
    script += energy_line(experiment, procedure, '', '', '', solvation)
    # Declare the name of the output file
    energyFile = "energy_values.txt"
    script += "\nener\n\nset tot ?ener\n\nopen write card unit 10 name "
    script += energyFile + "\n\nwrite title unit 10\n* @tot\n*\nclose unit 10\n"
    script += "stop\n"
    # Run CHARMM
    input, output = execute_CHARMM_script(script, procedure, gn)
    # Read in the calculated energy
    try:
        file = open(energyFile, "r")
        energy = float(file.readline())
        file.close()
    except IOError:
        text = "A CHARMM energy calculation failed to calculate an energy. "
        text += "Please review " + output + " to determine why."
        raise CHARMM_Error(text)
    # Clean up the folder
    #clean_up(molecules, input, output, procedure, energyFile)
    # Return the calculated energy
    return energy
