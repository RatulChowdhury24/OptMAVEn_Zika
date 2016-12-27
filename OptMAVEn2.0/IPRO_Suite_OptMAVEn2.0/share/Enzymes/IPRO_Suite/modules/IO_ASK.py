#!/usr/bin/env python

# The name of the file
__name__ = "IPRO Suite Functions to Ask for Attributes"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This file contains functions that ask a user for specific attribute
information."""

# Standard PYTHON modules
import os
import sys
# Other I/O modules
import IO_FUNCTIONS as FUNCTIONS
import IO_LOADING as LOADING
import IO_CHECK as CHECK
import IO_GET as GET
# Include all contents from the STANDARDS module
from STANDARDS import *

# Have functions to ask for basic pieces of information
def for_user(experiment):
    """Ask who the user is"""
    # Get the information from the system
    try:
        user = os.getlogin()
    except OSError:
        user = defaultUser
    # Confirm they are this person
    question = "It appears you are " + user + ". Is this correct?"
    default = FUNCTIONS.answer_question(question, "bool")
    # And use the FUNCTIONS.get_value function to store the value
    question = "What is your name?"
    FUNCTIONS.get_value("User", "string", experiment, [], default, user, \
                        question)

def for_type(experiment):
    """Find out what type of IPRO Suite Experiment this is"""
    question = "What type of IPRO Suite Experiment is this?"
    FUNCTIONS.get_value("Type", "string", experiment, supportedPrograms, False,\
                        None, question)

def for_name(experiment):
    """Ask what the Experiment should be named"""
    question = "What would you like to name this " + experiment["Type"] 
    question += " Experiment?"
    FUNCTIONS.get_value("Name", "string", experiment, [], False, None, question)

def for_format(experiment):
    """Get the Experiment's file format"""
    if len(supportedFormats) == 1:
        default = True
    else:
        default = False
    question = "What type of file format would you like to use in this "
    question += "experiment?"
    FUNCTIONS.get_value("File Format", "string", experiment, supportedFormats, \
                        default, defaultFormat, question)

def for_field(experiment):
    """Get the force field"""
    if len(supportedFields) == 1:
        default = True
    else:
        default = False
    question = "What force field would you like to use in this experiment?"
    FUNCTIONS.get_value("Force Field", "string", experiment, supportedFields, \
                        default, defaultField, question)

def basic_info(experiment):
    """Ask the user for the basic information for an IPRO Experiment"""
    # Clear the screen
    os.system("clear")
    # Tell the user what is happening
    message = """
Welcome to the Iterative Protein Redesign & Optimization (IPRO) Suite of
Programs from the Costas Maranas Laboratory in the Chemical Engineering
Department of the Pennsylvania State University.

Please make sure you have read the provided documentation. That is where the
methodologies, capabilities, and terminology of the IPRO Suite are defined and
explained. It is assumed you have a working knowledge of this information.

To begin, please provide some basic information about your experiment."""
    print screen_formatting(message[1:])
    # Make sure the current folder won't cause any later problems
    folder = os.getcwd() + "/"
    if CHECK.has_whitespace(folder):
        text = "The current folder is:\n" + folder + "\nThis contains "
        text += "whitespace, which will cause the IPRO Suite to fail. Sorry, "
        text += "but you must start over from a different folder."
        raise FUNCTIONS.IPRO_IOError(text)
    # Get the basic pieces of information
    for_user(experiment)
    for_type(experiment)
    for_name(experiment)
    # Store the Experiment's folder
    if "experiments" in os.listdir("./"):
        folder += "experiments/"
    experiment["Folder"] = folder + experiment["Name"] + "/"
    # Force field and file format
    for_field(experiment)
    for_format(experiment)

def for_CHARMM_topology_files(experiment, default):
    """Ask the user for the CHARMM Topology files"""
    # Just use the get_list function
    FUNCTIONS.get_list("CHARMM Topology Files", experiment, [], default, \
                       defaultCHARMMTopologies)

def for_CHARMM_parameter_files(experiment, default):
    """Ask for parameter files"""
    FUNCTIONS.get_list("CHARMM Parameter Files", experiment, [], default, \
                       defaultCHARMMParameters)

def for_CHARMM_energies(experiment, default):
    """Ask for the energy terms to calculate in CHARMM"""
    FUNCTIONS.get_list("CHARMM Energy Terms", experiment, [], default, \
                       defaultCHARMMEnergies)

def for_CHARMM_iterations(experiment, default):
    """Ask for the maximum number of iterations to run an energy minimization"""
    question = "What is the maximum number of iterations a CHARMM energy "
    question += "minimization may run for?"
    question += " The default is " + str(defaultCHARMMIterations) + "."
    FUNCTIONS.get_value("CHARMM Iterations", "integer", experiment, [], \
                        default, defaultCHARMMIterations, question)

def for_CHARMM_softening(experiment, default):
    """Ask if a Lennard-Jones softening term should be included in C++
    calculations"""
    accept = False
    while not accept:
        if default:
            useSoftening = defaultUseSoftening
        else:    
            question = "Would you like to use a Lennard-Jones softening "
            question += "term during C++ energy calculations?"
            question += " The default is " 
            if defaultUseSoftening:
                question += "'yes'."
            else:
                question += "'no'."
            useSoftening = FUNCTIONS.answer_question(question, "bool")
        accept = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")          
    if useSoftening:
        question = "What value should be used for the softening term?"
        question += " The default is " + str(defaultLjPhi) + "."
        FUNCTIONS.get_value("Lj Phi", "float", experiment, [], default, \
                            defaultLjPhi, question)
    else:
        experiment["Lj Phi"] = 1.0   

def CHARMM_info(experiment):
    """Ask the user how to use CHARMM"""
    # Clear the screen
    os.system("clear")
    # Create a message for the user
    message = """
You now need to specify information about how to use the CHARMM force field.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # find out if default settings should be used
    question = "Would you like to use the default CHARMM settings? Please note "
    question += "that they are only valid for amino and nucleic acids."
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the other values
    for_CHARMM_topology_files(experiment, default)
    for_CHARMM_parameter_files(experiment, default)
    for_CHARMM_energies(experiment, default)
    for_CHARMM_iterations(experiment, default)
    for_CHARMM_softening(experiment, default)

def for_docking_frequency(experiment, default):
    """Ask how often to run docking."""
    question = "Docking is not run every iteration of IPRO, but rather every N "
    question += "iterations. What should N be in your experiment?"
    question += " The default is " + str(defaultDockingFrequency) + "."
    FUNCTIONS.get_value("Docking Frequency", "integer", experiment, [], \
                        default, defaultDockingFrequency, question)

def for_docking_iterations(experiment, default):
    """Ask how long to run docking for"""
    question = "How many iterations of random structure movements should be "
    question += "carried out every time docking is run?"
    question += " The default is " + str(defaultDockingIterations) + "."
    FUNCTIONS.get_value("Docking Iterations", "integer", experiment, [], \
                        default, defaultDockingIterations, question)

def for_docking_sd_move(experiment, default):
    """Ask for the standard deviations on random movements"""
    question = "What should be the standard deviation of the random cartesian "
    question += "structure perturbations? The default is "
    question += format(defaultSDMove, '.3f') + " Angstroms."
    FUNCTIONS.get_value("Docking SD Movements", "float", experiment, [], \
                        default, defaultSDMove, question)

def for_docking_sd_rotate(experiment, default):
    """Ask for the standard deviation on rotations"""
    question = "What should be the standard deviation of the random structure "
    question += "rotations? The default is " + format(defaultSDRotate, '.3f')
    question += " degrees."
    FUNCTIONS.get_value("Docking SD Rotations", "float", experiment, [], \
                        default, defaultSDRotate, question)

def for_docking_start_temp(experiment, default):
    """Ask for the temperature at the start of docking."""
    question = "What should be the simulated annealing temperature at the start"
    question += " of docking? The default is " + format(defaultDockTempStart, \
                '.3f') + " Kelvin."
    FUNCTIONS.get_value("Docking Start Temp", "float", experiment, [], default,\
                        defaultDockTempStart, question)

def for_docking_end_temp(experiment, default):
    """Ask for the temperature at the end of docking."""
    question = "What should be the simulated annealing temperature at the end "
    question += "of docking? The default is " + format(defaultDockTempEnd,'.3f')
    question += " Kelvin."
    FUNCTIONS.get_value("Docking End Temp", "float", experiment, [], default, \
                        defaultDockTempEnd, question)

def docking_info(experiment):
    """Ask the user how to run docking"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You must now specify how to run docking.

Please refer to the provided documentation for additional information.

As a reference, simulated annealing temperatures are calculated as follows, with
units of Kelvin (energies have units of kcal / mol).

T = -(delta energy) / (Gas Constant * LN(fraction of retained positions))
3640 = -10 / (0.001986 * LN(1/4))"""
    print screen_formatting(message[1:])
    # Find out about using defaults
    question = "Would you like to use the default docking settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the values
    for_docking_frequency(experiment, default)
    for_docking_iterations(experiment, default)
    for_docking_sd_move(experiment, default)
    for_docking_sd_rotate(experiment, default)
    for_docking_start_temp(experiment, default)
    for_docking_end_temp(experiment, default)

def for_IPRO_iterations(experiment, default):
    """Ask for how many iterations to run IPRO"""
    question = "How many iterations of IPRO would you like to run? The default "
    question += "is " + str(defaultIterations) + "."
    FUNCTIONS.get_value("IPRO Iterations", "integer", experiment, [], default, \
                        defaultIterations, question)

def for_IPRO_temperature(experiment, default):
    """Ask the user for the simulated annealing temperature."""
    question = "What simulated annealing temperature would you like to use? The"
    question += " default is " + format(defaultAnnealingTemp, '.3f') + " Kelvin"
    FUNCTIONS.get_value("IPRO Annealing Temperature", "float", experiment, [], \
                        default, defaultAnnealingTemp, question)

def for_annealing_sharing(experiment, default):
    """Ask if the results of simulated annealing should be shared."""
    question = "Should results kept by simulated annealing be shared between "
    question += "different processors? The default is "
    if defaultShareAnnealing:
        question += "yes."
    else:
        question += "no."
    FUNCTIONS.get_value("Annealing Sharing", "bool", experiment, [], default, \
                        defaultShareAnnealing, question)

def for_energy_type(experiment, default):
    """Ask what type of energy calculations to do."""
    question = "Would you like to use 'Interaction' or 'Binding' energies to "
    question += "evaluate the results of the experiment? The default is "
    question += defaultEnergyCalc + ". Binding energies can only be used when "
    question += "structure refinements are used."
    FUNCTIONS.get_value("Energy Calculation", "string", experiment, ["Binding",\
                        "Interaction"], default, defaultEnergyCalc, question)

def IPRO_info(experiment):
    """Ask the user questions about how to run IPRO."""
    # Clear the screen and tell the user what's what
    os.system("clear")
    message = """
You now need to provide information about how to run iterations of IPRO during
your experiment.

Please refer to the provided documentation for additional information.

As a reminder, simulated annealing temperatures are calculated as follows and
have units of Kelvin (energies have units of kcal / mol).

T = -(delta energy) / (Gas Constant * LN(fraction of retained designs))
3640 = -10 / (0.001986 * LN(1/4))"""
    print screen_formatting(message[1:])
    # Find out if defaults should be used
    question = "Would you like to use the default IPRO settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the values
    for_IPRO_iterations(experiment, default)
    for_IPRO_temperature(experiment, default)
    for_annealing_sharing(experiment, default)
    for_energy_type(experiment, default)

def for_do_refinement(experiment, default):
    """Ask whether or not structure refinements should happen."""
    # If this is a mutator experiment, the answer is just yes
    if experiment["Type"] == "Mutator":
        experiment["Do Refinement"] = True
    else:
        question = "Would you like to do extensive structure refinements to "
        question += "improve the quality of results? The default is "
        if defaultDoRefinement:
            question += "yes."
        else:
            question += "no."
        FUNCTIONS.get_value("Do Refinement", "bool", experiment, [], default, \
                            defaultDoRefinement, question)

def for_refinement_iterations(experiment, default):
    """Ask how many iterations each structure should be refined for."""
    question = "How many iterations of IPRO should be carried out for each "
    question += "structure during refinements? The default is "
    question += str(defaultRefinementIterations) + "."
    FUNCTIONS.get_value("Refinement Iterations", "integer", experiment, [], \
                        default, defaultRefinementIterations, question)

def for_ensemble_number(experiment, default):
    """Ask how many structures to make during refinements."""
    question = "How many structures should be in each ensemble during "
    question += "structure refinements? The default is "
    question += str(defaultEnsembleNumber) + "."
    FUNCTIONS.get_value("Ensemble Number", "integer", experiment, [], default, \
                        defaultEnsembleNumber, question)

def refinement_info(experiment):
    """Ask the user how and if to do structure refinements"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
IPRO offers the option of doing extensive, ensemble-based structure refinements
to improve the quality of predictions. You must now provide the information on
how to do that.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Find out if defaults should be used
    question = "Would you like to use the default refinement settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the values
    for_do_refinement(experiment, default)
    for_refinement_iterations(experiment, default)
    for_ensemble_number(experiment, default)

def mutator_info(experiment):
    """Ask the user how to do structure refinements during Mutators"""
    # Clear the screen and tell the user what's what
    os.system("clear")
    message = """
During Mutator experiments, extensive structure refinements are carried out that
provide more accurate assesments of the generated mutants.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Find out if default settings should be used
    question= "Would you like to use the default structure refinement settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    for_do_refinement(experiment, default)
    for_refinement_iterations(experiment, default)
    for_ensemble_number(experiment, default)

def for_rotamer_library(experiment, default):
    """Ask where the rotamer libary is located."""
    question = "Where is the rotamer libary located? The default is:\n"
    question += RotamerLibraryPath
    FUNCTIONS.get_value("Rotamer Library", "string", experiment, [], default, \
                        RotamerLibraryPath, question)

def for_rotamer_window(experiment, default):
    """Ask what the rotamer window should be"""
    question = "In kcal / mol, what should the 'rotamer window' value be? The "
    question += "default is " + format(defaultRotamerWindow, '.3f') + "."
    FUNCTIONS.get_value("Rotamer Window", "float", experiment, [], default, \
                        defaultRotamerWindow, question)

def for_max_rotamers(experiment, default):
    """Ask for the maximum number of rotamers"""
    question = "For memory reasons, what is the maximum number of rotamers that"
    question += " may be considered by the rotamer selection MILP? The default "
    question += "is " + str(defaultMaxRotamers) + "."
    FUNCTIONS.get_value("Max Rotamer Number", "integer", experiment, [], \
                        default, defaultMaxRotamers)

def for_packing_method(experiment, default):
    """Ask how to repack rotamers - by sequence or distance"""
    question = "Would you like to only repack the sidechains of Residues that "
    question += "are close to perturbed regions in 'Sequence' or also by "
    question += "'Distance'? The default is " + defaultPackingMethod + "."
    FUNCTIONS.get_value("Packing Method", "string", experiment, ['Sequence', \
                        'Distance'], default, defaultPackingMethod, question)

def for_packing_selection(experiment, default):
    """Ask for which Residues to repack around"""
    question = "Would you like to repack the sidechains of Residues that are "
    question += "close only to the 'RESTRAINED' Residues during a peturbation, "
    question += "or also those that are close to the 'FREE' Residues? The "
    question += "default is " + defaultPackingSelection + "."
    if experiment["Packing Method"] == "Sequence":
        question += " This information is still needed even when the Packing "
        question += "Method is 'Sequence'."
    FUNCTIONS.get_value("Packing Selection", "string", experiment, ["FREE", \
                    "RESTRAINED"], default, defaultPackingSelection, question)

def for_packing_cutoff(experiment, default):
    """Ask for the distance cutoff for rotamer repacking"""
    question = "How close must a heavy Atom in a Residue be to a heavy Atom in "
    question += "the Packing Selection for the Residue to be repacked? The "
    question += "default is " + format(defaultPackingCutoff, '.3f')
    question += " Angstroms."
    if experiment["Packing Method"] == "Sequence":
        question += " This information is still needed even when the Packing "
        question += "Method is 'Sequence'."
    FUNCTIONS.get_value("Packing Cutoff", "float", experiment, [], default, \
                        defaultPackingCutoff, question)

def rotamer_info(experiment):
    """Ask the user how to use Rotamers"""
    # Clear the screen and give a message to the user
    os.system("clear")
    message = """
You must now provide information about how to use Rotamers.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Ask if they want to use defaults
    question = "Would you like to use the default rotamer settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the values
    for_rotamer_library(experiment, default)
    for_rotamer_window(experiment, default)
    for_max_rotamers(experiment, default)
    for_packing_method(experiment, default)
    for_packing_selection(experiment, default)
    for_packing_cutoff(experiment, default)

def for_use_solvation(experiment, default):
    """Ask the user if they want to use implicit solvation."""
    question = "Would you like to use implicit solvation during your "
    question += "experiment? The default is "
    if defaultUseSolvation:
        question += "yes."
    else:
        question += "no."
    FUNCTIONS.get_value("Use Solvation", "bool", experiment, [], default, \
                        defaultUseSolvation, question)
    # Store the settings for the specific IPRO functions 
    for procedure in ["Relaxation", "Perturbation", "Energy"]:
        experiment[procedure + " Solvation"] = experiment["Use Solvation"]

def for_solvation_type(experiment, default):
    """Ask the user what type of implicit solvation to use"""
    question = "What type of implicit solvation would you like to use? The "
    question += "default is " + defaultSolvationType + "."
    if len(supportedSolvations) == 1 or default:
        default2 = True
    FUNCTIONS.get_value("Solvation Type", "string", experiment, \
                supportedSolvations, default, defaultSolvationType, question)

def for_LK_files(experiment, default):
    """Ask for Lazaridis-Karplus solvation files"""
    FUNCTIONS.get_list("LK Solvation Files", experiment, [], default, \
                       defaultLKSolvationFiles)

def solvation_info(experiment):
    """Ask the user for information about using implicit solvation"""
    # Clear the screen and say what is happening
    os.system("clear")
    message = """
You must now provide information about how to use implicit solvation in your
experiment. 

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Ask about using the default settings
    question = "Would you like to use the default implicit solvation settings? "
    question += "Please note that they are only valid for amino acid systems."
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the solvation settings
    for_use_solvation(experiment, default)
    if experiment["Use Solvation"]:
        for_solvation_type(experiment, default)
        if experiment["Solvation Type"] == "Lazaridis-Karplus":
            for_LK_files(experiment, default)

# Now things get more complicated, and the asking for functions become more
# involved.
def Standard_Molecules(experiment):
    """Ask the user to specify Molecules"""
    # Clear the screen and say what is happening
    os.system("clear")
    message = """
You now need to select all of the Molecules you will be using in your
Experiment. First, identify all of the Design Molecules. Once that is done you
will be able to select the Target Molecules.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Store files of molecules here
    files = {}
    # And the list of Molecules used in the experiment here
    molecules = []
    # Get the Design Molecules
    GET.structures(molecules, files, "Design Molecule", "Design Molecules",\
                   experiment["Force Field"], experiment["File Format"])
    # Get the Target Molecules
    GET.structures(molecules, files, "Target Molecule", "Target Molecules", \
                   experiment["Force Field"], experiment["File Format"])
    # Store the Molecules in the Experiment
    experiment["Molecules"] = molecules
    # And search them to find out if any are Dimers
    GET.Dimers(experiment)

def DesignPositions(experiment):
    """Ask the user to specify Design Positions in a Molecule"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You must now identify the Design Positions for your experiment. These are amino
acid Residues in Design Molecules that are allowed to mutate. In other words,
they are the Residues that you are designing to improve the properties of your
system.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Make a list of the Design Molecules
    molecules = {}
    order = []
    for data in experiment["Molecules"]:
        if data[2].design:
            order.append(data[2].name)
            molecules[data[2].name] = data[2]
    if len(order) == 0:
        text = "Somehow there are no Design Molecules"
        raise FUNCTIONS.IPRO_IOError(text)
    # Store the Design Positions in this dictionary
    positions = {}
    # Use a while loop
    accept = False
    pick = True
    while not accept:
        # If Design Positions should be found
        if pick:
            # Get the Molecule
            if len(order) == 1:
                text = "\nMolecule " + order[0] + " was automatically selected "
                text += "because it is the only Design Molecule."
                print screen_formatting(text)
                mn = order[0]
            else:
                # Summarize the possible Design Molecules
                summary = ''
                for mn in order:
                    summary += "\nDesign Molecule " + mn + " contains "
                    summary += str(len(molecules[mn])) + " Residues"
                print screen_formatting(summary)
                # Ask which Molecule, giving an option of selecting None
                question = "In which Molecule would you like to specify Design "
                question += "Positions? You may answer 'none' if you would not "
                question += "like to pick any."
                answers = list(order)
                answers.extend(["NONE", 'None', 'none'])
                mn = FUNCTIONS.answer_question(question, "string", [], answers)
                if mn.upper() == "NONE":
                    pick = False
        # if a Molecule has been selected
        if pick:
            GET.DesignPositions(positions, molecules[mn], experiment)
            # if there are multiple Molecules, give the option of picking more
            # choices
            if len(order) > 1:
                question = "Would you like to specify Design Positions in "
                question += "another Molecule?"
                pick = FUNCTIONS.answer_question(question, "bool")
            else:
                pick = False
        # If there are no selected Design Positions, give the user an error
        if positions == {}:
            text = "\nYou must select at least one Design Position"
            print screen_formatting(text)
            pick = True
            continue
        # If a selection either shouldn't or "can't" be made, give the options
        # to be done
        if not pick:
            # Summarize the selected Design Positions
            summary = ''
            for mn in order:
                if mn in positions:
                    summary += "\nDesign Molecule " + mn + ": " + \
                        str(len(positions[mn])) + " Design Positions selected"
            print screen_formatting(summary)
            # Ask if the user is happy with this
            question = "Are you certain you have specified the correct Design "
            question += "Positions for your experiment?"
            right = FUNCTIONS.answer_question(question, "bool")
            # If the user is happy, be done
            if right:
                accept = True
            # if not, let them modify their choices
            else:
                text = "\nThen modify the selections you've made until you are."
                print screen_formatting(text)
                pick = True
    # Store the selected Design Positions in the Experiment
    experiment["Design Positions"] = positions

def PermittedKinds(experiment):
    """Ask the user to identify how Residues may mutate"""
    # clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You now have the OPTION of identifing how the Design Positions may mutate. By
default, positions that do not have explicitly listed kinds of amino acids they
may mutate to are permitted to mutate to any kind of amino acid.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Store the permitted kinds of amino acids in this dictionary
    permitted = {}
    # Use a while loop
    accept = False
    while not accept:
        question = "Would you like to select permitted kinds of amino acids for"
        question += " a Design Position?"
        do = FUNCTIONS.answer_question(question, "bool")
        if do:
            GET.PermittedKinds(permitted, experiment)
            continue
        # Confirm the user is happy with their selections
        question = "Are you certain you are happy with the permitted kinds of "
        question += 'amino acids that you have specified?'
        right = FUNCTIONS.answer_question(question, "bool")
        if right:
            accept = True
    # Store the permitted kinds of amino acids
    experiment["Permitted Kinds"] = permitted

def DesignGroups(experiment):
    """Identify how binding is being designed"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You now need to specify Design Groups for your experiment. A Design Group is a
group of Molecules that are all present in your system simultaneously. All
Design Molecules are automatically included, at least one Target Molecule must
be included, and each Design Group must be unique. You also specify how binding
is modified for the Design Group, with the options being:
1) improve - make sure binding continually gets better every time a design is
retained.
2) maintain - make sure binding does not get worse than it was initially.
3) reduce - make sure binding does not improve over its initial value.
4) eliminate - make sure binding continually gets worse every time a design is
retained.

Please note that the order you specify Design Groups in is VERY important. Make
sure the most important group is listed FIRST. The order of the groups after
that does not matter, but the first group has a MAJOR impact on the results of
every iteration of IPRO.

This is a very important, and potentially confusing, topic so please refer to
the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Store the Design Groups in this list
    groups = []
    # Make a list of the Experiment's Target Molecules
    molecules = []
    for data in experiment["Molecules"]:
        if not data[2].design:
            molecules.append(data[2])
    if len(molecules) == 0:
        text = "Somehow there are no Target Molecules"
        raise FUNCTIONS.IPRO_IOError(text)
    # Get the list of all possible Design Groups
    names = []
    for molecule in molecules:
        names.append(molecule.name)
    all = GET.all_combinations(names)
    # Use a while loop to get the Design Groups
    accept = False
    pick = True
    while not accept:
        # If a Design Group should be collected, do so
        if pick:
            GET.DesignGroup(groups, molecules)
        # Determine if every Target Molecule has been used in at least one
        # Design Group
        unused = []
        for name in names:
            have = False
            for group in groups:
                if name in group:
                    have = True
                    break
            if not have:
                unused.append(name)
        # If there are Molecules that haven't been used yet
        if len(unused) > 0:
            text = "\nThe following Target Molecules must be used in a "
            text += "Design Group:" + list_items(unused)
            print screen_formatting(text)
            pick = True
            continue
        # If it is possible to pick more Design Groups, ask for that option
        if len(groups) < len(all):
            question = "Would you like to specify another Design Group?"
            pick = FUNCTIONS.answer_question(question, "bool")
            if pick:
                continue
        else:
            pick = False
        # At this point, all Molecules have been used in a Design Group and
        # there are either no more or the user doesn't want to make a choice
        if not pick:
            # Summarize the Design Groups
            summary = ''
            for i, group in enumerate(groups):
                summary += "\n" + str(i+1) + ") " + group[0].capitalize()
                summary += " binding to" + list_items(group[1:])
            print screen_formatting(summary)
            # Make sure the user is happy
            question = "Are you certain you are happy with the Design Groups "
            question += "you have specified?"
            right = FUNCTIONS.answer_question(question, "bool")
            if right:
                accept = True
                continue
            # If they aren't happy, give them options
            phrases = ["'start over' at picking Design Groups"]
            answers = ["START OVER", "Start Over", 'Start over', 'start over']
            if len(groups) < len(all):
                phrases.append("select 'more' Design Groups")
                answers.extend(["MORE", "More", 'more'])
            if len(groups) > 1:
                phrases.append("'remove' a specific Design Group")
                answers.extend(["REMOVE", "Remove", "remove"])
            # It is guaranteed that there is at least one phrase
            if len(phrases) == 1:
                what = answers[0].upper()
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1]
                else:
                    question += ", " + phrases[1] + ', or ' + phrases[2]
                question += "?"
                what = FUNCTIONS.answer_question(question, "string", [], \
                                                 answers)
            # Respond appropriately
            if what.upper() == "MORE":
                pick = True
                continue
            elif what.upper() == "START OVER":
                pick = True
                groups = []
                continue
            # Otherwise, ask which one to delete
            accept2 = False
            while not accept2:
                answers = range(1, len(groups) + 1)
                question = "Which Design Group should be removed from the list?"
                i = FUNCTIONS.answer_question(question, "integer", [], answers)
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept2 = True
            # Delete that Design Group
            del groups[i-1]
    # Store the Design Groups in the Experiment
    experiment["Design Groups"] = groups

def OptzymeGroups(experiment):
    """Identify how binding is being optimized during an OptZyme experiment."""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You must now specify the Optzyme Groups used in your experiment. An Optzyme
Group is a collection of molecules that is involved in optimization of a given
catalytic property. Based on the catalytic property, the number of Binding
Assemblies varies. A Binding Assembly is a group of Molecules that are
simultaneously present in your system. In each Optzyme Group, all Design
Molecules are automatically included, at least one Target Molecule must be
included, and each Optzyme Group must be unique. The catalytic properties that
may be modified within OptZyme are:
1) km - modify binding of the substrate to form the Michaelis complex. One
        Binding Assembly with ground state molecule(s) must be provided. 
2) kcat/km - modify binding of the transition state analogue (TSA). One Binding
             Assembly with the TSA molecule(s) must be provided.
3) kcat - modify a weighted difference of binding affinity between the TSA and
          the substrate. Two Binding Assemblies (one with the substrate(s) and
          one with the TSA molecule(s)) must be provided. Additionally, a
          weighting factor (RT(TSA)/RT(S)) to weight the difference of the two
          binding energies must be included. Thus, the objective function for a
          kcat calculation is: obj = BE(TSA) - (RT(TSA)/RT(S))*BE(S)

You must also specify how binding is modified for the Optzyme Group, with the
options being:
1) improve - make sure binding continually gets better every time a design is
retained.
2) maintain - make sure binding does not get worse than it was initially.
3) reduce - make sure binding does not improve over its initial value.
4) eliminate - make sure binding continually gets worse every time a design is
retained.

Please note that the order you specify Optzyme Groups is VERY important. Ensure
that the most important group is listed FIRST. The order of the other groups is
irrelevant, but the first group has a MAJOR impact on the results of every
iteration of OptZyme.

This is a very important, and potentially confusing, topic so please refer to
the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Store the Optzyme Groups in this list
    groups = []
    # Make a list of the Experiment's Target Molecules
    molecules = []
    for data in experiment["Molecules"]:
        if not data[2].design:
            molecules.append(data[2])
    # If no target molecules exist, raise an error
    if len(molecules) == 0:
        text = "Somehow there are no Target Molecules"
        raise FUNCTIONS.IPRO_IOError(text)
    # Get the list of all possible Optzyme Groups
    names = []
    for molecule in molecules:
        names.append(molecule.name)
    # Multiply the total number by 2 since the set of molecules can be used in
    # multiple states    
    all = GET.all_combinations(names) * 2
    # Use a while loop to acquire the Optzyme Groups
    accept = False
    pick = True
    while not accept:
        # If an Optzyme Group should be collected, do so
        if pick:
            GET.OptzymeGroup(groups, molecules)
        # Determine if every Target Molecule has been used in at least one
        # Optzyme Group
        unused = []
        for name in names:
            have = False
            for group in groups:
                for i in range(len(group[3])):
                    if name in group[3][i][1:]:
                        have = True
                        break
            if not have:
                unused.append(name)
        # If there are Target Molecules that have not been used yet
        if len(unused) > 0:
            text = "\nThe following Target Molecules must be used in an "
            text += "Optzyme Group:" + list_items(unused)
            print screen_formatting(text)
            pick = True
            continue
        # If it is possible to pick more Optzyme Groups, ask for that option
        if len(groups) < len(all):
            question = "Would you like to specify another Optzyme Group?"
            pick = FUNCTIONS.answer_question(question, "bool")
            if pick:
                continue
        else:
            pick = False
        # At this point, all Molecules have been used in an Optzyme Group and 
        # the user does not want to make another Optzyme Group
        if not pick:
            # Summarize the Optzyme Groups
            summary = ''
            for i, group in enumerate(groups):
                summary += "\n" + str(i+1) + ") " + group[1].capitalize()
                summary += " " + group[0].capitalize() 
                if group[0] == "kcat":
                    summary += " with a weight of " + str(group[2])
                for j in [0, 1]:
                    if j == 1:
                        if len(group[3]) == 1:
                            break
                        elif len(group[3]) == 2:
                            summary += "and"
                    summary += " by binding to" + list_items(group[3][j][1:])
                    summary += " in the "
                    if group[3][j][0] == "improve":
                        summary += "TSA state "
                    else:
                        summary += "ground state "
            print screen_formatting(summary)
            # Ensure that the user is satisfied
            question = "Are you certain you are happy with the Optzyme Groups "
            question += "that you have specified?"
            right = FUNCTIONS.answer_question(question, "bool")
            if right:
                accept = True
                continue
            # If the user is unhappy, give them additional options
            phrases = ["'start over' at selecting Optzyme Groups"]
            answers = ["START OVER", "Start Over", "Start over", "start over"]
            if len(groups) < len(all):
                phrases.append("select 'more' Optzyme Groups")
                answers.extend(["MORE", "More", "more"])
            if len(groups) > 1:
                phrases.append("'remove' a specific Optzyme Group")
                answers.extend(["REMOVE", "Remove", "remove"])
            # It is guaranteed that there is at least one phrase stored
            if len(phrases) == 1:
                what = answers[0].upper()
            # If there are more than 1, list all of the options
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1]
                else:
                    question += ", " + phrases[1] + ", or " + phrases[2]
                question += "?"
                what = FUNCTIONS.answer_question(question, "string", [], \
                                                 answers)
            # Perform the desired function
            if what.upper() == "MORE":
                pick = True
            elif what.upper() == "START OVER":
                pick = True
                groups = []
            # Otherwise, ask which Optzyme Group should be deleted
            else:
                accept2 = False
                while not accept2:
                    answers = range(1, len(groups) + 1)
                    question = "Which Optzyme Group should be removed from the"
                    question += " list?"
                    i = FUNCTIONS.answer_question(question, "integer", [], \
                                                  answers)
                    # Confirm that this is the correct Optzyme Group
                    right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                    if right:
                        accept2 = True
                # Delete that Optzyme Group
                del groups[i-1]
    # Store the Optzyme Groups in the Experiment
    experiment["Optzyme Groups"] = groups

def Mutants(experiment):
    """Identify the mutants during a mutator experiment"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You must now identify the mutants that you would like to examine in this
experiment. You may either do this manually, one mutation at a time, or load
them from a prepared file. If you do load them from a file, the formatting must
match the standard output for mutations for Mutator Experiments:

Mutation:     Mutant #: Mutate Residue RN in Molecule MN to AA

The line must start exactly with the header 'Mutation:'. The number of spaces
after that is irrelevant. # is the number of the mutant, RN is the name of the
Residue being mutated, MN is the name of the Molecule Residue RN is in, and AA
is the kind of amino acid it is mutating to (3 letter, capitalized code).

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Ask how the user would like to identify the mutations
    question = "Would you like to input the mutations by 'hand' or 'load' them "
    question += "from a file?"
    answers = ['HAND', 'Hand', 'hand', 'LOAD', 'Load', 'load']
    what = FUNCTIONS.answer_question(question, "string", [], answers)
    # Respond appropriately
    if what.upper() == "LOAD":
        # use a while loop
        accept = False
        while not accept:
            question = "What is the name of the file that contains the "
            question += "mutations? It must be located in either the current "
            question += 'folder or an input_files folder.'
            fileName = FUNCTIONS.answer_question(question, "string")
            # Get the file
            try:
                f = CHECK.for_file(fileName, "./", True)
            except FUNCTIONS.IPRO_IOError as error:
                print str(error)
                continue
            lines = f.readlines()
            f.close()
            # Store the information in this dictionary
            data = {}
            for line in lines:
                try:
                    i = line.index(":")
                except ValueError:
                    continue
                attribute = line[:i].strip()
                if attribute not in data:
                    data[attribute] = []
                data[attribute].append(line[i+1:].strip())
            # Try to load that data into the Experiment
            errors = LOADING.Mutants(experiment, data)
            if len(errors) > 0:
                print screen_formatting(errors)
            else:
                accept = True
    # If the user wants to specify them by hand
    else:
        # Store the mutants here
        mutants = []
        # Use a while loop
        accept = False
        pick = True
        while not accept:
            if pick:
                GET.Mutant(mutants, experiment)
                # Ask if they want another
                question = "Would you like to pick another Mutant?"
                pick = FUNCTIONS.answer_question(question, "bool")
                if pick:
                    continue
            # If no mutant has been selected, try again
            if len(mutants) == 0:
                text = "\nYou must specify at least one mutant"
                print screen_formatting(text)
                continue
            # if the user is done, summarize their choices
            if not pick:
                summary = ''
                for i, mutant in enumerate(mutants):
                    N = len(mutant)
                    summary += "\n" + str(i+1) + ") " + str(N) + " mutation"
                    if N > 1:
                        summary += "s"
                print screen_formatting(summary)
                # Find out if the user is happy
                question = "Are you certain you have specified the correct "
                question += "mutants?"
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    continue
                # Otherwise, give options
                phrases = ["'start over' at picking mutants", \
                           "choose 'more' mutants"]
                answers=["START OVER", 'Start Over', 'Start over','start over',\
                         "MORE", "More", 'more']
                if len(mutants) > 1:
                    phrases.append("'remove' a specific mutant")
                    answers.extend(['REMOVE', 'Remove', 'remove'])
                # Assemble the question and ask
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1]
                else:
                    question += ", " + phrases[1] + " or " + phrases[2]
                question += "?"
                what = FUNCTIONS.answer_question(question, "string", [],answers)
                # Respond appropriately
                if what.upper() == "MORE":
                    pick = True
                    continue
                elif what.upper() == "START OVER":
                    pick = True
                    mutants = []
                    continue
                accept2 = False
                while not accept2:
                    question = "Which Mutant should be removed from the list?"
                    answers = range(1, len(mutants) + 1)
                    i = FUNCTIONS.answer_question(question, "integer", [], \
                                                  answers)
                    right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                    if right:
                        accept2 = True
                # Delete that mutant
                del mutants[i-1]
        # Store the mutants in the Experiment
        experiment["Mutants"] = mutants
        # Store every Residue that gets mutated as a Design Position in the
        # experiment
        positions = {}
        for mutant in mutants:
            for mutation in mutant:
                mn = mutation[0]
                rn = mutation[1]
                if mn not in positions:
                    positions[mn] = []
                if rn not in positions[mn]:
                    positions[mn].append(rn)
        experiment["Design Positions"] = positions

def fixedAtoms(experiment):
    """Ask the user about Atoms that may never move"""
    # Store the information here
    fixed = {}
    # use a while loop
    accept = False
    while not accept:
        question = "Would you like to fix one or more Atoms in place?"
        do = FUNCTIONS.answer_question(question, "bool")
        if do:
            GET.fixedAtoms(fixed, experiment)
            continue
        # Make certain the user is happy
        question = "Are you certain you are happy with the Atoms you have "
        question += "permanently fixed in place?"
        accept = FUNCTIONS.answer_question(question, "bool")
    return fixed

def position_restraints(experiment):
    """Ask the user for position restraints"""
    # Store the restraints here
    restraints = []
    # Use a while loop
    accept = False
    while not accept:
        question = "Would you like to specify a position restraint?"
        do = FUNCTIONS.answer_question(question, "bool")
        if do:
            GET.position_restraint(restraints, experiment)
            continue
        question = "Are you certain you are happy with the position restraints "
        question += "you have specified?"
        accept = FUNCTIONS.answer_question(question, "bool")
    return restraints

def distance_restraints(experiment):
    """Ask the user for distance restraints"""
    # Store them here
    restraints = []
    # Use a while loop
    accept = False
    while not accept:
        question = "Would you like to specify a distance restraint?"
        do = FUNCTIONS.answer_question(question, "bool")
        if do:
            GET.distance_restraint(restraints, experiment)
            continue
        question = "Are you certain you are happy with the distance restraints "
        question += "you have specified?"
        accept = FUNCTIONS.answer_question(question, "bool")
    return restraints

def dihedral_restraints(experiment):
    """Ask the user for dihedral restraints"""
    # Store them here
    restraints = []
    # Use a while loop
    accept = False
    while not accept:
        question = "Would you like to specify a dihedral restraint?"
        do = FUNCTIONS.answer_question(question, "bool")
        if do:
            GET.dihedral_restraint(restraints, experiment)
            continue
        question = "Are you certain you are happy with the dihedral restraints "
        question += "you have specified?"
        accept = FUNCTIONS.answer_question(question, "bool")
    return restraints

def Restraints(experiment):
    """Ask the user for structure restraints"""
    # Clear the screen and say what is happening
    os.system("clear")
    message = """
You now have the OPTION of specifying various types of structure restraints to
ensure that your system behaves how it should. Please note that this is an
expert-level application of the IPRO Suite. Due to their complexity, it is very
possible to cause significant problems in your experiment by creating
inappropriate restraints. Please use extreme caution and careful consideration.

You will have four options:
1) Permanently fixing Atoms in place. When you do this, those Atoms will never,
ever be modified by IPRO in any way.
2) Position restraints - make sure that one or more Atoms stay close to their
initial positions / structure
3) Distance restraints - restrain the distance between two Atoms
4) Dihedral restraints - restrain the dihedral angle formed by four Atoms

Please refer to the provided documentation or additional information."""
    print screen_formatting(message[1:])
    # Get fixed atoms 
    fixed = fixedAtoms(experiment)
    position = position_restraints(experiment)
    distance = distance_restraints(experiment)
    dihedral = dihedral_restraints(experiment)
    # Store that information in a dictionary
    restraints = {}
    if fixed != {}:
        restraints["Fixed Atoms"] = fixed
    if position != []:
        restraints["Position"] = position
    if distance != []:
        restraints["Distance"] = distance
    if dihedral != []:
        restraints["Dihedral"] = dihedral
    # Only store the information if any restraints were collected
    if restraints != {}:
        experiment["Restraints"] = restraints

def antigen_molecules(experiment):
    """Ask the user to specify antigen molecules"""
    # Clear the screen and say what is happening
    os.system("clear")
    message = """
    Identify the antigen molecules you want to use
    Please refer to the provided documentation for additional information."""
    # And the list of Molecules used in the experiment here
    print screen_formatting(message[1:])
    # Store files of molecules here
    files = {}
    molecules = []
    GET.structures(molecules, files, "Target Molecule", "Target Molecules", \
                   experiment["Force Field"], experiment["File Format"])
    experiment["Molecules"] = molecules

def Epitope(experiment):
    """Ask the user to specify the epitope in a Molecule"""
    # Clear the screen and tell the user what is happening
    os.system("clear")
    message = """
You must now identify the Epitope of the antigen for your experiment. These are amino
acid Residues in Target Molecules that are binding to the antibody.

Please refer to the provided documentation for additional information."""
    print screen_formatting(message[1:])
    # Make a list of the Target  Molecules
    molecules = {}
    order = []
    for data in experiment["Molecules"]:
        if not data[2].design:
            order.append(data[2].name)
            molecules[data[2].name] = data[2]
    if len(order) == 0:
        text = "Somehow there are no Target Molecules"
        raise FUNCTIONS.IPRO_IOError(text)
    # Store the epitope in this dictionary
    positions = {}
    # Use a while loop
    accept = False
    pick = True
    while not accept:
        # If Design Positions should be found
        if pick:
            # Get the Molecule
            if len(order) == 1:
                text = "\nMolecule " + order[0] + " was automatically selected "
                text += "because it is the only Target Molecule."
                print screen_formatting(text)
                mn = order[0]
            else:
                # Summarize the possible Target Molecules
                summary = ''
                for mn in order:
                    summary += "\nTarget Molecule " + mn + " contains "
                    summary += str(len(molecules[mn])) + " Residues"
                print screen_formatting(summary)
                # Ask which Molecule, giving an option of selecting None
                question = "In which Molecule would you like to specify epitope?"
                question += "You may answer 'none' if you would not "
                question += "like to pick any."
                answers = list(order)
                answers.extend(["NONE", 'None', 'none'])
                mn = FUNCTIONS.answer_question(question, "string", [], answers)
                if mn.upper() == "NONE":
                    pick = False
        # if a Molecule has been selected
        if pick:
            GET.EpitopePositions(positions, molecules[mn], experiment)
            # if there are multiple Molecules, give the option of picking more
            # choices
            if len(order) > 1:
                question = "Would you like to specify epitope in "
                question += "another Molecule?"
                pick = FUNCTIONS.answer_question(question, "bool")
            else:
                pick = False
        # If there are no selected Design Positions, give the user an error
        if positions == {}:
            text = "\nYou must select at least one epitope residue"
            print screen_formatting(text)
            pick = True
            continue
        # If a selection either shouldn't or "can't" be made, give the options
        # to be done
        if not pick:
            # Summarize the selected Design Positions
            summary = ''
            for mn in order:
                if mn in positions:
                    summary += "\nTarget Molecule " + mn + ": " + \
                        str(len(positions[mn])) + " epitope selected"
            print screen_formatting(summary)
            # Ask if the user is happy with this
            question = "Are you certain you have specified the correct epitope "
            question += "residues for your experiment?"
            right = FUNCTIONS.answer_question(question, "bool")
            # If the user is happy, be done
            if right:
                accept = True
            # if not, let them modify their choices
            else:
                text = "\nThen modify the selections you've made until you are."
                print screen_formatting(text)
                pick = True
    # Store the selected Design Positions in the Experiment
    experiment["Epitope Positions"] = positions

def for_OPTMAVEN_human_sequences_file(experiment, default):
    """Ask the user for the human 9mer sequences"""
    # Just use the get_list function
    FUNCTIONS.get_list("Human Sequences File", experiment, [], default, \
                       defaultHumanSequencesFile)

def for_OPTMAVEN_integer_cuts_file(experiment, default):
    """Ask the user for the Maps integer cuts file"""
    # Just use the get_list function
    FUNCTIONS.get_list("Maps Integer Cuts File", experiment, [], default, \
                       defaultIntegerCutsFile)

def for_OPTMAVEN_antibody_parts_files(experiment, default):
    """Ask the user for the antibody parts files"""
    # Just use the get_list function
    FUNCTIONS.get_list("Antibody Parts Files", experiment, [], default, \
                       defaultAntibodyParts)

def optmaven_info(experiment):
    """Ask the user questions about how to run OptMAVEn."""
    # Clear the screen and tell the user what's what
    os.system("clear")
    message = """
You now need to provide information about how to run OptMAVEn during
your experiment.

Please refer to the provided documentation for additional information.
    """
    print screen_formatting(message[1:])
    # Find out if defaults should be used
    question = "Would you like to use the default OptMAVEn settings?"
    default = FUNCTIONS.answer_question(question, "bool")
    # Get the values
    for_OPTMAVEN_human_sequences_file(experiment, default)
    for_OPTMAVEN_integer_cuts_file(experiment, default)
    for_OPTMAVEN_antibody_parts_files(experiment, default)
    for_OPTMAVEN_antigen_position(experiment, default)
    for_OPTMAVEN_position_refinement(experiment, default)
    for_OPTMAVEN_(experiment, default)
    for_annealing_sharing(experiment, default)
    for_energy_type(experiment, default)

