#!/usr/bin/env python

# The name of this file
__name__ = "IPRO Suite Functions to Get Attributes"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This file contains functions that get specific attributes using a user's
input."""

# Inlcude standard PYTHON modules
import os
import sys
# Include other I/O modules
import IO_FUNCTIONS as FUNCTIONS
import IO_CHECK as CHECK
# Include those modules
from STANDARDS import *
import MOLECULES

# Start by getting information that is related to the use of Molecule structures
def remove_extra_residues(molecule, hetatms):
    """Give the user the option of removing Residues from a Molecule"""
    # First, begin by making a new Molecule out of the hetatms Residues
    text = ''
    for residue in hetatms:
        text += format(residue, residue.fileFormat)
    new = MOLECULES.Molecule(text, None, None, molecule.name, False, \
                             molecule.forceField, molecule.fileFormat)
    # Figure out if the user wants to exclude any Residues from this list
    exclude = []
    # Use a while loop to get all of the appropriate information
    accept = False
    while not accept:
        # Summarize the new Molecule, leaving out excluded Residues
        print screen_formatting(MOLECULES.summarize(new, exclude))
        # Get the names of the Residues that may still be excluded
        options = []
        for residue in new:
            if residue.name not in exclude:
                options.append(residue.name)
        # If there are no options, don't allow the ability to remove more
        if len(options) == 0:
            do = False
        # Otherwise, ask
        elif len(options) == 1:
            question = "Would you like to exclude Residue " + options[0] 
            question += " from use in the Experiment?"
            do = FUNCTIONS.answer_question(question, "bool")
        else:
            question = "Would you like to exclude one or more of these Residues"
            question += " from use in the Experiment?"
            do = FUNCTIONS.answer_question(question, "bool")
        # If the user wants to exclude Residues
        if do:
            # If there is only one option, remove it automatically
            if len(options) == 1:
                residues = options
            # Otherwise, ask
            else:
                answers = list(options)
                answers.extend(['RANGE', 'Range', 'range'])
                answers.extend(['NONE', 'None', 'none'])
                question = "Which Residue would you like to exclude? You may "
                question += "answer 'range' to exclude a sequential range or "
                question += "'none' if this was a mistake."
                # Get the residue specification
                what = FUNCTIONS.answer_question(question, "string", [],answers)
                # if nothing should be done, change do to False
                if what.upper() == "NONE":
                    do = False
                # If a specific Residue should not be considered
                elif what.upper() != "RANGE":
                    residues = [what]
                # Otherwise ask for the range
                elif len(options) == 2:
                    residues = options
                else:
                    # Get the first Residues
                    answers = options[:-1]
                    question = "What is the first Residue in the sequential "
                    question += "range to be excluded?"
                    rn1 = FUNCTIONS.answer_question(question, "string", [], \
                                                    answers)
                    # Get the index of that Residue
                    i1 = options.index(rn1)
                    # Get the remaining Residues
                    answers = options[i1+1:]
                    if len(answers) == 1:
                        rn2 = answers[0]
                    else:
                        question = "What is the last Residue in the "
                        question += "sequential range to be excluded?"
                        rn2 = FUNCTIONS.answer_question(question, "string", \
                                                        [], answers)
                    # Get that second index
                    i2 = options.index(rn2)
                    # Get the Residues
                    residues = options[i1:i2+1]
        # confirm that those Residues really should be excluded
        if do:
            # Make the summary differently depending on how many Residues are in
            # the residues list
            if len(residues) == 1:
                question = "Are you certain Residue " + residues[0] + " should "
                question += "excluded from use in your experiment?"
            else:
                question = "Are you certain the " + len(residues) + " Residues "
                question += "between and inclusive of " + residues[0] + ' and '
                question += residues[1] + " should be excluded from the "
                question += "experiment?"
            # Get the answer
            use = FUNCTIONS.answer_question(question, "bool")
            # If they should be excluded
            if use:
                # Make a new list so everything stays organized
                newExclude = []
                for residue in new:
                    if residue.name in exclude or residue.name in residues:
                        newExclude.append(residue.name)
                # Use this new list
                exclude = newExclude
        # If the user didn't want to specify residues to exclude, allow them
        # other options
        else:
            # Make a summary of the excluded Residues
            summary = "\nResidues to exclude from use in Molecule " + new.name
            # list objects is from STANDARDS
            summary += ":" + list_items(exclude)
            print screen_formatting(summary)
            # Find out if the user is happy with that
            right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            # If the user is happy, change accept and continue the while loop
            if right:
                accept = True
                continue
            # Otherwise, ask how to proceed
            # Store answers here
            answers = []
            # The types of answers the user may give
            phrases = []
            # If there are more options, allow the user to pick more
            if len(options) > 0:
                phrases.append("select 'more' Residues to exclude")
                answers.extend(['MORE', 'More', 'more'])
            # If there's at least one excluded Residue, the user can start over
            if len(exclude) > 0:
                phrases.append("'start over' at excluding Residues")
                answers.extend(['START OVER', 'Start Over', 'Start over', \
                                'start over'])
            # If there are multiple answers, a single one may be removed
            if len(exclude) > 1:
                phrases.append("'remove' a single Residue from the excluded " \
                               + "list")
                answers.extend(['REMOVE', 'Remove', 'remove'])
            # If there are none of these options, just break the while loop
            # because something WEIRD is happening
            if len(phrases) == 0:
                break
            # Assemble the question
            question = "Would you like to " + phrases[0]
            if len(phrases) == 2:
                question += " or " + phrases[1]
            elif len(phrases) == 3:
                question += ", " + phrases[1] + ", or " + phrases[2]
            question += "?"
            # Get the answer
            what = FUNCTIONS.answer_question(question, "string", [], answers)
            # If the choice is just to do more, continue the while loop
            if what.upper() == "MORE":
                continue
            # If the choice is to start over, reset exclude to empty
            if what.upper() == "START OVER":
                exclude = []
                continue
            # Otherwise, ask which Residue to exclude
            question = "Which Residue should be removed from the excluded list?"
            rn = FUNCTIONS.answer_question(question, "string", [], exclude)
            # Remove that Residue from the list
            i = exclude.index(rn)
            del exclude[i]
    # Now that all excluded Residues are known, assemble a final molecule
    possibles = []
    for residue in new:
        possibles.append(residue.name)
    # Store the text of the final molecule here
    text = ''
    for residue in molecule:
        # If it isn't in the possibles list OR it isn't in the excluded list,
        # than it should be included
        if residue.name not in possibles or residue.name not in exclude:
            text += format(residue, residue.fileFormat)
    # Make a new Molecule
    final = MOLECULES.Molecule(text, None, None, molecule.name, False, \
                               molecule.forceField, molecule.fileFormat)
    # Return that final Molecule
    return final

def structure(molecules, files, tag, forceField, fileFormat, none = False):
    """Retrieve a specific Molecule for use in an Experiment"""
    # Use a while loop so the user can start over if there is a problem
    accept = False
    while not accept:
        # Get the file that contains the relevant Molecule
        question = "What is the name of the file containing the " + tag + "?"
        # If the user is allowed to say this was a mistake
        if none:
            question += " You may answer 'none' if this was a mistake."
        fileName = FUNCTIONS.answer_question(question, "string")
        # If the user said 'none'
        if none and fileName.upper() == "NONE":
            break
        # Try to access the specified file
        if fileName not in files:
            try:
                f = MOLECULES.MoleculeFile(fileName, "./", forceField, \
                                           fileFormat)
            except MOLECULES.MoleculeError as error:
                print str(error)
                continue
            # Store the file in the files dictionary
            files[fileName] = f
        # Access the appropriate structure
        FILE = files[fileName]
        # Determine if any previous Molecules from this file have been specified
        exclude = []
        for data in molecules:
            if data[0] == fileName:
                exclude.append(data[1])
        # Get the structures from that file that CAN be used
        answers = []
        for mn in FILE.moleculeNames:
            if mn not in exclude:
                answers.append(mn)
        # If there are no possible answers, tell the user and continue the while
        # loop
        if len(answers) == 0:
            text = "\nAll Molecules in " + fileName + " have already been used"
            # Screen formatting is from STANDARDS
            print screen_formatting(text)
            continue
        # Summarize the contents of the file
        print screen_formatting(MOLECULES.summarize(FILE, exclude))
        # Figure out which Molecule the user wants
        if len(answers) == 1:
            mn = answers[0]
            text = "\nMolecule " + mn + " was the only remaining Molecule in "
            text += fileName + ", so it was chosen automatically."
            print screen_formatting(text)
        # If there are multiple possible answers, ask the user
        else:
            question = "Which Molecule would you like to select as the " + tag
            question += "?"
            mn = FUNCTIONS.answer_question(question, "string", [], answers)
        # Get a duplicate of the File's Molecule so that it isn't modified
        molecule = FILE[mn].duplicate()
        # Check that molecule's structure
        missing, gaps, hetatms = CHECK.structure(molecule)
        # Determine if there is a problem
        problem = ''
        for rn in missing:
            problem += "\nResidue " + rn + " is incomplete"
        for line in gaps:
            problem += "\n" + line
        # If there is a problem, give it an appropriate header
        if problem != '':
            problem = "\nThe following structural problems were identified:" + \
                      problem
            print screen_formatting(problem)
        # Confirm that the user wants to use this Molecule anyway
        question = "Are you certain you want to use this Molecule?"
        use = FUNCTIONS.answer_question(question, "bool")
        # If the user doesn't want to use it, continue the while loop
        if not use:
            continue
        # If there are Residues that maybe shouldn't be included, ask the user
        # about it
        if hetatms != []:
            molecule = remove_extra_residues(molecule, hetatms)
        # Get a name for the Molecule
        accept2 = False
        # Get the names of the Molecules that have already been selected
        usedNames = []
        for data in molecules:
            usedNames.append(data[2].name)
        while not accept2:
            question = "What would you like to name this Molecule in your IPRO "
            question += "Suite Experiment?"
            name = FUNCTIONS.answer_question(question, "string", [], [], \
                                             usedNames)
            # Try to use that name
            try:
                molecule.name = name
            except MOLECULES.MoleculeError as error:
                print str(error)
                continue
            # Confirm this is right
            right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            if right:
                accept2 = True
        # Modify the Molecule's Design status based on the tag
        if tag in ['Design Molecule', 'Design', 'Enzyme', 'Protein', \
                   'Antibody']:
            molecule.design = True
        # Store this Molecule and end the while loop
        molecules.append([fileName, mn, molecule, tag])
        accept = True
    # Nothing needs to be returned, because of pass by reference

def structures(molecules, files, single, plural, forceField, fileFormat):
    """Retrieve multiple structures of a particular type."""
    # Use a while loop
    accept = False
    # Have a tag to say whether or not to select a structure at the start of the
    # while loop
    pick = True
    while not accept:
        # If a structure should be selected, do so!
        if pick:
            # Determine if there are any existing Molecules of this type
            none = False
            for data in molecules:
                if data[3] == single:
                    none = True
                    break
            # Get a structure
            structure(molecules, files, single, forceField, fileFormat, none)
            # Ask if they want to pick another
            question = "Would you like to select another " + single + "?"
            another = FUNCTIONS.answer_question(question, "bool")
            if another:
                continue
        # Create a summary of the Structures selected so far
        summary = "\n" + plural + " selected so far:"
        # Make a list of the structures already selected that match the tag
        items = []
        for data in molecules:
            if data[3] == single:
                text = data[2].name + " (" + data[1] + " from " + data[0] + ")"
                items.append(text)
        # list items is from STANDARDS
        summary += list_items(items)
        print screen_formatting(summary)
        # If no Molecules of this type have been selected, try again
        if len(items) == 0:
            text = "\nYou must select at least one " + single
            print screen_formatting(text)
            pick = True
            continue
        # Confirm that the user is happy with their choices
        question = "Are you certain you have selected the correct "+plural+"?"
        right = FUNCTIONS.answer_question(question, "bool")
        # If the user is happy, be done with this while loop
        if right:
            accept = True
            continue
        # Otherwise, ask how to modify things
        answers = ['MORE', 'More', 'more']
        phrases = ["select 'more' " + plural]
        # If it is possible to start over at selecting the Molecules
        if len(items) > 0:
            answers.extend(['START OVER', 'Start Over', 'Start over', \
                            'start over'])
            phrases.append("'start over' at selecting " + plural)
        # If there are multiple selections, one may be removed
        if len(items) > 1:
            answers.extend(['REMOVE', 'Remove', 'remove'])
            phrases.append("'remove' a specific " + single)
        # Make the question
        question = "Would you like to " + phrases[0]
        if len(phrases) == 2:
            question += " or " + phrases[1]
        elif len(phrases) == 3:
            question += ", " + phrases[1] + ", or " + phrases[2]
        question += "?"
        if len(phrases) == 1:
            what = answers[0]
        else:
            what = FUNCTIONS.answer_question(question, "string", [], answers)
        # Respond appropriately. If more are wanted, do that
        if what.upper() == "MORE":
            pick = True
            continue
        if what.upper() == "START OVER":
            pick = True
            # Delete every instance in the Molecules of this type being selected
            sorted = False
            while not sorted:
                for i in range(len(molecules)):
                    if molecules[i][3] == single:
                        del molecules[i]
                        sorted = False
                        break
            # continue the while loop
            continue
        # Otherwise, ask which Molecule to remove
        question = "Which Molecule should be removed from the Experiment?"
        answers = []
        for data in molecules:
            if data[3] == single:
                answers.append(data[2].name)
        # Get that Molecule's name
        mn = FUNCTIONS.answer_question(question, "string", [], answers)
        # Remove that Molecule
        for i in range(len(molecules)):
            if molecules[i][2].name == mn:
                del molecules[i]
                break
        # In this instance only do not automatically select another Molecule the
        # next time through this for loop
        pick = False
    # Nothing needs to be returned thanks to pass by reference

def Dimers(experiment):
    """Search through an Experiment's Molecules to find Dimers"""
    # Store the calculated dimers here
    dimers = []
    # go through the Molecules
    for i in range(len(experiment["Molecules"]) - 1):
        # If it is a Design Molecule, use it
        if experiment["Molecules"][i][2].design:
            mol1 = experiment["Molecules"][i][2]
        # If it isn't, skip it
        else:
            continue
        # Check it against every other Molecule
        for j in range(i+1, len(experiment["Molecules"])):
            if experiment["Molecules"][j][2].design:
                mol2 = experiment["Molecules"][j][2]
            else:
                continue
            # Determine which Molecule is longer
            if len(mol1) >= len(mol2):
                molA = mol1
                molB = mol2
            else:
                molA = mol2
                molB = mol1
            # Calculate the sequence similarity between these Design Molecules
            same = 0
            for residue in molA:
                if residue.name in molB and residue.kind == \
                molB[residue.name].kind:
                    same += 1
            # Calculate a percent
            percent = 100.0 * float(same)/float(len(molA))
            # Only do anything if the percent sequence similarity is > 95%
            if percent >= 95.0:
                question = "Design Molecules " + mol1.name + " and " + mol2.name
                question += " have a sequence similarity of "
                question += format(percent, '.2f') + "%. Are they dimers that "
                question += "should always have the same sequence?"
                # Use a while loop
                accept = False
                while not accept:
                    same = FUNCTIONS.answer_question(question, "bool")
                    right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                    if right:
                        accept = True
                # If they are dimers, store that
                if same:
                    dimers.append([mol1.name, mol2.name])
    # Store the Dimers in the experiment
    experiment["Dimers"] = dimers

def DesignPositions(positions, molecule, experiment):
    """Select Design Positions in a Molecule"""
    # If there are previously selected Design Positions for this Molecule, get
    # them
    if molecule.name in positions:
        chosen = positions[molecule.name]
    else:
        chosen = []
    # Use a while loop
    accept = False
    # Have a tag to say whether or not to make more picks
    pick = True
    while not accept:
        # Make a list of the acceptable Residues as well as the ones to ignore
        options = []
        ignore = []
        for residue in molecule:
            if residue.name in chosen:
                ignore.append(residue.name)
            elif residue.kind not in aminoAcids[residue.fileFormat]:
                ignore.append(residue.name)
            else:
                options.append(residue.name)
        # If appropriate, pick one or more Residues as Design Positions
        if pick and len(options) == 0:
            text = "\nThere are no remaining Residues in Molecule "
            text += molecule.name + " that may be Design Positions."
            print screen_formatting(text)
            pick = False
        elif pick and len(options) == 1:
            text = "\nResidue " + options[0] + " is the only remaining Residue "
            text += "in Molecule " + molecule.name + ", so it was chosen "
            text += "automatically."
            print screen_formatting(text)
            residues = options
        elif pick:
            # Summarize the Molecule
            print screen_formatting(MOLECULES.summarize(molecule, ignore))
            # Ask the user
            answers = list(options)
            question = "What Residue is a Design Position? You may answer "
            question+= "'none' if you would not like to pick a Design Position"
            answers.extend(['NONE', 'None', 'none'])
            if len(options) > 1:
                question += " or you may answer 'range' if you would like to "
                question += "select a sequential set of Residues."
                answers.extend(['RANGE', 'Range', 'range'])
            choice = FUNCTIONS.answer_question(question, "string", [], answers)
            # Respond to this choice appropriately. If they don't want to pick
            # anything say that is what happened
            if choice.upper() == "NONE":
                pick = False
            # If they picked a specific Residue, store that
            elif choice.upper() != "RANGE":
                residues = [choice]
            # Otherwise they picked a range. If there's only two choices, use
            # them
            elif len(options) == 2:
                text = "\nResidues " + options[0] + " and " + options[1]
                text += " are the only options, so they were chosen "
                text += "automatically"
                print screen_formatting(text)
                residues = options
            # Otherwise, they must be asked
            else:
                # Get the first Residue
                answers = options[:-1]
                question = "What is the first Residue in the sequential range "
                question += "of Design Positions?"
                rn1 = FUNCTIONS.answer_question(question, "string", [], answers)
                # Get the index of that Residue
                i1 = options.index(rn1)
                # Get the second Residue
                answers = options[i1+1:]
                if len(answers) == 1:
                    rn2 = answers[0]
                else:
                    question = "What is the last Residue in the sequential "
                    question += "range of Design Positions?"
                    rn2 = FUNCTIONS.answer_question(question, "string", [], \
                                                    answers)
                # And get that index
                i2 = options.index(rn2)
                # Get the full list of Residues
                residues = options[i1:i2+1]
        # If a selection was made, confirm it is right
        if pick:
            # Get the question differently depending on how many Residues were
            # picked
            if len(residues) == 1:
                question = "Are you certain Residue " + residues[0] + " should "
                question += "used as a Design Position?"
            else:
                question = "Are you certain the " + str(len(residues))
                question += " Residues between and inclusive of Residues "
                question += residues[0] + " and " + residues[-1] + " should be "
                question += "used as Design Positions?"
            # Confirm this is correct
            right = FUNCTIONS.answer_question(question, "bool")
            # If the choice is correct, store the new Residues in the chosen
            # list - but do so in a way that maintains the order of the Residues
            # from the Molecule
            if right:
                newChosen = []
                for residue in molecule:
                    # If the Residue has ever been selected, store it
                    if residue.name in chosen or residue.name in residues:
                        newChosen.append(residue.name)
                # Replace the chosen list
                chosen = newChosen
            # Find out if the user would like to pick more Design Positions
            if len(residues) != len(options):
                question = "Would you like to select more Design Positions from"
                question += " Molecule " + molecule.name + "?"
                more = FUNCTIONS.answer_question(question, "bool")
                if more:
                    # pick must be True, so it'll restart that process
                    continue
                else:
                    pick = False
            # If there are no remaining options, don't go that route
            else:
                pick = False
        # If Residues either weren't selected or couldn't be
        if not pick:
            # Summarize the selected Design Positions
            summary = "\nSelected Design Positions in Molecule " + molecule.name
            summary += ":" + list_items(chosen)
            print screen_formatting(summary)
            # Find out if the user is happy with this
            question = "Are you certain these are the correct Design Positions "
            question += "in Molecule " + molecule.name + "?"
            right = FUNCTIONS.answer_question(question, "bool")
            # If the user is certain this is correct, be done
            if right:
                accept = True
                continue
            # Assemble options for what to do
            answers = []
            phrases = []
            # If there are remaining options
            if len(options) > 0:
                phrases.append("select 'more' Design Positions")
                answers.extend(['MORE', 'More', 'more'])
            # If there is at least one position selected
            if len(chosen) > 0:
                phrases.append("'start over' at picking Design Positions in " \
                               + "Molecule " + molecule.name)
                answers.extend(['START OVER', 'Start Over', 'Start over', \
                                'start over'])
            # If there are multiple selections
            if len(chosen) > 1:
                phrases.append("'remove' a specific Design Position")
                answers.extend(['REMOVE', 'Remove', 'remove'])
            # If there are no choices, just break the while loop because that is
            # WEIRD
            if len(phrases) == 0:
                break
            elif len(phrases) == 1:
                what = answers[0]
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1] + "?"
                else:
                    question += ", " + phrases[1] + ", or " + phrases[2] + "?"
                what = FUNCTIONS.answer_question(question, "string", [],answers)
            # Respond appropriately to those options
            if what.upper() == "MORE":
                pick = True
                continue
            if what.upper() == "START OVER":
                chosen = []
                pick = True
                continue
            # Otherwise, a particular Residue should be removed
            question = "Which Residue should NOT be a Design Position?"
            rn = FUNCTIONS.answer_question(question, "string", [], chosen)
            # Get the index and delete that choice
            i = chosen.index(rn)
            del chosen[i]
            # In this case, don't automatically make a pick the next time
            # through the while loop
    # now that the Design Positions in this Molecule have been chosen, modify
    # the positions dictionary
    if len(chosen) == 0:
        # If there are no selected positions, but there USED to be, delete the
        # previous entry
        if molecule.name in positions:
            del positions[molecule.name]
    else:
        positions[molecule.name] = chosen
    # If there are Dimers, make sure their information matches
    if "Dimers" in experiment:
        for pair in experiment["Dimers"]:
            if molecule.name in pair:
                if molecule.name == pair[0]:
                    mn2 = pair[1]
                else:
                    mn2 = pair[0]
                if molecule.name in positions:
                    positions[mn2] = positions[molecule.name]
                elif mn2 in positions:
                    del positions[mn2]
    # Positions does not need to be returned because it was passed by reference

def EpitopePositions(positions, molecule, experiment):

    """Select Epitope Positions in a Molecule"""
    # If there are previously selected Epitope Positions for this Molecule, get
    # them
    if molecule.name in positions:
        chosen = positions[molecule.name]
    else:
        chosen = []
    # Use a while loop
    accept = False
    # Have a tag to say whether or not to make more picks
    pick = True
    while not accept:
        # Make a list of the acceptable Residues as well as the ones to ignore
        options = []
        ignore = []
        for residue in molecule:
            if residue.name in chosen:
                ignore.append(residue.name)
            elif residue.kind not in aminoAcids[residue.fileFormat]:
                ignore.append(residue.name)
            else:
                options.append(residue.name)
        # If appropriate, pick one or more Residues as Epitope Positions
        if pick and len(options) == 0:
            text = "\nThere are no remaining Residues in Molecule "
            text += molecule.name + " that may be Epitope Positions."
            print screen_formatting(text)
            pick = False
        elif pick and len(options) == 1:
            text = "\nResidue " + options[0] + " is the only remaining Residue "
            text += "in Molecule " + molecule.name + ", so it was chosen "
            text += "automatically."
            print screen_formatting(text)
            residues = options
        elif pick:
            # Summarize the Molecule
            print screen_formatting(MOLECULES.summarize(molecule, ignore))
            # Ask the user
            answers = list(options)
            question = "What Residue is a Epitope Position? You may answer "
            question+= "'none' if you would not like to pick a Design Position"
            answers.extend(['NONE', 'None', 'none'])
            if len(options) > 1:
                question += " or you may answer 'range' if you would like to "
                question += "select a sequential set of Residues."
                answers.extend(['RANGE', 'Range', 'range'])
            choice = FUNCTIONS.answer_question(question, "string", [], answers)
            # Respond to this choice appropriately. If they don't want to pick
            # anything say that is what happened
            if choice.upper() == "NONE":
                pick = False
            # If they picked a specific Residue, store that
            elif choice.upper() != "RANGE":
                residues = [choice]
            # Otherwise they picked a range. If there's only two choices, use
            # them
            elif len(options) == 2:
                text = "\nResidues " + options[0] + " and " + options[1]
                text += " are the only options, so they were chosen "
                text += "automatically"
                print screen_formatting(text)
                residues = options
            # Otherwise, they must be asked
            else:
                # Get the first Residue
                answers = options[:-1]
                question = "What is the first Residue in the sequential range "
                question += "of Epitope Positions?"
                rn1 = FUNCTIONS.answer_question(question, "string", [], answers)
                # Get the index of that Residue
                i1 = options.index(rn1)
                # Get the second Residue
                answers = options[i1+1:]
                if len(answers) == 1:
                    rn2 = answers[0]
                else:
                    question = "What is the last Residue in the sequential "
                    question += "range of Epitope Positions?"
                    rn2 = FUNCTIONS.answer_question(question, "string", [], \
                                                    answers)
                # And get that index
                i2 = options.index(rn2)
                # Get the full list of Residues
                residues = options[i1:i2+1]
        # If a selection was made, confirm it is right
        if pick:
            # Get the question differently depending on how many Residues were
            # picked
            if len(residues) == 1:
                question = "Are you certain Residue " + residues[0] + " should "
                question += "used as a Epitope Position?"
            else:
                question = "Are you certain the " + str(len(residues))
                question += " Residues between and inclusive of Residues "
                question += residues[0] + " and " + residues[-1] + " should be "
                question += "used as Epitope Positions?"
            # Confirm this is correct
            right = FUNCTIONS.answer_question(question, "bool")
            # If the choice is correct, store the new Residues in the chosen
            # list - but do so in a way that maintains the order of the Residues
            # from the Molecule
            if right:
                newChosen = []
                for residue in molecule:
                    # If the Residue has ever been selected, store it
                    if residue.name in chosen or residue.name in residues:
                        newChosen.append(residue.name)
                # Replace the chosen list
                chosen = newChosen
            # Find out if the user would like to pick more Epitope Positions
            if len(residues) != len(options):
                question = "Would you like to select more Epitope Positions from"
                question += " Molecule " + molecule.name + "?"
                more = FUNCTIONS.answer_question(question, "bool")
                if more:
                    # pick must be True, so it'll restart that process
                    continue
                else:
                    pick = False
            # If there are no remaining options, don't go that route
            else:
                pick = False
        # If Residues either weren't selected or couldn't be
        if not pick:
            # Summarize the selected Epitope Positions
            summary = "\nSelected Epitope Positions in Molecule " + molecule.name
            summary += ":" + list_items(chosen)
            print screen_formatting(summary)
            # Find out if the user is happy with this
            question = "Are you certain these are the correct Epitope Positions "
            question += "in Molecule " + molecule.name + "?"
            right = FUNCTIONS.answer_question(question, "bool")
            # If the user is certain this is correct, be done
            if right:
                accept = True
                continue
            # Assemble options for what to do
            answers = []
            phrases = []
            # If there are remaining options
            if len(options) > 0:
                phrases.append("select 'more' Epitope Positions")
                answers.extend(['MORE', 'More', 'more'])
            # If there is at least one position selected
            if len(chosen) > 0:
                phrases.append("'start over' at picking Epitope Positions in " \
                               + "Molecule " + molecule.name)
                answers.extend(['START OVER', 'Start Over', 'Start over', \
                                'start over'])
            # If there are multiple selections
            if len(chosen) > 1:
                phrases.append("'remove' a specific Design Position")
                answers.extend(['REMOVE', 'Remove', 'remove'])
            # If there are no choices, just break the while loop because that is
            # WEIRD
            if len(phrases) == 0:
                break
            elif len(phrases) == 1:
                what = answers[0]
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1] + "?"
                else:
                    question += ", " + phrases[1] + ", or " + phrases[2] + "?"
                what = FUNCTIONS.answer_question(question, "string", [],answers)
            # Respond appropriately to those options
            if what.upper() == "MORE":
                pick = True
                continue
            if what.upper() == "START OVER":
                chosen = []
                pick = True
                continue
            # Otherwise, a particular Residue should be removed
            question = "Which Residue should NOT be a Design Position?"
            rn = FUNCTIONS.answer_question(question, "string", [], chosen)
            # Get the index and delete that choice
            i = chosen.index(rn)
            del chosen[i]
            # In this case, don't automatically make a pick the next time
            # through the while loop
    # now that the Epitope Positions in this Molecule have been chosen, modify
    # the positions dictionary
    if len(chosen) == 0:
        # If there are no selected positions, but there USED to be, delete the
        # previous entry
        if molecule.name in positions:
            del positions[molecule.name]
    else:
        positions[molecule.name] = chosen
    # If there are Dimers, make sure their information matches
    if "Dimers" in experiment:
        for pair in experiment["Dimers"]:
            if molecule.name in pair:
                if molecule.name == pair[0]:
                    mn2 = pair[1]
                else:
                    mn2 = pair[0]
                if molecule.name in positions:
                    positions[mn2] = positions[molecule.name]
                elif mn2 in positions:
                    del positions[mn2]
    # Positions does not need to be returned because it was passed by reference

def all_combinations(items):
    """Make a list of all possible combinations of a set of items"""
    # Store the list here
    all = []
    # Start with each item by itself
    for item in items:
        all.append([item])
    # Do this iteratively based on the number of total items
    levels = {1:all}
    for i in range(2, len(items) + 1):
        # Initialize the list for this number of items
        levels[i] = []
        # Go through the previous lists
        for old in levels[i-1]:
            # Loop through the items
            for item in items:
                # Skip it if it is already in this list
                if item in old:
                    continue
                # Make a new list and store the item in it
                new = list(old)
                new.append(item)
                # Sort it
                new.sort()
                # If it is a unique list, store it
                if new not in levels[i]:
                    levels[i].append(new)
        # Store the lists with this many items
        all.extend(levels[i])
    return all

def DesignGroup(groups, molecules):
    """Identify a Design Group"""
    # Use a while loop
    accept = False
    # Store the chosen Molecules' names here
    chosen = []
    # Have a flag to say whether or not to ask for another Molecule
    pick = True
    while not accept:
        # Make a list of the Molecules that may be selected, as well as those
        # that should be ignored
        options = []
        ignore = []
        for molecule in molecules:
            if molecule.design:
                ignore.append(molecule.name)
            elif molecule.name in chosen:
                ignore.append(molecule.name)
            else:
                options.append(molecule.name)
        # If a Molecule should be selected
        if pick:
            # If there are no Molecules, say that
            if len(options) == 0:
                text = "\nThere are no remaining Molecules that may be part "
                text += "of the Design Group."
                print screen_formatting(text)
                pick = False
            # if there is only one, choose it automatically
            elif len(options) == 1:
                text = "\nMolecule " + options[0] +" was selected automatically"
                text += " because it was the only option."
                print screen_formatting(text)
                mn = options[0]
            else:
                # Summarize the choices
                for molecule in molecules:
                    if molecule.name not in ignore:
                        print screen_formatting(MOLECULES.summarize(molecule))
                # Ask the user
                question = "Which Target Molecule should be included in the "
                question += "Design Group?"
                mn = FUNCTIONS.answer_question(question, "string", [], options)
        # If a choice was made, confirm that it is correct
        if pick:
            question = "Are you certain Molecule " + mn + " should be used in "
            question += "the Design Group?"
            right = FUNCTIONS.answer_question(question, "bool")
            # If they're certain, store the Molecule
            if right:
                chosen.append(mn)
                chosen.sort()
            # Ask if another should be selected
            if len(options) > 1:
                question = "Would you like to select another Target Molecule "
                question += "for the Design Group?"
                more = FUNCTIONS.answer_question(question, "bool")
                if more:
                    continue
                else:
                    pick = False
            else:
                pick = False
        # If a selection wasn't made or couldn't be made or shouldn't be made
        # next, check the overall Design Group
        if not pick:
            # Summarize the contents of the Design Group
            summary = "\nSelected Molecules for the Design Group:"
            summary += list_items(chosen)
            print screen_formatting(summary)
            # Confirm that this summary is correct
            question = "Are you certain the list of Molecules in this Design "
            question += "Group is correct?"
            right = FUNCTIONS.answer_question(question, "bool")
            if right:
                accept = True
                continue
            # If the user isn't happy, start accumulating options
            answers = []
            phrases = []
            # If there are more possible options
            if len(options) > 0:
                answers.extend(['MORE', 'More', 'more'])
                phrases.append("select 'more' Molecules for the Design Group")
            # If there is at least one choice made so far
            if len(chosen) > 0:
                answers.extend(['START OVER', 'Start Over', 'Start over', \
                                'start over'])
                phrases.append("'start over' at selecting Molecules for the " \
                               + "Design Group")
            # If there are at least two choices so far
            if len(chosen) > 1:
                answers.extend(['REMOVE', 'Remove', 'remove'])
                phrases.append("'remove' a specific Molecule from the Group")
            # If there are no options, break the while loop because that is
            # WEIRD
            if len(phrases) == 0:
                break
            # If there's only one option, pick it
            elif len(phrases) == 1:
                what = answers[0]
            # Otherwise, ask
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1] + "?"
                else:
                    question += ", " + phrases[1] + ", or " + phrases[2] + "?"
                what = FUNCTIONS.answer_question(question, "string", [], \
                                                 answers)
            # Respond appropriately to the request
            if what.upper() == "MORE":
                pick = True
                continue
            if what.upper() == "START OVER":
                chosen = []
                pick = True
                continue
            # Otherwise the choice is to remove a Molecule. In that case DON'T
            # pick next time through. For now, find out which Molecule to remove
            question = "Which Molecule should be removed from the Design Group?"
            mn = FUNCTIONS.answer_question(question, "string", [], chosen)
            i = chosen.index(mn)
            del chosen[i]
    # At this point, the user is happy with the Design Group they've specified.
    # Evaluate it
    if len(chosen) != 0:
        # When the molecules are stored in chosen, they are sorted
        # alphabetically, so this check can be done directly
        have = False
        for group in groups:
            if chosen == group[1:]:
                have = True
                break
        # If the group has been previously specified
        if have:
            text = "\nThat selection of Molecules has already been picked as a "
            text += "Design Group and cannot be reused."
            print screen_formatting(text)
        else:
            # Ask how to modify the binding
            accept = False
            while not accept:
                question = "Would you like to 'improve', 'maintain', 'reduce', "
                question += "or 'eliminate' binding to this Design Group?"
                binding = FUNCTIONS.answer_question(question, "string", [], \
                          ['improve', 'maintain', 'reduce', 'eliminate'])
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept = True
            # Store the Design Group
            group = [binding]
            group.extend(chosen)
            groups.append(group)

def OptzymeGroup(groups, molecules):
    """Identify an Optzyme Group"""
    # Ask the catalytic property being modified 
    accept = False
    while not accept:
        question = "Would you like to modify 'km', 'kcat/km', or 'kcat' in "
        question += "this Optzyme Group?"
        property = FUNCTIONS.answer_question(question, "string", [], ['km', \
                   'kcat/km', 'kcat'])
        # Confirm that the user has answered properly
        right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
        if right:
            accept = True
    # Store the binding assemblies in a list
    assemblies = []
    # Ask the user the Target Molecules which molecules should be included in
    # the ground state and TSA state calculations, whichever are applicable
    for state in ["ground", "TSA"]:
        # If this is a kcat/km calculation, there are no ground state
        # calculations. Similarly, there are no TSA state calculations in a km
        # calculation
        if [property, state] in [["kcat/km", "ground"], ["km", "TSA"]]:
            continue
        # Use a while loop for the binding assemblies
        accept = False
        # Store the chosen Molecules' names here
        chosen = []
        # Have a flag to say whether or not to ask for another Molecule
        pick = True
        while not accept:
            # Make a list of the Molecules that may be selected, as well as
            # those that should be ignored
            options = []
            ignore = []
            for molecule in molecules:
                if molecule.design:
                    ignore.append(molecule.name)
                elif molecule.name in chosen:
                    ignore.append(molecule.name)
                else:
                    options.append(molecule.name)
            # If a Molecule should be selected
            if pick:
                # If there are no Molecules, say that
                if len(options) == 0:
                    text = "\nThere are no remaining Molecules that may be part"
                    text += " of this " + state + " state Binding Assembly."
                    print screen_formatting(text)
                    pick = False
                # If there is only one, choose it automatically      
                elif len(options) == 1:
                    text = "\nMolecule " + options[0] + " was selected "
                    text += "automatically because it was the only option "
                    text += "remaining for the " + state + " state Binding "
                    text += "Assembly."
                    print screen_formatting(text)
                    mn = options[0]
                # Otherwise, list the options for the user
                else:
                    for molecule in molecules:
                        if molecule.name not in ignore:
                            print screen_formatting( \
                                              MOLECULES.summarize(molecule))
                    # Ask the user which Target Molecule they would like to
                    # select
                    question = "Which Target Molecule should be included in the"
                    question += " " + state + " state Binding Assembly?"
                    mn = FUNCTIONS.answer_question(question, "string", [], \
                                                   options)
            # If a choice was made, verify it by the user
            if pick:                                     
                question = "Are you certain you would like to use Molecule " 
                question += mn + " in the " + state + " state Binding Assembly?"
                right = FUNCTIONS.answer_question(question, "bool")
                # If they are certain, store the Molecule
                if right:
                    chosen.append(mn)
                    chosen.sort()
                # Ask the user if another Molecule should be selected
                if len(options) > 1:
                    question =  "Would you like to select another Target "
                    question += "Molecule for the " + state + " state "
                    question += "Binding Assembly?"
                    more = FUNCTIONS.answer_question(question, "bool")
                    if more:
                        continue
                    else:
                        pick = False
                else:
                    pick = False
            # If a selection wasn't or couldn't be made, check the overall
            # Binding Assembly
            if not pick:
                # Summarize the contents of the Binding Assembly
                summary = "\nSelected Molecules for the " + state + " state "
                summary += "Binding Assembly:"
                summary += list_items(chosen)
                print screen_formatting(summary)
                # Confirm that this summary is correct
                question = "Are you sure that the list of Target Molecules in "
                question += "this Binding Assembly for the " + state + " state"
                question += " is correct?"
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    continue
                # If the user isn't satisfied, start accumulating options
                answers = []
                phrases = []
                # If there are more possible options
                if len(options) > 0:
                    answers.extend(['MORE', 'More', 'more'])
                    phrases.append("select 'more' Molecules for the " + \
                                  state + " state Binding Assembly")
                # If there is at least one choice made so far
                if len(chosen) > 0:
                    answers.extend(['START OVER', 'Start Over', 'Start over', \
                                   'start over'])
                    phrases.append("'start over' at selecting Molecules " + \
                              "for the " + state + " state Binding Assembly")
                # If there are at least two choices so far
                if len(chosen) > 1:
                    answers.extend(['REMOVE', 'Remove', 'remove'])
                    phrases.append("'remove' a specific Molecule from the " + \
                                   state + " state Binding Assembly")
                # If there are no options, break the while loop because that is
                # bizarre
                if len(phrases) == 0:
                    break
                # If there is only one option, pick it
                elif len(phrases) == 1:
                    what = answers[0]
                # Otherwise, ask
                else:
                    question = "Would you like to " + phrases[0]
                    if len(phrases) == 2:
                        question += " or " + phrases[1] + "?"
                    else:
                        question += ", " + phrases[1] + ", or " + phrases[2]
                        question += "?"
                    # Ask the user this question
                    what = FUNCTIONS.answer_question(question, "string", [], \
                                                     answers)
                # Respond appropriately to the request
                # Add more Target Molecules to the Binding Assembly
                if what.upper() == "MORE":
                    pick = True
                # Start listing the Target Molecules from scratch    
                elif what.upper() == "START OVER":
                    chosen = []
                    pick = True
                # Otherwise, the choice it to remove a Target Molecule. In that
                # case, DON'T pick next time through. For now, find out which
                # Target Molecule to remove
                else:
                    question = "Which Target Molecule should be removed from "
                    question += "the " + state + " state Binding Assembly?"
                    mn = FUNCTIONS.answer_question(question, "string", [], \
                                                   chosen)
                    i = chosen.index(mn)
                    del chosen[i]
        # At this point, the user is happy with the Binding Assembly they've
        # specified. Evaluate it
        if len(chosen) != 0:
            # When the molecules are stored in chosen, they are stored
            # alphabetically so a check for duplicates can be performed directly
            have = False
            for assembly in assemblies:
                if chosen == assembly[1:]:
                    have = True
                    break
            # If both the ground and TSA states used the same set of Target
            # Molecules during kcat calculations
            if have:
                text = "\nThat selection of Target Molecules has already been "
                text += "chosen for the ground state Binding Assembly for "
                text += "kcat calculation and it cannot be reused."
                print screen_formatting(text)
            else:
                # Since the objective of a Binding Assembly is irrelevant within
                # OptZyme, use the Binding Assembly objective to distinguish
                # between ground and TSA states 
                if state == "ground":
                    binding = "eliminate"
                else:
                    binding = "improve"
                # Store the Binding Assembly
                assembly = [binding]
                assembly.extend(chosen)
                assemblies.append(assembly)
    # Evaluate the validity of the inputs provided for the Optzyme Group(s)
    if len(assemblies) != 0:
        # Check the number of Binding Assemblies in the Optzyme Group
        if (property == "kcat" and len(assemblies) != 2) or (property in \
            ["km", "kcat/km"] and len(assemblies) != 1):
                text = "An improper number of Binding Assemblies have been "
                text += "selected for this Optzyme Group based on the "
                text += "catalytic property that is being modified."
                print screen_formatting(text)
        # When the molecules are stored, they are retained alphabetically so
        # this check can be done correctly
        have = False
        for group in groups:
            if len(assemblies) in [1, 2]:           
                # Check for a duplicate Optzyme Group
                if group[3] == assemblies:
                    have = True
                    break
                # The indices in the assemblies list must be 0 or 1
                for i in [0, 1]:
                    # The same goes for the assmblies list within 'group'
                    for j in [0, 1]:
                        # If the index does not exist in either list, skip it
                        if (len(group[3]) == 1 and j == 1) or \
                           (len(assemblies) == 1 and i == 1):
                            continue
                        # Check for a duplicate set of Target Molecules and a
                        # different state (ground or TSA) for the set of
                        # molecules. The same set of Target Molecules can be
                        # used as a single state if used in both kcat and
                        # either km or kcat/km calculations
                        elif (assemblies[i][1:] ==  group[3][j][1:]) and \
                             (assemblies[i][0] != group[3][j][0]):
                                have = True
                                break
        # Notify the user if an error was found
        if have:
            text = "This Optzyme Group has already been specified or the set "
            text += "of Target Molecules within it was specified as a ground "
            text += "state in one OptzymeGroup and a TSA state in another."
            print screen_formatting(text)
        else:      
            # If this is a kcat experiment, ask the weighting factor
            if property == "kcat":
                # Keep asking the question until an acceptable answer is
                # provided
                while True:
                    question = "What value should be used to weight the "
                    question += "difference in ground and TSA state interaction" 
                    question += " energies, (RT)TSA/(RT)S? The default value is"
                    question += " " + str(defaultWeight) + "."
                    weight = FUNCTIONS.answer_question(question, "float")
                    # The value of 'weight' must be a positive number
                    if float(weight) > 0.0:
                        question = "Are you sure you would like to use a "
                        question += "weight of " + str(weight) + " "
                        question += "for this Optzyme Group?"
                        # Confirm that the user is satisfied with the answer
                        right = FUNCTIONS.answer_question(question, "bool")
                        if right:
                            break
                    # Make the user choose a positive number
                    else:
                        text = "You must answer with a positive value. Please "
                        text += "try again."
                        print screen_formatting(text)
            # Otherwise, just store a dummy value since it is not used
            else:
                weight = defaultWeight
            # Ask how to modify the binding
            accept = False
            while not accept:
                question = "Would you like to 'improve', 'maintain', 'reduce',"
                question += " or 'eliminate' " + property + " for this Optzyme"
                question += " group. Remember that improving " + property 
                question += " indicates that the numerical value is "
                # Remind the user that an improvement in km means a lower value
                if property == "km":
                    question += "lower."
                # However, improving kcat/km or kcat indicates a higher value    
                else:
                    question += "higher."
                objective = FUNCTIONS.answer_question(question, "string", [], \
                            ['improve', 'maintain', 'reduce', 'eliminate'])
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept = True
            # Store the Optzyme group         
            group = [property, objective, weight, assemblies]
            groups.append(group)

def PermittedKinds(permitted, experiment):
    """Identify Permitted Kinds of amino acids for a Design Position."""
    # Make a summary of the Molecules that have Design Positions
    answers = []
    summary = "\nMolecules with Design Positions"
    moleculeKeys = experiment["Design Positions"].keys()
    moleculeKeys.sort()
    for mn in moleculeKeys:
        if len(experiment["Design Positions"][mn]) > 0:
            answers.append(mn)
            summary += "\nMolecule " + mn + " contains " + \
                    str(len(experiment["Design Positions"][mn])) + " Residues"
    # Figure out which Molecule to mutate
    mn = None
    if len(answers) == 1:
        text = "\nMolecule " + answers[0] + " is the only Molecule that "
        text += "contains Design Positions, so it was automatically chosen."
        print screen_formatting(text)
        mn = answers[0]
    else:
        # Print the summary
        print screen_formatting(summary)
        question = "Which Molecule contains a Design Position you would "
        question += "like to specify permitted kinds of amino acids for?"
        mn = FUNCTIONS.answer_question(question, "string", [], answers)
    # Confirm that this is correct
    right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
    if not right:
        mn = None
    # Only proceed if there is a selected Molecule
    if mn != None:
        answers = experiment["Design Positions"][mn]
        # Make a list of those residues
        summary = "\nDesign Positions in Molecule " + mn + ":"
        summary += list_items(answers)
        # if there is only one choice, make it automatically
        if len(answers) == 1:
            text = "\nResidue " + answers[0] + " is the only Design Position in"
            text += " Molecule " + mn + ", so it was automatically selected."
            print screen_formatting(text)
            rn = answers[0]
        # otherwise ask
        else:
            print screen_formatting(summary)
            question = "For which Design Position would you like to specify "
            question += "permitted kinds of amino acids?"
            rn = FUNCTIONS.answer_question(question, "string", [], answers)
        # Confirm this is what the user wants to do
        right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
        if not right:
            mn = None
    # If the user has provided a specification, get the information
    if mn != None:
        # If there are already specifications for the Residue, get them
        if mn in permitted and rn in permitted[mn]:
            chosen = permitted[mn][rn]
        else:
            chosen = []
        # We're going to need the file format a lot to access various amino acid
        # related information from STANDARDS. Store it here
        ff = experiment["File Format"]
        # use a while loop
        accept = False
        pick = True
        while not accept:
            # Determine what amino acids are available for selection
            options = []
            for aa in aminoAcids[ff]:
                if ff == "PDB" and aa == "HSD":
                    pass
                elif aa in chosen:
                    pass
                else:
                    options.append(aa)
            # If it is possible to make a selection and that has been requested
            if pick:
                # Make a summary of the currently selected amino acids for the
                # Residue
                summary = "\nCurrently permitted kinds of amino acids for "
                summary += "Residue " + rn + " in Molecule " + mn + ":"
                summary += list_items(chosen)
                # If there's no choice, say no pick was made
                if len(options) == 0:
                    pick = False
                # If there's only one option, choose it automatically
                elif len(options) == 1:
                    print screen_formatting(summary)
                    text = "\n" + options[0] + " is the only remaining amino "
                    text += "acid, so it was automatically chosen."
                    aas = options
                # Otherwise, ask
                else:
                    print screen_formatting(summary)
                    # Assemble all possible answers
                    answers = list(options)
                    # Get the lower case and capitalized versions of those AAs,
                    # too, along with the one letter codes
                    for aa in options:
                        answers.extend([aa.lower(), aa.capitalize(), \
                                convertAA[ff][aa], convertAA[ff][aa].lower()])
                    # Finally, provide some additional options
                    answers.extend(['ALL', 'All', 'all'])
                    answers.extend(['NONE', 'None', 'none'])
                    # Start making the question
                    question = "What kind of amino acid is the Residue "
                    question += "permitted to mutate to? You may answer 'all' "
                    question += "to include all of them or 'none' if you would "
                    question += "not like to pick any."
                    # Go through the different groups
                    words = propertiesAA[ff].keys()
                    words.sort()
                    for word in words:
                        available = []
                        for aa in propertiesAA[ff][word]:
                            if aa in options:
                                available.append(aa)
                        if len(available) > 0:
                            question += "\nYou may answer '" + word.capitalize()
                            question += "' to use:" + list_items(available)
                            answers.extend([word.upper(), word, word.lower()])
                    # Get the answer
                    choice = FUNCTIONS.answer_question(question, "string", [], \
                                                       answers)
                    # Get the appropriate information from the answer
                    if choice.upper() == "NONE":
                        pick = False
                    elif choice.upper() == "ALL":
                        aas = options
                    elif choice.capitalize() in words:
                        aas = []
                        for aa in propertiesAA[ff][choice.capitalize()]:
                            if aa in options:
                                aas.append(aa)
                    elif len(choice) == 1:
                        aas = [convertAA[ff][choice.upper()]]
                    else:
                        aas = [choice.upper()]
            # if a pick was made, either automatically or manually, confirm 
            # that it is correct
            if pick:
                question = "Are you certain Residue " + rn + " should be "
                question += "allowed to mutate to "
                if len(aas) == 1:
                    question += "this amino acid?"
                else:
                    question += "these amino acids?"
                right = FUNCTIONS.answer_question(question, "bool")
                # If this is correct, store the information
                if right:
                    new = []
                    for aa in aminoAcids[ff]:
                        if ff == "PDB" and aa == 'HSD':
                            continue
                        if aa in aas or aa in chosen:
                            new.append(aa)
                    chosen = new
                else:
                    pick = False
            # If a choice was made, allow the option of making another
            if pick:
                question = "Would you like to specify more amino acids that "
                question += "Residue " + rn + " may mutate to?"
                more = FUNCTIONS.answer_question(question, "bool")
                if not more:
                    pick = False
            # If the user either can't or won't make another choice, give them
            # options
            if not pick:
                # Make a summary of the Residues
                summary = "\nCurrently permitted kinds of amino acids for "
                summary += "Residue " + rn + " in Molecule " + mn + ":"
                summary += list_items(chosen)
                print screen_formatting(summary)
                # Ask if the user is happy with this
                question = "Are you certain this list of amino acid kinds is "
                question += "correct?"
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    continue
                # If the user isn't happy, give them options
                answers = []
                phrases = []
                if len(options) > 0:
                    phrases.append("select 'more' permitted amino acids")
                    answers.extend(['MORE', 'More', 'more'])
                if len(chosen) > 0:
                    phrases.append("'start over' at selecting amino acids")
                    answers.extend(['START OVER', 'Start Over', 'Start over', \
                                    'start over'])
                if len(chosen) > 1:
                    phrases.append("'remove' a specific amino acid")
                    answers.extend(['REMOVE', 'Remove', 'remove'])
                # If there are no answers, break the while loop
                if len(phrases) == 0:
                    break
                # If there's only one, do that
                elif len(phrases) == 1:
                    what = answers[0]
                else:
                    question = "Would you like to " + phrases[0]
                    if len(phrases) == 2:
                        question += " or " + phrases[1]
                    else:
                        question += ", " + phrases[1]+", or "+phrases[2]
                    question += "?"
                    what = FUNCTIONS.answer_question(question, "string", [], \
                                                     answers)
                # Respond appropriately
                if what.upper() == "MORE":
                    pick = True
                    continue
                elif what.upper() == "START OVER":
                    pick = True
                    chosen = []
                    continue
                # Otherwise ask which one to remove
                answers = list(chosen)
                # Get the full list of possible AA names
                for aa in chosen:
                    answers.extend([aa.capitalize(), aa.lower(), \
                        convertAA[ff][aa], convertAA[ff][aa].lower()])
                # Ask the user
                question = "Which kind of amino acid should Residue " + rn 
                question += " NOT be permitted to change to?"
                aa = FUNCTIONS.answer_question(question, "string", [], answers)
                # Make sure it is a 3 letter code
                if len(aa) == 1:
                    aa = convertAA[ff][aa.upper()]
                # Get the index
                i = chosen.index(aa.upper())
                # Delete that aa from the list
                del chosen[i]
        # Now make sure that permitted is properly formatted
        if len(chosen) == 0:
            if mn in permitted and rn in permitted[mn]:
                del permitted[mn][rn]
            if permitted[mn] == {}:
                del permitted[mn]
        else:
            if mn not in permitted:
                permitted[mn] = {}
            permitted[mn][rn] = chosen
        # If there are Dimers, make certain the permitted kinds of amino acid
        # information matches
        if "Dimers" in experiment:
            for pair in experiment["Dimers"]:
                if mn in pair:
                    if mn == pair[0]:
                        mn2 = pair[1]
                    else:
                        mn2 = pair[0]
                    if mn in permitted:
                        permitted[mn2] = permitted[mn]
                    elif mn2 in permitted:
                        del permitted[mn2]
    # Nothing needs to be returned thanks to pass by reference

def Mutant(mutants, experiment):
    """Identify all mutations in a mutant"""
    # We're going to need file format information a lot (later on), so store
    # that in an easily written variable
    ff = experiment["File Format"]
    # Start by getting a list of the Design Molecules
    molecules = []
    for molecule in experiment["Molecules"]:
        if molecule[2].design:
            molecules.append(molecule[2])
    if len(molecules) == 0:
        text = "Somehow there are no Design Molecules, so no mutants may be "
        text += "found."
        raise FUNCTIONS.IPRO_IOError(text)
    # Store the information about which Residues may be mutated in each
    residues = {}
    for molecule in molecules:
        for residue in molecule:
            if residue.kind in aminoAcids[residue.fileFormat]:
                if molecule.name not in residues:
                    residues[molecule.name] = []
                residues[molecule.name].append(residue.name)
    # Store the mutant's information here
    mutant = []
    # Use a while loop
    accept = False
    pick = True
    while not accept:
        # If a choice should be made
        if pick:
            # Create a summary of the Molecules that may be chosen
            mns = []
            summary = ''
            for molecule in molecules:
                if molecule.name in residues:
                    mn = molecule.name
                    mns.append(mn)
                    summary += "\nMolecule " + mn + " contains "
                    summary += str(len(residues[mn])) + " Residues that may be "
                    summary += "selected for mutation"
            # if there are no more possible choices
            if mns == []:
                text = "\nThere are no more Residues eligible for mutation"
                print screen_formatting(text)
                pick = False
            # Otherwise, a choice must be made
            else:
                print screen_formatting(summary)
                # If there's only one option, pick it automatically
                if len(mns) == 1:
                    text = "\nMolecule " + mns[0] + " is the only possible "
                    text += "choice, so it was selected automatically."
                    print screen_formatting(text)
                    mn = mns[0]
                # Otherwise ask
                else:
                    question = "In which Molecule is the Residue that should be"
                    question += " mutated?"
                    mn = FUNCTIONS.answer_question(question, "string", [], mns)
        # If a Molecule has been found, get a Residue
        if pick:
            # Summarize the Molecule, excluding Residues that cannot mutate
            ignore = []
            molecule = None
            for mol in molecules:
                if mol.name == mn:
                    molecule = mol
                    break
            for residue in molecule:
                if residue.name not in residues[mn]:
                    ignore.append(residue.name)
            print screen_formatting(MOLECULES.summarize(molecule, ignore))
            # identify the Residue that should be mutated
            question = "Which Residue should be mutated?"
            rn = FUNCTIONS.answer_question(question, "string", [], residues[mn])
            # Confirm that the Residue should be mutated
            question = "Are you certain you want to mutate Residue " + rn
            question += " in Molecule " + mn + "?"
            pick = FUNCTIONS.answer_question(question, "bool")
        # If a mutation has been requested, get it
        if pick:
            question = "Residue " + rn + " in Molecule " + mn + " is currently "
            question += molecule[rn].kind + ". What should it be mutated to?"
            # Assemble all of the different answers
            options = []
            for aa in aminoAcids[ff]:
                if ff == "PDB" and aa == "HSD":
                    continue
                if aa != molecule[rn].kind:
                    options.append(aa)
            answers = list(options)
            for aa in options:
                answers.extend([aa.lower(), aa.capitalize(), \
                                convertAA[ff][aa], convertAA[ff][aa].lower()])
            # get the amino acid
            aa = FUNCTIONS.answer_question(question, "string", [],\
                                           answers).upper()
            # Convert it to a three letter code
            if len(aa) == 1:
                aa = convertAA[ff][aa]
            # Confirm it
            pick = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            # Store the mutation in the mutant
            if pick:
                mutant.append([mn, rn, aa])
                # Modify the residues dictionary
                i = residues[mn].index(rn)
                del residues[mn][i]
                if residues[mn] == []:
                    del residues[mn]
                # If there are Dimers in the Experiment, account for that, too
                if "Dimers" in experiment:
                    for pair in experiment["Dimers"]:
                        if mn in pair:
                            if mn == pair[0]:
                                mn2 = pair[1]
                            else:
                                mn2 = pair[0]
                            if rn in residues[mn2]:
                                mutant.append(residues[mn2, rn, aa])
                                i = residues[mn2].index(rn)
                                del residues[mn2][i]
                                if residues[mn2] == []:
                                    del residues[mn2]
                # Ask if the user wants to select another mutation
                question = "Would you like to select another mutation for this "
                question += "mutant?"
                pick = FUNCTIONS.answer_question(question, "bool")
                if pick:
                    continue
        # If a pick either shouldn't be made or can't be
        if not pick:
            # Summarize the mutants
            mutationNumbers = []
            summary = "\nMutations"
            for i, mutation in enumerate(mutant):
                summary += "\n" + str(i+1) + ") Residue " + mutation[1] + " in "
                summary += "Molecule " + mutation[0] + " to " + mutation[2]
                mutationNumbers.append(i+1)
            print screen_formatting(summary)
            # Find out if the user is happy with this
            question = "Are you certain these are the exact mutations for this "
            question += "mutant?"
            right = FUNCTIONS.answer_question(question, "bool")
            if right:
                accept = True
                continue
            # If the user isn't happy, give options
            phrases = []
            answers = []
            if residues != {}:
                phrases.append("select 'more' mutations")
                answers.extend(["MORE", "More", 'more'])
            if len(mutant) > 0:
               phrases.append("'start over' at picking this mutant's mutations")
               answers.extend(["START OVER", 'Start Over', 'Start over', \
                               'start over'])
            if len(mutant) > 1:
                phrases.append("'remove' a specific mutant")
                answers.extend(["REMOVE", "Remove", 'remove'])
            # If there are no options, break the while loop
            if len(phrases) == 0:
                break
            # If there's only one, do that
            if len(phrases) == 1:
                what = answers[0]
            else:
                question = "Would you like to " + phrases[0]
                if len(phrases) == 2:
                    question += " or " + phrases[1]
                else:
                    question += ", " + phrases[1] + ", or " + phrases[2]
                question += "?"
                what = FUNCTIONS.answer_question(question, "string", [], \
                                                 answers)
            # Respond appropriately
            if what.upper() == "MORE":
                pick = True
                continue
            if what.upper() == "START OVER":
                mutant = []
                pick = True
                # Reset the Residues dictionary
                residues = {}
                for molecule in molecules:
                    for residue in molecule:
                        if residue.kind in aminoAcids[ff]:
                            if molecule.name not in residues:
                                residues[molecule.name] = []
                            residues[molecule.name].append(residue.name)
                continue
            # Otherwise, ask which Mutant to delete
            question = "Which mutation should be removed from the mutant?"
            i = FUNCTIONS.answer_question(question, "integer", [], \
                                          mutationNumbers) - 1
            # Store that mutation in this list
            mutations = [mutant[i]]
            # Check for Dimer matching information
            if "Dimers" in experiment:
                for pair in experiment["Dimers"]:
                    if mutant[i][0] in pair:
                        if mutant[i][0] == pair[0]:
                            mn2 = pair[1]
                        else:
                            mn2 = pair[0]
                        # Search through the other mutations
                        for j in range(len(mutant)):
                            if j == i:
                                continue
                            if mutant[j][0] == mn2 and mutant[j][1] == \
                            mutant[i][1] and mutant[j] not in mutations:
                                mutations.append(mutant[j])
            # Fix residues as each mutation is removed
            for mutation in mutations:
                if mutation[0] not in residues:
                    residues[mutation[0]] = []
                residues[mutation[0]].append(mutation[1])
                i = mutant.index(mutation)
                del mutant[i]
    # If there is a mutant, make sure it is new
    if mutant != []:
        unique = True
        for other in mutants:
            if len(other) == len(mutant):
                same = True
                for mutation in mutant:
                    if mutation not in other:
                        same = False
                        break
                if same:
                    unique = False
                    break
        if unique:
            mutants.append(mutant)
        else:
            text = "\nThat mutant has been previously specified, so your "
            text += "selection was dicarded"
            print screen_formatting(text)
    # Thanks to pass by reference, nothing needs to be returned

def fixedAtoms(fixed, experiment):
    """Identify Atoms that are never allowed to move."""
    # This is the information that this function will try to find
    gn = None
    mn = None
    rn = None
    atoms = []
    # Identify the Molecule
    answers = ["NONE", 'None', 'none']
    summary = ''
    for data in experiment["Molecules"]:
        answers.append(data[2].name)
        summary += "\nMolecule " + data[2].name + " has " + str(len(data[2]))
        summary += " Residues"
    print screen_formatting(summary)
    question = "Which Molecule has Atoms that you want to fix in place? You "
    question += "may answer 'none' if this is a mistake."
    mn = FUNCTIONS.answer_question(question, "string", [], answers)
    if mn.upper() == "NONE":
        mn = None
    else:
        molecule = None
        for data in experiment["Molecules"]:
            if data[2].name == mn:
                molecule = data[2]
    # If a Molecule has been selected, get a Residue
    if mn != None:
        answers = ["NONE", 'None', 'none']
        for residue in molecule:
            answers.append(residue.name)
        # Summarize the MOLECULE
        print screen_formatting(MOLECULES.summarize(molecule))
        # Ask for the Residue
        question = "Which Residue has Atoms that should be fixed in place? You"
        question += " may answer 'none' if this is a mistake."
        rn = FUNCTIONS.answer_question(question, "string", [], answers)
        if rn.upper() == "NONE":
            rn = None
        else:
            residue = molecule[rn]
    # If a Residue was selected
    if rn != None:
        # Find the Atoms that are not permitted to be used
        exclude = []
        if mn in experiment["Design Positions"] and rn in \
        experiment["Design Positions"][mn]:
            for atom in residue:
                if atom.name not in backboneAtoms[residue.fileFormat]:
                    exclude.append(atom.name)
        # Use a while loop to get the atoms
        accept = False
        pick = True
        while not accept:
            # Determine what (other) Atoms may be fixed in place
            ignore = []
            options = []
            for atom in residue:
                if atom.name in exclude:
                    ignore.append(atom.name)
                elif atom.name in atoms:
                    ignore.append(atom.name)
                else:
                    options.append(atom.name)
            # If a choice should be made
            if pick:
                # If there are no options, say that
                if len(options) == 0:
                    text = "\nThere are no other Atoms in this Residue that "
                    text += "may be fixed in place"
                    print screen_formatting(text)
                    pick = False
                # If there's only one, pick it automatically
                elif len(options) == 1:
                    chosen = options
                    text = "\nAtom " + options[0] + " was the only remaining "
                    text += "choice, so it was picked automatically"
                    print screen_formatting(text)
                else:
                    # Print a summary
                    print screen_formatting(MOLECULES.summarize(residue,ignore))
                    # Ask which Atom
                    question = "Which Atom should be fixed in place? You may "
                    question += "answer 'all' to freeze all of them or 'none' "
                    question += "to select none of them."
                    answers = list(options)
                    answers.extend(["ALL", 'All', 'all'])
                    answers.extend(["NONE", 'None', 'none'])
                    choice = FUNCTIONS.answer_question(question, "string", [], \
                                                       answers)
                    # Get the appropriate atoms
                    if choice.upper() == "ALL":
                        chosen = options
                    elif choice.upper() == "NONE":
                        pick = False
                    else:
                        chosen = [choice]
            # If a selection was made
            if pick:
                new = []
                for atom in residue:
                    if atom.name in chosen or atom.name in atoms:
                        new.append(atom.name)
                atoms = new
                # Find out if more are wanted
                if len(options) > 1 and choice.upper() != "ALL":
                    question = "Would you like to freeze more Atoms in this "
                    question += "Residue in place?"
                    pick = FUNCTIONS.answer_question(question, "bool")
                else:
                    pick = False
            # If a choice either shouldn't be made or can't be
            if not pick:
                # Summarize the selections
                summary = "\nAtoms in Residue " + rn + " that may never move:"
                summary += list_items(atoms)
                print screen_formatting(summary)
                # Ask if the user is happy
                question = "Are you certain these are the Atoms in Residue "
                question += rn + " that may never move?"
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    continue
                # Otherwise, provide options
                phrases = []
                answers = []
                if len(options) > 0:
                    phrases.append("select 'more' Atoms")
                    answers.extend(['MORE', 'More', 'more'])
                if len(atoms) > 0:
                    phrases.append("'start over' at picking Atoms in Residue " \
                                   + rn)
                    answers.extend(["START OVER", 'Start Over', 'Start over', \
                                    "start over"])
                if len(atoms) > 1:
                    phrases.append("'remove' a specific Atom from the list")
                    answers.append(['REMOVE', 'Remove', 'remove'])
                # If there are no options, break the while loop
                if len(phrases) == 0:
                    break
                elif len(phrases) == 1:
                    what = answers[0]
                else:
                    question = "Would you like to " + phrases[0]
                    if len(phrases) == 2:
                        question += " or " + phrases[1]
                    else:
                        question += ", " + phrases[1] + ", or " + phrases[2]
                    question += "?"
                    what = FUNCTIONS.answer_question(question, "string", [], \
                                                     answers)
                # Respond appropriately
                if what.upper() == "MORE":
                    pick = True
                    continue
                elif what.upper() == "START OVER":
                    pick = True
                    atoms = []
                    continue
                question = "Which Atom should be removed from the list of "
                question += "frozen Atoms?"
                an = FUNCTIONS.answer_question(question, "string", [], atoms)
                i = atoms.index(an)
                del atoms[i]
    # If Atoms were selected to freeze in place, identify the Design Group
    if atoms != []:
        # Get the Design Groups this Molecule is in
        if molecule.design:
            numbers = range(1, len(experiment["Design Groups"]) + 1)
        else:
            numbers = []
            for i, group in enumerate(experiment["Design Groups"]):
                if molecule.name in group:
                    numbers.append(i+1)
        # If it is only a single group, use that
        if len(numbers) == 1:
            gn = 'all'
        # Otherwise, ask
        else:
            answers = []
            for i in numbers:
                answers.append(str(i))
            question = "Molecule " + mn + " appears in the following Design "
            question += "Groups:" + list_items(answers)
            question += ". In which Design Group should the Atoms be frozen? "
            question += "You may answer 'all' to use all of them."
            answers.extend(["ALL", "All", 'all'])
            gn = FUNCTIONS.answer_question(question, "string", [], answers)
            # Make sure the answer is correct
            if gn.upper() == "ALL":
                gn = 'all'
            else:
                gn = int(gn)
    # If a Group was selected, confirm that this should be done
    if gn != None:
        question = "Are you certain the specified Atoms in Residue " + rn+" in "
        question += "Molecule " + mn + " in "
        if gn == 'all':
            question += "all Design Groups"
        else:
            question += "Design Group " + str(gn)
        question += " should always be fixed in place?"
        right = FUNCTIONS.answer_question(question, "bool")
        # if so, store the information
        if right:
            if gn not in fixed:
                fixed[gn] = {}
            if mn not in fixed[gn]:
                fixed[gn][mn] = {}
            fixed[gn][mn][rn] = atoms
    # Fixed is passed by reference, so it doesn't need to be returned

def position_restraint(restraints, experiment):
    """Generate a position restraint"""
    # The parameters that will need to be found
    gn = None
    mn = None
    rn = None
    an = None
    gn2 = None
    # Force field related parameters
    if experiment["Force Field"] == "CHARMM":
        fc = None
    else:
        text = "The IO GET position restraint function does not support the "
        text += str(experiment["Force Field"]) + " force field."
        raise FUNCTIONS.IPRO_IOError(text)
    # Get the Molecule
    molecules = []
    for molecule in experiment["Molecules"]:
        molecules.append(molecule[2])
    # Use a while loop
    chosen = []
    accept = False
    pick = True
    while not accept:
        # Make a list of the possible options
        options = []
        for molecule in molecules:
            if molecule.name not in chosen:
                options.append(molecule.name)
        # If a Molecule should be picked
        if pick:
            # If there are no options, make no choice
            if len(options) == 0:
                text = "\nThere are no more Molecules to select"
                print screen_formatting(text)
                pick = False
            elif len(options) == 1:
                text = "\nMolecule " + options[0] + " was the only choice, so "
                text += "it was picked automatically"
                print screen_formatting(text)
                mols = options
            else:
                summary = ''
                for molecule in molecules:
                    if molecule.name in options:
                        summary += "\nMolecule " + molecule.name + " contains "
                        summary += str(len(molecule)) + " Residues"
                print screen_formatting(summary)
                # Ask which Molecule to use
                question = "Which Molecule would you like to select for the "
                question += "position restraint? You may answer 'all' to use "
                question += "all of them or 'none' if this was a mistake."
                answers = list(options)
                answers.extend(['ALL', 'All', 'all', "NONE", "None", 'none'])
                # Get the answer
                choice = FUNCTIONS.answer_question(question, "string", [], \
                                                   answers)
                # Respond appropriately
                if choice.upper() == "ALL":
                    mols = options
                elif choice.upper() == "NONE":
                    pick = False
                else:
                    mols = [choice]
        if pick:
            # Store the choice
            chosen.extend(mols)
            # If possible, ask if they want to include more
            if choice.upper() != 'ALL' and len(options) > 1:
                question = "Would you like to include more Molecules in the "
                question += "position restraint?"
                pick = FUNCTIONS.answer_question(question, "bool")
                if pick:
                    continue
            else:
                pick = False
        # If a choice wasn't or couldn't be made
        if not pick:
            # Summarize the Molecules
            summary = "\nMolecules of the position restraint:"
            summary += list_items(chosen)
            print screen_formatting(summary)
            # Confirm this is right
            right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            if right:
                accept = True
            else:
                text ="\nThe list of selected Molecules has been reset to empty"
                print screen_formatting(text)
                pick = True
                chosen = []
    # Store the chosen Molecule information properly
    if chosen == []:
        pass
    elif len(chosen) == 1:
        mn = chosen[0]
        molecule = None
        for mol in molecules:
            if mol.name == mn:
                molecule = mol
                break
    elif len(chosen) == len(molecules):
        mn = 'all'
        rn = 'all'
        an = 'all'
    else:
        mn = chosen
        rn = 'all'
        an = 'all'
    # If it is appropriate to make a Residue selection
    if mn not in [None, 'all'] and rn == None:
        # Store the chosen Residues in this list
        chosen = []
        # Use a while loop
        accept = False
        pick = True
        while not accept:
            # Make a list of the Residues that may be selected
            options = []
            for residue in molecule:
                if residue.name not in chosen:
                    options.append(residue.name)
            # If a choice should be made
            if pick:
                # If there are no more options
                if len(options) == 0:
                    text = "\nThere are no more Residues in Molecule " + mn
                    text += " to select."
                    print screen_formatting(text)
                    pick = False
                elif len(options) == 1:
                    reses = options
                    text = "\nResidue " + options[0] + " was automatically "
                    text += "selected."
                    print screen_formatting(text)
                # Otherwise ask,
                else:
                    print screen_formatting(MOLECULES.summarize(molecule, \
                                            chosen))
                    question = "Which Residue should be used in the position "
                    question += "restraint? You may answer 'all' to use them "
                    question += "all or 'none' for none."
                    answers = list(options)
                    answers.extend(['ALL', 'All', 'all', 'NONE', 'None','none'])
                    choice = FUNCTIONS.answer_question(question, "string", [], \
                                                       answers)
                    # Respond appropriately
                    if choice.upper() == "ALL":
                        reses = options
                    elif choice.upper() == "NONE":
                        pick = False
                    else:
                        reses = [choice]
            # If a choice was made, store it
            if pick:
                chosen.extend(reses)
                # If appropriate, ask about more residues
                if len(options) > 1 and choice.upper() != "ALL":
                    question = "Would you like to include more Residues in this"
                    question += " position restraint?"
                    pick = FUNCTIONS.answer_question(question, "bool")
                    if pick:
                        continue
                else:
                    pick = False
            # If a pick can't be made or isn't asked for
            if not pick:
                # Summarize the selected Residues
                summary = "\nResidues in Molecule " + mn + " in the position "
                summary += "restraint:" + list_items(chosen)
                print screen_formatting(summary)
                # Confirm this is correct
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept = True
                else:
                    text = "\nThe list of selected Residues has been reset to "
                    text += "empty."
                    print screen_formatting(text)
                    chosen = []
                    pick = True
        # If Residues were selected, store them appropriately
        if len(chosen) == 0:
            pass
        elif len(chosen) == 1:
            rn = chosen[0]
        else:
            rn = chosen
            an = 'all'
            if len(rn) == len(molecule):
                rn = 'all'
    # If specific Atoms can be selected
    if rn not in [None, 'all'] and an == None:
        # Get the Residue
        residue = molecule[rn]
        exclude = []
        # If it is a Design Position, exclude its side chain
        if mn in experiment["Design Positions"] and rn in \
        experiment["Design Positions"][mn]:
            for atom in residue:
                if atom.name not in backboneAtoms[residue.fileFormat]:
                    exclude.append(atom.name)
        # Use a while loop
        accept = False
        chosen = []
        pick = True
        while not accept:
            # Find the Atoms that may be selected, as well as those that should
            # be ignored
            options = []
            ignore = []
            for atom in residue:
                if atom.name in exclude:
                    ignore.append(atom.name)
                elif atom.name in chosen:
                    ignore.append(atom.name)
                else:
                    options.append(atom.name)
            # If a selection should be made
            if pick:
                if len(options) == 0:
                    text = "\nThere are no more Atoms in Residue " + rn
                    text += " that may be selected."
                    print screen_formatting(text)
                    pick = False
                elif len(options) == 1:
                    atoms = options
                    text = "\nAtom " + atoms[0] + " was chosen automatically"
                    print screen_formatting(text)
                else:
                    print screen_formatting(MOLECULES.summarize(residue,ignore))
                    question = "Which Atom should be used in the position "
                    question += "restraint? You may answer 'all' to use all of "
                    question += "them or 'none' for none."
                    answers = list(options)
                    answers.extend(['ALL', 'All', 'all', 'NONE', 'None','none'])
                    choice = FUNCTIONS.answer_question(question, "string", [], \
                                                       answers)
                    if choice.upper() == "NONE":
                        pick = False
                    elif choice.upper() == "ALL":
                        atoms = options
                    else:
                        atoms = [choice]
            # If a pick was made, store it
            if pick:
                chosen.extend(atoms)
                if len(options) > 1 and choice.upper() != "ALL":
                    question = "Would you like to include more Atoms in the "
                    question += "position restraint?"
                    pick = FUNCTIONS.answer_question(question, "bool")
                    if pick:
                        continue
                else:
                    pick = False
            # If a choice shouldn't be or can't be made
            if not pick:
                summary = "\nAtoms in the position restraint:"
                summary += str(chosen)
                print screen_formatting(summary)
                # Find out if the user is happy
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept = True
                else:
                    text ="\nThe list of selected Atoms has been reset to empty"
                    print screen_formatting(text)
                    pick = True
                    chosen = []
        # If Atoms were selected, store that information
        if len(chosen) == 0:
            pass
        elif len(chosen) == 1:
            an = chosen[0]
        elif len(chosen) == len(residue):
            an = 'all'
        else:
            an = chosen
    # If Atoms were selected, get the force field information
    if an != None:
        # Store the text for the final summary here
        fieldText = ''
        # Do this by force field
        if experiment["Force Field"] == "CHARMM":
            accept = False
            while not accept:
                question = "What force constant would you like to use for this "
                question += "position restraint?"
                fc = FUNCTIONS.answer_question(question, "float")
                if fc <= 0:
                    text = "\nYou must answer with a positive number"
                    print screen_formatting(text)
                    continue
                right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                if right:
                    accept = True
            # Create the field text
            fieldText += " with a force constant of " + format(fc, '.3f') 
            items = [fc]
    # Figure out the Design Group(s) to apply this to
    if an != None:
        # If there are multiple Molecules, apply this to all Design Groups
        if mn == 'all' or isinstance(mn, list):
            gn = 'all'
        # Otherwise, figure out which Design Groups are possible
        else:
            groups = []
            for i, group in enumerate(experiment['Design Groups']):
                if molecule.design or mn in group:
                    groups.append(str(i+1))
            # If there's only one option, use 'all'
            if len(groups) == 1:
                gn = 'all'
            # Otherwise, ask
            else:
                summary = "\nMolecule " + mn + " appears in the following "
                summary += "Design Groups:" + list_items(groups)
                print screen_formatting(summary)
                groups.extend(['ALL', 'All', 'all'])
                # Use a while loop
                accept = False
                while not accept:
                    question = "Which Design Group should the position "
                    question += "restraint be applied to? You may answer 'all' "
                    question += "to use it on all of them."
                    gn = FUNCTIONS.answer_question(question, "string", [], \
                                                   groups)
                    if gn.upper() == "ALL":
                        gn = 'all'
                    else:
                        gn = int(gn)
                    right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
                    if right:
                        accept = True
    # Give the user a final option of rejecting this restraint
    if gn != None:
        question = "Are you certain the position restraint you've specified "
        question += "is correct?"
        right = FUNCTIONS.answer_question(question, "bool")
        if right:
            # Make the restraint
            restraint = [gn, mn, rn, an]
            restraint.extend(items)
            # We're not actually providing the option of restraining to another
            # design group
            restraints.append(restraint)

def distance_restraint(restraints, experiment):
    """Get a distance restraint"""
    # Make sure the Experiment's force field will be supported later in the
    # function
    if experiment["Force Field"] not in ['CHARMM']:
        text = "The IO GET distance restraint function does not support the "
        text += str(experiment["Force Field"]) + " force field."
        raise FUNCTIONS.IPRO_IOError(text)
    # Get the two Atoms that will be used in the restraint
    atoms = []
    while len(atoms) < 2:
        # Get the Molecules
        molecules = {}
        answers = []
        for data in experiment["Molecules"]:
            molecules[data[2].name] = data[2]
            answers.append(data[2].name)
        # If there are no Molecules, throw an error
        if len(answers) == 0:
            text = "Somehow there are no Molecules"
            raise FUNCTIONS.IPRO_IOError(text)
        # If there's only one Molecule, choose it automatically
        elif len(answers) == 1:
            text = "\nMolecule " + answers[0] + " was automatically chosen "
            text += "for the distance restraint."
            mn = answers[0]
        else:
            # Make a summary of the Molecules
            summary = ''
            for mn in answers:
                summary += "\nMolecule " + mn + " contains " + \
                           str(len(molecules[mn])) + " Residues"
            print screen_formatting(summary)
            # Ask which Molecule to use
            question = "Which Molecule contains Atom " + str(len(atoms) + 1)
            question += " for the distance restraint?"
            mn = FUNCTIONS.answer_question(question, "string", [], answers)
        # Get the Residue
        molecule = molecules[mn]
        print screen_formatting(MOLECULES.summarize(molecule))
        if len(molecule) == 1:
            rn = molecule[0].name
            text = "\nResidue " + rn + " was automatically selected for the "
            text += 'distance restraint'
            print screen_formatting(text)
        else:
            question = "Which Residue contains Atom " + str(len(atoms) + 1)
            question += " for the distance restraint? You may answer 'none' "
            question += "if you did not mean to select Molecule " + mn
            answers = ["NONE", 'None', 'none']
            for residue in molecule:
                answers.append(residue.name)
            rn = FUNCTIONS.answer_question(question, "string", [], answers)
            if rn.upper() == "NONE":
                continue
        # Get the Atom
        residue = molecule[rn]
        # If this Residue is a Design Position, don't do things with its side
        # chain
        exclude = []
        if mn in experiment["Design Positions"] and rn in \
        experiment["Design Positions"][mn]:
            for atom in residue:
                if atom.name not in backboneAtoms[residue.fileFormat]:
                    exclude.append(atom.name)
        answers = ['NONE', 'None', 'none']
        options = []
        for atom in residue:
            if atom.name not in exclude:
                options.append(atom.name)
        if len(options) == 1:
            text = "\nAtom " + options[0] + " was chosen automatically"
            an = options[0]
            print screen_formatting(text)
        else:
            question = "Which Atom should be used in the distance restraint? "
            question += "You may answer 'none' if you did not mean to select "
            question += "Residue " + rn + " in Molecule " + mn 
            answers.extend(options)
            print screen_formatting(MOLECULES.summarize(residue, exclude))
            an = FUNCTIONS.answer_question(question, "string", [], answers)
            if an.upper() == "NONE":
                continue
        # Confirm that the selection is unique
        if [mn, rn, an] in atoms:
            text = "\nThat atom was already specified. Please try again."
            print screen_formatting(text)
        # Confirm the user wants to use this Atom
        question = "Are you certain you want to use Atom " + an + " in Residue "
        question += rn + " in Molecule " + mn + " in the distance restraint?"
        right = FUNCTIONS.answer_question(question, "bool")
        if right:
            atoms.append([mn, rn, an])
    # Figure out what Design Group(s) has both those Molecules
    # Start by figuring out if either / both Molecules is a design molecule
    d1 = False
    d2 = False
    for mn in molecules:
        if molecules[mn].design:
            if atoms[0][0] == mn:
                d1 = True
            if atoms[1][0] == mn:
                d2 = True
    groups = []
    for i, group in enumerate(experiment["Design Groups"]):
        if (d1 or atoms[0][0] in group) and (d2 or atoms[1][0] in group):
            groups.append(str(i+1))
    # If there are no groups, tell the user
    if len(groups) == 0:
        text = "\nThere are no Design Groups that contain both Atoms"
        print screen_formatting(text)
        gn = None
    elif len(groups) == 1:
        gn = 'all'
    else:
        summary = "\nDesign Groups containing both Atoms:" + list_items(groups)
        print screen_formatting(summary)
        accept = False
        while not accept:
            question = "Which Design Group should the distance restraint be "
            question += "used in? You may answer 'all' for all of them."
            answers = list(groups)
            answers.append(['ALL', 'All', 'all'])
            gn = FUNCTIONS.answer_question(question, "string", [], answers)
            if gn.upper() == "ALL":
                gn = 'all'
            else:
                gn = int(gn)
            right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            if right:
                accept = True
    # If there is a selected Design Group, get the force field parameters
    if gn != None:
        # Get the two Atoms, and calculate the distance between them
        atom1 = molecules[atoms[0][0]][atoms[0][1]][atoms[0][2]]
        atom2 = molecules[atoms[1][0]][atoms[1][1]][atoms[1][2]]
        dis = MOLECULES.calculate_distance(atom1, atom2)
        # CHARMM
        if experiment["Force Field"] == "CHARMM":
            # use a while loop to get all of the needed items
            accept = False
            while not accept:
                text = "\nThe two Atoms are currently " + format(dis, '.3f') 
                text += " Angstroms apart."
                print screen_formatting(text)
                question = "What is the minimum permissible distance between "
                question += "the Atoms?"
                accept2 = False
                while not accept2:
                    rmin = FUNCTIONS.answer_question(question, "float")
                    if rmin <= 0:
                        print "\nYou must answer with a positive number"
                    else:
                        accept2 = True
                question = "What force constant should be used if this minimum "
                question += "distance is violated?"
                accept2 = False
                while not accept2:
                    kmin = FUNCTIONS.answer_question(question, "float")
                    if kmin <= 0:
                        print "\nYou must answer with a positive number"
                    else:
                        accept2 = True
                question = "What is the maximum permissible distance between "
                question += "the Atoms?"
                accept2 = False
                while not accept2:
                    rmax = FUNCTIONS.answer_question(question, "float")
                    if rmax < rmin:
                        text = "\nYou must answer with a number greater than "
                        text += "the minimum distance"
                        print screen_formatting(text)
                    else:
                        accept2 = True
                question = "What force constant should be used if the maximum "
                question += "distance is violated?"
                accept2 = False
                while not accept2:
                    kmax = FUNCTIONS.answer_question(question, "float")
                    if kmax <= 0:
                        print "\nYou must answer with a positive number"
                    else:
                        accept2 = True
                question = "What is the maximum force the distance restraint "
                question += "may generate?"
                accept2 = False
                while not accept2:
                    fmax = FUNCTIONS.answer_question(question, "float")
                    if fmax <= 0:
                        print "\nYou must answer with a positive number"
                    else:
                        accept2 = True
                # Confirm that these values are correct
                question = "Are you certain you want to use the following "
                question += "settings:\nKmin: " + format(kmin, '.3f')
                question += "\nRmin: " + format(rmin, '.3f')
                question += "\nKmax: " + format(kmax, '.3f')
                question += "\nRmax: " + format(rmax, '.3f')
                question += "\nFmax: " + format(fmax, '.3f')
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    items = [kmin, rmin, kmax, rmax, fmax]
        # Provide one last opportunity to be certain the restraint is correct
        question = "Are you certain the distance restraint you have specified "
        question += "is correct?"
        right = FUNCTIONS.answer_question(question, "bool")
        if right:
            restraint = [gn, atoms[0], atoms[1]]
            restraint.extend(items)
            restraints.append(restraint)

def dihedral_restraint(restraints, experiment):
    """Generate a dihedral restraint"""
    # Make sure the force field will be supported
    if experiment["Force Field"] not in ["CHARMM"]:
        text = "The IO GET dihedral restraint function does not support the "
        text += str(experiment["Force Field"]) + " force field."
        raise FUNCTIONS.IPRO_IOError(text)
    # Make a dictionary of all of the Molecules and have a list of their names
    # in order
    molNames = []
    molecules = {}
    for data in experiment["Molecules"]:
        molecules[data[2].name] = data[2]
        molNames.append(data[2].name)
    if len(molNames) == 0:
        text = "Somehow there are no Molecules"
        raise FUNCTIONS.IPRO_IOError(text)
    # Store the collected Atoms in this list
    atoms = []
    while len(atoms) < 4:
        # Get the Molecule
        if len(molNames) == 1:
            text = "\nMolecule " + molNames[0] + " was automatically selected"
            print screen_formatting(text)
            mn = molNames[0]
        # Ask the user
        else:
            summary = ''
            for mn in molNames:
                summary += "\nMolecule " + mn + " contains " + \
                           str(len(molecules[mn])) + " Residues"
            print screen_formatting(summary)
            question = "Which Molecule contains Atom " + str(len(atoms) + 1)
            question += " for the dihedral restraint?"
            mn = FUNCTIONS.answer_question(question, "string", [], molNames)
        # Get that Molecule
        molecule = molecules[mn]
        # Get the Residue
        if len(molecule) == 1:
            rn = molecule[0].name
            text = "\nResidue " + rn + " was automatically selected"
            print screen_formatting(text)
        else:
            print screen_formatting(MOLECULES.summarize(molecule))
            question = "Which Residue contains Atom " + str(len(atoms) + 1)
            question += " for the dihedral restraint?"
            answers = []
            for residue in molecule:
                answers.append(residue.name)
            rn = FUNCTIONS.answer_question(question, "string", [], answers)
        residue = molecule[rn]
        # And get the Atom
        if len(residue) == 1:
            an = residue[0].name
            text = "\nAtom " + an + " was automatically selected"
            print screen_formatting(text)
        else:
            # If the Residue is a Design Position, don't consider side chain
            # options
            exclude = []
            if mn in experiment["Design Positions"] and rn in \
            experiment["Design Positions"][mn]:
                for atom in residue:
                    if atom.name not in backboneAtoms[residue.fileFormat]:
                        exclude.append(atom.name)
            options = []
            for atom in residue:
                if atom.name not in exclude:
                    options.append(atom.name)
            print screen_formatting(MOLECULES.summarize(residue, exclude))
            # Get the answer
            an = FUNCTIONS.answer_question(question, "string", [], options)
        # If the Atom has already been used
        atom = [mn, rn, an]
        if atom in atoms:
            text = "\nThat Atom is already specified in the restraint. Please "
            text += "try again."
            print screen_formatting(text)
        else:
            question = "Are you certain you want to use Atom " + an + " in "
            question += "Residue " + rn + " in Molecule " + mn + " in the "
            question += "dihedral restraint?"
            right = FUNCTIONS.answer_question(question, "bool")
            if right:
                atoms.append(atom)
    # Get the Design Groups that are possible
    groups = []
    for i in range(len(experiment["Design Groups"])):
        groups.append(i + 1)
    # Test each of the Molecules to see what groups they are in
    for atom in atoms:
        # Design Molecules are in every Design Group
        if molecules[atom[0]].design:
            continue
        for i in range(len(experiment["Design Groups"])):
            if atom[0] not in experiment["Design Groups"][i]:
                if i+1 in groups:
                    j = groups.index(i+1)
                    del groups[j]
    # If there are no groups, tell the user
    if len(groups) == 0:
        text = "\nThere are no Design Groups that all 4 Atoms are in "
        text += "simultaneously"
        print screen_formatting(text)
        gn = None
    elif len(groups) == 1:
        gn = 'all'
    else:
        answers = []
        for i in groups:
            answers.append(str(i))
        summary = "\nPossible Design Groups:" + list_items(answers)
        print screen_formatting(summary)
        question = "Which Design Group would you like to apply the dihedral "
        question += "restraint to? You may answer 'all' to use all of them"
        answers.extend(["ALL", "All", 'all'])
        accept = False
        while not accept:
            gn = FUNCTIONS.answer_question(question, "string", [], answers)
            if gn.upper() == "ALL":
                gn = 'all'
            else:
                gn = int(gn)
            right = FUNCTIONS.answer_question(FUNCTIONS.confirm, "bool")
            if right:
                accept = True
    # If there's a Design Group, get the other values
    if gn != None:
        # CHARMM
        if experiment["Force Field"] == "CHARMM":
            accept = False
            while not accept:
                question = "What should the force constant for the dihedral "
                question += "restraint be?"
                accept2 = False
                while not accept2:
                    fc = FUNCTIONS.answer_question(question, "float")
                    if fc <= 0:
                        print "\nYou must answer with a positive number"
                    else:
                        accept2 = True
                question = "In degrees, to what angle should the dihedral "
                question += "angle be minimized?"
                accept2 = False
                while not accept2:
                    angle = FUNCTIONS.answer_question(question, "float")
                    if not -180 <= angle <= 180:
                        text = "\nThe dihedral angle must be between -180 and "
                        text += "180."
                        print screen_formatting(text)
                    else:
                        accept2 = True
                # Confirm these values are correct
                question = "Are you certain you want to use a force constant of"
                question += " " + format(fc, '.3f') + " to minimize the "
                question += "dihedral angle " + format(angle, '.3f') 
                question += " degrees?"
                right = FUNCTIONS.answer_question(question, "bool")
                if right:
                    accept = True
                    items = [fc, angle]
        # Make extra certain this is what the user wants to do
        question = "Are you certain you have specified the proper dihedral "
        question += "restraint?"
        right = FUNCTIONS.answer_question(question, "bool")
        if right:
            restraint = [gn, atoms[0], atoms[1], atoms[2], atoms[3]]
            restraint.extend(items)
            restraints.append(restraint)

