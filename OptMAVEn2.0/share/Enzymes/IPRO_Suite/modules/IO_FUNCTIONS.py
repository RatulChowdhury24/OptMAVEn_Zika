#!/usr/bin/env python

# The name of this module
__name__ = "IPRO Suite Basic Input / Output Functions"
# Documentation
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University.

This file contains functions for the input and output of information for setting
up IPRO Suite Experiments and for running them."""

# Include needed PYTHON modules
import os
import sys
# Include the CHECK module of the Input / Output folder
import IO_CHECK as CHECK
# Include the contents of the STANDARDS module, which has standard IPRO settings
# and information
from STANDARDS import *
# Also include the EXPERIMENT and MOLECULES modules
import EXPERIMENT
import MOLECULES

# Some generic information for answering questions
yesNo = ['YES', 'Y', 'NO', 'N', 'yes', 'y', 'no', 'n', 'Yes', 'No', "TRUE", \
"True", "true", "t", "T", "FALSE", "False", "false", "f", "F"]
yes = ['YES', 'Yes', 'yes', 'Y', 'y', "TRUE", "True", "true", "t", "T"]
confirm = "Are you sure this is correct?"

class IPRO_IOError(IPRO_Error):
    """An error class for Input / Output problems in IPRO."""
    def __init__(self, error = ''):
        """The initialization of the IPRO I/O Error class."""
        IPRO_Error.__init__(self, error)

def answer_question(question, kind, printAnswers = [], permitAnswers = [],
                    ignoreAnswers = []):
    """Get user input in answer to a question."""
    # If there is only a single permitted answer, just return that
    if len(printAnswers) + len(permitAnswers) == 1:
        if len(printAnswers) == 1:
            answer = printAnswers[0]
        else:
            answer = permitAnswers[0]
        text = str(answer) + " was the only permitted answer, so it was chosen "
        text += "automatically."
        print screen_formatting("\n" + text)
        return answer
    # Otherwise ask the user
    # Modify the question so that it prints to the screen properly
    question = "\n" + screen_formatting(question)
    # Keep track of whether or not the answer is allowed
    accept = False
    while not accept:
        # Ask the question and get the answer
        print question
        answer = raw_input()
        # Validate the the type of answer is acceptable
        error = False
        if kind.lower() in ['int', 'integer']:
            # Try to make an integer, which won't work if the string has a
            # decimal place in it (i.e. is a float)
            try:
                answer = int(answer)
            except ValueError:
                error = True
        # If it is a floating point number, try to make that
        elif kind.lower() in ['float', 'floating']:
            try:
                answer = float(answer)
            except ValueError:
                error = True
        # If it should be a boolean answer, the answer should have been yes or
        # no
        elif kind.lower() in ['bool', 'boolean']:
            if answer not in yesNo:
                error = True
            elif answer in yes:
                # Just return a value, as boolean questions don't do the answer
                # checking
                return True
            else:
                return False
        # If there's a problem, print an error message
        if error:
            text = "\nYou must answer with a "+kind+" value. Please try again."
            print text
            continue
        # If the answer has been explicitly stated as something to ignore
        elif answer in ignoreAnswers:
            text = "\nThat answer is not permitted at this time. Please try "
            text += "again."
            print text
            continue
        # If there are no lists of the answers that are allowed, accept this
        # answer. Also do so if the answer is in one of the allowed lists
        if (len(printAnswers) == 0 and len(permitAnswers) == 0) or answer in \
        printAnswers or answer in permitAnswers:
            accept = True
        # If the answer isn't allowed, create an error message
        else:
            text = "\nThat answer is not permitted. Please try again."
            if len(printAnswers) > 0:
                text += " Possible answers are:"
                for term in printAnswers:
                    text += " " + term
            print screen_formatting(text)
    return answer

def get_value(attribute, kind, experiment = None, answers = [], default=False, \
              defaultValue = None, question = None):
    """Retrieve a value to use for a particular attribute in an Experiment."""
    # If the default value should be used, try to do so
    if default:
        # If this is an Experiment input, try to store that value (dictionaries
        # also OK)
        if isinstance(experiment, (EXPERIMENT.Experiment, dict)):
            try:
                experiment[attribute] = defaultValue
            except EXPERIMENT.ExperimentError as error:
                text = "\nUsing the default " + attribute+" has failed because:"
                text += str(error)
                print screen_formatting(text)
                default = False
        # If it's not an experiment, just return the value
        else:
            return defaultValue
    # If user input is needed, ask for it
    if not default:
        # Use a while loop
        accept = False
        while not accept:
            # If the question isn't a string, generate a default question
            if not isinstance(question, str):
                question = "What should the " + attribute + " value be?"
                if defaultValue != None:
                    question += " The default value is " + str(defaultValue)+"."
            # Get the answer to the question
            answer = answer_question(question, kind, answers)
            # Try to store that answer to make sure it is useable
            if isinstance(experiment, (EXPERIMENT.Experiment, dict)):
                try:
                    experiment[attribute] = answer
                except IPRO_Error as error:
                    print str(error)
                    continue
            # Confirm that this value is correct
            accept = answer_question(confirm, "bool")
        # If the answer needs to be returned
        if not isinstance(experiment, (EXPERIMENT.Experiment, dict)):
            return answer

def modify_list(items, single, plural):
    """Determine how to modify a set of selected items"""
    # If there are no items, just return the empty list and an indication that
    # values should be asked for
    if len(items) == 0:
        return items, True
    # Start assembling a question
    question = "Would you like to pick 'more' " + plural
    # Store answers here
    answers = ['MORE', 'More', 'more']
    # If there is more than one item, ask if the user wants to remove a specific
    # one
    if len(items) > 1:
        question += ", 'remove' a specific " + single + ","
        answers.extend(['REMOVE', 'Remove', 'remove'])
    question += " or 'start over' at picking them?"
    answers.extend(['START OVER', 'Start Over', 'Start over', 'start over'])
    # Get the answer
    answer = answer_question(question, "string", [], answers)
    # If the user wants more, be done
    if answer.upper() == "MORE":
        return items, True
    # if they want to start over, return an empty list 
    elif answer.upper() == 'START OVER':
        return [], True
    # Otherwise, ask which value to delete
    if isinstance(items[0], str):
        answers = items
        kind = "string"
    else:
        answers = range(1, len(items))
        kind = "integer"
    # Figure out which value to remove
    question = "Which " + single + " should be removed?"
    n = answer_question(question, kind, [], answers)
    # Get the index of the answer
    if isinstance(n, int):
        n -= 1
    else:
        n = items.index(n)
    del items[n]
    # Return the modified items as well as an indication that more values should
    # not be asked for.
    return items, False

def get_list(attribute, experiment = None, answers = [], default = False, \
             defaultValue = None):
    """Retrieve a LIST of values for an IPRO Suite Experiment attribute."""
    # If the default value should be used, do so
    if default:
        # If the default value can be stored, do so
        if isinstance(experiment, (EXPERIMENT.Experiment, dict)):
            try:
                experiment[attribute] = defaultValue
            except EXPERIMENT.ExperimentError as error:
                text = "\nUsing the default " + attribute + " has failed "
                text += "because:" + str(error)
                print screen_formatting(text)
                default = False
        # Otherwise, just return the value
        else:
            return defaultValue
    # If user input is needed to get the value
    if not default:
        # Store the values in this list
        values = []
        # Use a while loop
        accept = False
        # And have it set up to ask the user for a value
        ask = True
        while not accept:
            # If a value should be asked for
            if ask:
                # All of the attributes that use this function end with 's'. Cut
                # that off when asking the question
                question = "What is a " + attribute[:-1] + "? The " + attribute
                question += " specified so far are:" + list_items(values)
                # Determine if there are default values remaining
                allowedDefaults = []
                if isinstance(defaultValue, list):
                    for value in defaultValue:
                        if value not in values:
                            allowedDefaults.append(value)
                # Make a list of the acceptable answers right now
                options = list(answers)
                # If there are default values that may be used
                if len(allowedDefaults) > 0:
                    question += "\nDefault " + attribute + " that have not been"
                    question += " used yet are:" + list_items(allowedDefaults)
                if len(allowedDefaults) > 1:
                    question += "\nYou may answer 'all' to use all of these "
                    question += "default values."
                    # If the options list is not empty, modify it to allow all
                    # answers
                    if len(options) > 0:
                        options.extend(['ALL', 'All', 'all'])
                # Get an answer, making sure a value isn't duplicated
                answer = answer_question(question, "string", options, [], values)
                # If the answer is 'all', use all remaining default values
                if len(allowedDefaults) > 1 and answer.upper() == "ALL":
                    values.extend(allowedDefaults)
                # Otherwise check the value
                else:
                    try:
                        CHECK.list_check(attribute, [answer])
                        values.append(answer)
                    except IPRO_IOError as error:
                        print str(error)
            # If nothing has been selected, continue the while loop
            if len(values) == 0:
                ask = True
                continue
            # Summarize the selected items
            summary = "\nSelected " + attribute + ":" + list_items(values)
            print screen_formatting(summary)
            # Determine if another value should be gotten
            question = "Would you like to specify another " + attribute[:-1]+"?"
            another = answer_question(question, "bool")
            if another:
                ask = True
                continue
            # If the user does not want another value, try to store the values
            right = True
            if isinstance(experiment, (EXPERIMENT.Experiment, dict)):
                try:
                    experiment[attribute] = values
                except EXPERIMENT.ExperimentError as error:
                    print str(error)
                    right = False
            # Find out if the user is happy with the values
            if right:
                question = "Are you certain you have selected the correct " 
                question += attribute + "?"
                right = answer_question(question, "bool")
            # If the results are correct and were easily stored, be done
            if right:
                accept = True
                continue
            # If everything is not correct, allow the user the option of
            # modifying the values
            values, ask = modify_list(values, attribute, attribute[:-1])
        # If the answer needs to be returned
        if not isinstance(experiment, (EXPERIMENT.Experiment, dict)):
            return values
