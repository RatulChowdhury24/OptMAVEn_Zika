#!/usr/bin/env python

# Name of this file
__name__ = "The IPRO Suite Submitter Module"
# Documentation string
__doc__ = """
Originally written in 2013 by Robert Pantazes of the Costas Maranas Lab in the
Chemical Engineering Department of the Pennsylvania State University.

This file contains the functions needed to generate and execute a script that
can be run on a computing cluster (i.e. not manually). It needs to be edited to
match the requirements of every new group to install the IPRO Suite."""

# Include PYTHON Modules
import os
import sys
# Include everything from the STANDARDS module, which in this case is useful for
# the functions and parameters used in answering questions
from STANDARDS import *
#mfa from IO_FUNCTIONS import answer_question
#mfa from IO_FUNCTIONS import confirm

def submit_script(fileName):
    """A function that submits a script to run on a cluster of machines"""
    print "submitting ", fileName
    i = os.system("qsub " + fileName)

def script_generator(experiment):
    """Create one or more scripts to run the experiment."""
    # First, find out if the user would like to generate any scripts
    question = "Would you like to generate scripts for the automatic running "
    question += "of the " + experiment["Name"] + " " + experiment["Type"]
    question += " experiment?"
    # This question returns a boolean true or false value
    do = answer_question(question, "bool")
    # If the user doesn't want scripts made, don't do anything
    if not do:
        return do
    # Find out how many processors to have run the experiment
    question = "Currently, the IPRO Suite cannot run truly in parellel, due to "
    question += "incompatibilities of various programs. However, it is written "
    question += "such that multiple, independent processors can work on the "
    question += "experiment at the same time. How many processors would you "
    question += "like to have working on this experiment?"
    # Use a short while loop to make sure they give an appropriate answer
    accept = False
    while not accept:
        # This time the question will return an integer
        N = answer_question(question, "int")
        # If this integer is less than 1, ask again
        if N < 1:
            comment = "The number of processors must be a positive integer"
            print comment
            continue
        # Confirm that this is right. confirm is a built in question in the
        # STANDARDS module that just asks "Are you sure this is right?"
        right = answer_question(confirm, "bool")
        if right:
            accept = True
    # Find out the maximum length of time the job may run for. For the CDM group
    # this is 240 hours
    question = "What is the maximum number of hours a job may run for?"
    question += "\nFor the CDM group, the answer must be between 1 and 240."
    # Here we are expecting an integer answer between 1 and 240. We don't want
    # those possibilities output (there's a lot) if the answer isn't acceptable,
    # so put it in the second set of possible answers
    hours = answer_question(question, "int", [], range(1, 241))
    # Now comes parts where things will likely have to be edited. For the CDM
    # group at Penn State, we may be running these jobs on either Lion-XJ or
    # Lion-XF. Ask which cluster the user is on.
    question = "Which cluster are you using, 'J' or 'F'?"
    # Here we are expecting a string answer, and because we know the only
    # possibilities are J or F we can explicitly state those so the answer
    # question function will make sure one of those answers is used.
    C = answer_question(question, "string", ['J', 'F'])
    # All of the CDM scripts will have the same content, so only create that
    # once. However, the starting information is dependent on what cluster the
    # job is submitted to
    text = ''
    if C == 'F':
        text += "#PBS -l pmem=4gb\n"
    else:
        text += "#PBS -l nodes=1:ppn=1\n"
    # Make the rest of the script
    text += "#PBS -l walltime=" + str(hours) + ":00:00\n"
    text += "#PBS -j oe\n"
    text += "#PBS -o /dev/null\n"
    text += "#PBS -q lionx" + C.lower() + "-cdm\n\n"
    text += "set -u\n"
    text += "cd " + experiment["Folder"] + "\n\n"
    text += 'echo " "\n'
    text += 'echo " "\n'
    text += 'echo "Job started on $(hostname -s) at $(date)"\n'
    text += "python " + experiment["Type"] + ".py\n"
    text += 'echo " "\n'
    text += 'echo "Job ended at $(date)"\n'
    text += 'echo " "\n'
    # Write this to the proper number of files
    for I in range(1, N+1):
        f = open("job_" + experiment["Name"] + "_" + str(I), "w")
        f.write(text)
        f.close()
        # Make the file executable
        os.system("chmod u+x job_" + experiment["Name"] + "_" + str(I))
    # Find out if the user would like these jobs submitted
    question ="Would you like the calculations to be automatically started now?"
    do = answer_question(question, "bool")
    if do:
        for I in range(1, N+1):
            submit_script("job_" + experiment["Name"] + "_" + str(I))

#mfa#
def experiment_script(experiment, name, command, directory=None, hours=23, time_=None, self_destruct=False):
    if time_ is None:
        qscript_name = name
        qscript_command = command
    else:
        qscript_name = os.path.splitext(name)[0] + ".timer.sh"
        slavescript_name = name
        time_file = os.path.splitext(name)[0] + ".time.txt"
        qscript_command = "(time -p {}) 2> {}".format(slavescript_name, time_file)
        qscript_command += "\ncd {}".format(experiment["Folder"])
        qscript_command += "\npython {}/modules/PERFORMANCE.py {} {}".format(
                InstallFolder, time_file, time_)
        if self_destruct:
            qscript_command += "\nrm {} {}".format(slavescript_name, time_file)
        open(slavescript_name, "w").write(command)
        os.system("chmod u+x {}".format(slavescript_name))
    if self_destruct:
        qscript_command += "\nrm " + qscript_name
    text = "#PBS -l pmem=4gb\n"
    text += "#PBS -l walltime=" + str(hours) + ":00:00\n"
    text += "#PBS -j oe\n"
    text += "#PBS -o /dev/null\n"
    text += "#PBS -q lionxf\n\n"
    text += "set -u\n"
    text += 'echo "Job started on $(hostname -s) at $(date)"\n'
    if directory is not None:
        text += "cd {}\n".format(directory)
    text += qscript_command + "\n"
    text += 'echo "Job ended at $(date)"\n'
    open(qscript_name, "w").write(text)
    os.system("chmod u+x {}".format(qscript_name))
    submit_script(qscript_name)
