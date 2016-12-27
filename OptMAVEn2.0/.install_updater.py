#!/usr/bin/env python

# This program updates all programs and the STANDARDS module with the current
# folder of the installation. That way the programs and functions all know where
# to look to access information.

import os
# Get the current folder
current = os.getcwd() + "/"
# Update the appropriate files
for fileName in ['programs/IPRO.py', 'modules/STANDARDS.py', \
                 'programs/Start_Experiment.py', 'programs/Mutator.py']:
    # Store the contents of the file
    lines = []
    f = open(fileName, "r")
    for line in f:
        # If it is the line that needs to be replaced, do so
        if line.startswith("InstallFolder"):
            lines.append('InstallFolder = "' + current + '"\n')
        # Otherwise just store the line
        else:
            lines.append(line)
    f.close()
    # Rewrite the file with the updated text
    f = open(fileName, "w")
    for line in lines:
        f.write(line)
    f.close()
