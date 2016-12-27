""" This file controls submitting jobs to the queue.
"""

import os
from subprocess import PIPE, Popen


# Write a shell script that will run the program.
def generate_script(command, directory, script_name, system):
	# Ensure the system is recognized.
	if system.lower not in ["xf"]:
		raise NotImplementedError("Cannot run script on system other than LionXF.")
	# Ensure the script does not already exist.
	if os.path.exists(script_name):
		raise OSError("Script already exists: {}".format(script_name))
	# Write the script.
	open(script_name, "w").write("""#PBS -l pmem=8gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf
set -u
echo "Started at $(date)"
cd {d}
{c}
echo "Ended at $(date)"
rm {s}""".format(c=command, d=directory))


# Run a shell script on the queue.
def submit_script(command, directory, script_name, system):
	generate_script(command, directory, script_name, system)
	os.system("qsub {}".format(script_name))


# Run a program directly.
def run_command(command):
	process = Popen(command, stdout=PIPE)
	output, err = process.communicate()
	exit_code = process.wait()
	return exit_code
