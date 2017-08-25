#!/usr/bin/env python

InstallFolder = "/gpfs/work/m/mfa5147/OptMAVEn2.0"

# The name of this program
__name__ = "The IPRO Suite Experiment Starter"

# Documentation string
__doc__ = """
Written in 2013 by Robert Pantazes of the Costas Maranas Lab in the Chemical
Engineering Department of the Pennsylvania State University

This program creates an Experiment class object, which asks the user for all
information needed to run an experiment. Once that has been gathered, it goes to
the Experiment's folder, puts the proper program there, and optionally starts
the Experiment."""

# Include PYTHON information
import os
import sys
sys.path.append(os.path.join(InstallFolder, "modules"))
# Include the EXPERIMENT module
import EXPERIMENT

# Make the experiment, using user input
experiment = EXPERIMENT.Experiment(False)
# Move to the experiment's folder and output it's information
os.chdir(experiment["Folder"])
experiment.output()
# Start running OptMAVEn from the experiment's folder.
starter_script = os.path.join(experiment["Folder"], "start.sh")
time_file = os.path.join(experiment["Folder"], "time.out")
hours = 2
open(starter_script, "w").write("""#PBS -l walltime={h}:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u
cd {f}
echo "Job started on $(hostname -s) at $(date)"
(time -p python {i}/programs/Optmaven.py) 2> {t}
python {i}/modules/PERFORMANCE.py {t} antigen_relaxation_and_positioning
echo "Job ended at $(date)"
""".format(h=hours, f=experiment["Folder"], i=InstallFolder, t=time_file))
os.system("chmod u+x {}".format(starter_script))
# Start the experiment, and time the startup.
os.system("qsub {}".format(starter_script))
