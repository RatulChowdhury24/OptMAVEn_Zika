__doc__ = """ Check the status of all of your OptMAVEn experiments. """

from collections import defaultdict, OrderedDict
import os
import sys

InstallFolder = "/gpfs/work/m/mfa5147/OptMAVEn2.0"

if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = None

status_counts = defaultdict(int)
experiment_status = OrderedDict()
exps_dir = os.path.join(InstallFolder, "experiments")
if os.path.isdir(exps_dir):
    for experiment in sorted(os.listdir(exps_dir)):
        exp_dir = os.path.join(exps_dir, experiment)
        details_file = os.path.join(exp_dir, "Experiment_Details.txt")
        if (name is None or name in experiment) and os.path.isfile(details_file):
            status_file = os.path.join(exp_dir, "status")
            if os.path.isfile(status_file):
                with open(status_file) as f:
                    status = f.read().strip()
            else:
                status = "initializing"
            status_counts[status] += 1
            experiment_status[experiment] = status
print
print "STATUS OF EXPERIMENTS:"
print "\n".join(["{}: {}".format(exp, status) for exp, status in experiment_status.iteritems()])
print
print "SUMMARY:"
print "\n".join(["{}: {}".format(status, count) for status, count in status_counts.iteritems()])
print "TOTAL EXPERIMENTS: {}".format(sum(status_counts.values()))
