""" This script calculates the optimal combination of MAPs parts for
one antigen position using CPLEX. """

import os
import re
import sys

import STANDARDS
from OPTMAVEN import select_parts_cplex

# Get the position file from the arguments.
if len(sys.argv) != 2:
    raise OSError("Usage: args: position/file.dat")
pos_file = sys.argv[1]
if not os.path.isfile(pos_file):
    raise OSError("Position file does not exist: {}".format(pos_file))

# Get the position from the position file name.
zRot_, zRot, x_, x, y_, y, z_, z = re.match(
        "(zr)([0-9-.]+)(x)([0-9-.]+)(y)([0-9-.]+)(z)([0-9-.]+)", 
        os.path.splitext(os.path.basename(pos_file))[0]).groups()
position_str = "zr{}x{}y{}z{}".format(zRot, x, y, z)
# Make a results file for this part.
results_fields = ("zRot", "x", "y", "z", "cut", "energy", "HV", "HJ", "HCDR3",
        "KV", "KJ", "KCDR3", "LV", "LJ", "LCDR3")
results_file = os.path.join(os.path.dirname(pos_file), "parts_{}.csv".format(
        position_str))
if os.path.isfile(results_file):
    os.remove(results_file)

# Go through five rounds of optimal part selection.
solution_cuts = list()
solution_num = 5
lines = [",".join(results_fields)]
for i in range(1, solution_num + 1):
    # Solve the MILP with CPLEX.
    parts, energy = select_parts_cplex(pos_file, solution_cuts)
    # Turn the parts results into a dictionary.
    line_dict = {field: "" for field in results_fields}
    line_dict.update({name: number for name, number in zip(parts[::2],
            parts[1::2])})
    line_dict.update({"zRot": zRot, "x": x, "y": y, "z": z, "cut": i, "energy":
            energy})
    # Write the results to the results file.
    line = ",".join(map(str, [line_dict[field] for field in results_fields]))
    lines.append(line)
    # Don't find the same solution again.
    solution_cuts.append(parts)

open(results_file, "w").write("\n".join(lines))
