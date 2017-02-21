""" This program finds the series of import calls to custom modules help resolve import hell. """

import os
import sys
from collections import OrderedDict

# The program at which to start.
main = sys.argv[1]
# Other directories in which to search for modules.
module_dirs = sys.argv[2:]

# Make a list of modules.
modules = list()
for d in [os.getcwd()] + module_dirs:
    modules.extend([os.path.join(d, module) for module in os.listdir(d) if module.endswith(".py")])

module_names = {os.path.basename(module)[: -3]: module for module in modules}

print module_names
def read_program(program):
    # Read the contents of a program and follow import statements.
    path = module_names.get(program)
    raw_input(program)
    if path is None:
       return
    imports = list()
    with open(path) as f:
        for line in f:
            if "import" in line:
                if line.strip().startswith("import"):
                    imports.extend([m.strip().split()[0] for m in line[6:].split(",")])
                elif line.strip().startswith("from"):
                    imports.append(line.strip().split()[1])
    return OrderedDict([(program, read_program(program)) for program in imports])

read_program(main[: -3])
