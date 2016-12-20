""" Import C++ standards into Python. """

import os
import re

# Returns tuple of (data type (string), parameter name)
parameter_def = re.compile("(string )([A-Za-z0-9])")


def get_standards():
	with open("STANDARDS2.h") as s:
		for line in s:
			if line.startswith("string"):
				name = parameter_def.match(line).groups()[1]
			elif line.strip().startswith("return"):
