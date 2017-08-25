""" Generate a parameter file for OptMAVEn 2.0. """

from collections import OrderedDict


# Get a number of the specific type.
def get_number(t, what):
	if t is int:
		need = "an integer"
	elif t is float:
		need = "a float"
	else:
		raise ValueError("Type of numeric must be int or float.")
	answer = None
	while answer is None:
		try:
			answer = t(raw_input("Please enter {} value for {}: "
					.format(need, what)))
		except ValueError:
			print("{} is needed. Please try again.".format(
					need.capitalize()))
	return answer


# Get a path to a file that exists.
def get_file_path(what):
	file_path = None
	while file_path is None or not os.path.isfile(file_path):
		file_path = raw_input("Please enter the path to the {} file:\n"
				.format(what))
		if not os.path.isfile(file_path):
			print("That file does not exist. Please try again.")
	return file_path


# Get a range of numbers.
def get_range(what):
	begin = get_number(float, "first number in the range of {}"
			.format(what))
	end = None
	while end is None or end < begin:
		end = get_number(float, "last number in the range of {}"
				.format(what))
		if end < begin:
			print("The last number cannot be less than the first. "
					"Please try again.")
	if end > begin:
		step = 0.0
		while step <= 0.0:
			step = get_number(float, "step between consecutive "
					"numbers in the range of {}".format(
					what))
			if step <= 0.0:
				print("The step must be positive. Please try "
						"again.")
	else:
		step = 1.0
	return range(begin, end, step)


def OptMAVEn2(parameter_file_path):
	info = OrderedDict()
	# Get the files of the antigen and H and K chains.
	files = [("Antigen:", "antigen"), ("AntibodyH:", "heavy chain"),
			("AntibodyK:", "kappa chain")]
	for field, what in files:
		info[field] = get_file_path(what)
	# Get the translations and rotation angles.
	moves = [("ZAngles:", "z rotation angles"),
			("XTrans:", "x translations"),
			("YTrans:", "y translations"),
			("ZTrans:", "z translations")]
	for field, what in moves:
		info[field] = " ".join(["{:>5f}".format(x) for x in get_range(
				what)])
	# Write the parameter file.
	open(parameter_file_path, "w").write("\n".join(["{} {}".format(field,
			value) for field, value in info.items()]))
