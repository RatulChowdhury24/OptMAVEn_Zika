import os
import shutil
import time


#mfa#
# The install folder should have subdirectories called "databases",
# "input_files", "modules", "programs", and "share".
subdirs = ("databases", "input_files", "modules", "programs", "share")

#mfa#
def validInstallFolder(path):
	""" Check to see if a directory is a valid install folder. """
	return all([os.path.isdir(os.path.join(path, sd)) for sd in subdirs])


#mfa#
# The syntax that defines the install folder varies depending on the file type.
install_syntax = {".py": "InstallFolder =", ".tcl": "set InstallFolder"}

def install_line(ext, line, install_folder):
	""" Format a line to define the install folder. """
	# Lines that define the install folder contain this string.
	try:
		syntax = install_syntax[ext]
	except KeyError:
		raise IOError('Unknown file extension: "{}"'.format(ext))
	split_line = line.split(syntax)
	# If the line defines the install folder, then define the install
	# folder on that line.	
	if len(split_line) > 1:
		return '{}{} "{}"\n'.format(split_line[0], syntax, install_folder)
	# Otherwise, return the line unaltered.
	else:
		return line


def install_file(file_, install_folder):
	""" Install a file by initializing its InstallFolder. """
	# Determine the type of file.
	path, ext = os.path.splitext(file_)
	# Format all of the lines.
	lines = "".join([install_line(ext, line, install_folder) for line in
			open(file_).readlines()])
	# Write all of the formatted lines back into the file.
	open(file_, "w").write(lines)
	

# Remind the user to copy the appropriate files.
os.system("clear")
message = """REMINDER: Before running the experiment, please copy all force
field inputs to "input_files" and all initial structures to "structures"."""
print message

#mfa#
# Get the actual install folder. Test the current directory first.
InstallFolder = os.getcwd()
parent, current = os.path.split(InstallFolder)
# If the current directory is not valid, test every directory above it.
while not validInstallFolder(InstallFolder) and parent != current:
	InstallFolder = parent
	parent, current = os.path.split(InstallFolder)
# Exit if no valid install folder was found.
if not validInstallFolder(InstallFolder):
	message = "Cannot locate a valid install folder, which should have "\
			"subdirectories called {}. Please enter the directory "\
			"of the IPRO Suite and rerun install.py.".format(
			", ".join(subdirs))
	print message
	exit()

# Edit all of the programs.
ProgramsFolder = os.path.join(InstallFolder, "programs")
edits = [os.path.join(ProgramsFolder, prgm) for prgm in os.listdir(
		ProgramsFolder)]
# Also edit the STANDARDS.py module and the TCL modules.
ModulesFolder = os.path.join(InstallFolder, "modules")
modules = [module for module in os.listdir(ModulesFolder) if module.endswith(
        ".tcl") or module in ("STANRADRS.py",)]
for module in modules:
	edits.append(os.path.join(ModulesFolder, module))
# Change the InstallFolder to the actual folder for each file.
for file_ in edits:
    install_file(file_, InstallFolder)

# Copy Start_Experiment.py to the InstallFolder.
shutil.copy(os.path.join(ProgramsFolder, "Start_Experiment.py"), InstallFolder)

# Allow the user to read the previous reminder.
raw_input("Press Enter to continue.")
os.system("clear")  # Clear the message.

# Confirm that the user has successfully finished installation.
message = "The IPRO Suite has been installed successfully into {}. Press Enter"\
		" to continue.".format(InstallFolder)
raw_input(message)
