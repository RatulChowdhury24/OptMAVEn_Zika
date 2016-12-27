__doc__ = """ This module is the interface between IPRO and VMD and
contains a built-in function for performing each task with VMD. """

import os
import time

import MOLECULES
import STANDARDS


# Define how to run VMD.
ModulesFolder = os.path.join(STANDARDS.InstallFolder, "modules")
InputsFolder = os.path.join(STANDARDS.InstallFolder, "input_files")
vmd_command = "vmd"
no_display_flag = "-dispdev none"
execute_flag = "-e"
molecules_flag = "-m"
frames_flag = "-f"
args_flag = "-args"

# Define how to run NAMD.
namd_command = "namd"


def split_file(file_path):
    """ Return a tuple of directory, file base, and file extension. """
    directory, file_name = os.path.split(file_path)
    file_base, extension = os.path.splitext(file_name)
    return directory, file_base, extension
    

def extracted_chains_name(molecule, chains):
    """ Generate a name for the chains extracted from a molecule. """
    # Get the name and directory of the molecule.
    directory, molecule_name, ext = split_file(molecule)
    # Add a "chains prefix."
    return os.path.join(directory, "chains{}_{}.pdb".format("".join([c.replace(
            " ", "") for c in chains]), molecule_name))


def relaxed_name(molecule):
    """ Generate a name for the molecule after it has undergone an
    energy minimization. """
    # Get the name and directory of the molecule.
    directory, molecule_name, ext = split_file(molecule)
    # Add a "relaxation prefix."
    return os.path.join(directory, "relaxed_{}.pdb".format(molecule_name))


def psf_name(molecule):
    """ Generate a name for the PSF of a molecule. """
    # Get the name and directory of the molecule.
    directory, molecule_name, ext = split_file(molecule)
    # Change the file extension.
    return os.path.join(directory, "{}.psf".format(molecule_name))


def run_vmd_script(script, molecules=None, frames=None, args=None):
    """ Run VMD using a script and optional arguments. """
    command = "{} {} {} {}".format(vmd_command, no_display_flag, execute_flag,
            os.path.join(STANDARDS.InstallFolder, "modules", script))
    if isinstance(molecules, list) or isinstance(molecules, tuple):
        molecules = " ".join(molecules)
    if isinstance(molecules, str):
        command += " {} {}".format(molecules_flag, molecules)
    if isinstance(frames, list) or isinstance(frames, tuple):
        frames = " ".join(frames)
    if isinstance(frames, str):
        command += " {} {}".format(frames_flag, frames)
    if isinstance(args, list) or isinstance(args, tuple):
        args = " ".join(args)
    if isinstance(args, str):
        command += " {} {}".format(args_flag, args)
    i = os.system(command) #FIXME: use subprocess.Popen instead
    print command
    if i != 0:
        raise Exception("Running VMD with this command has failed:\n{}".format(
                command))


def run_namd(configuration_file):
    """ Run NAMD using a configuration file. """
    # FIXME: use subprocess.Popen instead
    command = "{} {}".format(namd_command, configuration_file)
    i = os.system(command)
    if i != 0:
        raise Exception("Running NAMD with this command has failed:\n{}".format(
                command))


def extract_chains(input_coords, chains, output_coords):
    """ Extract specific chains from a coordinate file. """
    # The molecule's existence should have been verified, but make sure.
    if not os.path.isfile(input_coords):
        raise IOError("File does not exist: {}".format(input_coords))
    if isinstance(chains, (list, tuple)):
        chains = " ".join(chains)
    run_vmd_script(os.path.join(ModulesFolder, "extract_chains.tcl"), molecules=
            input_coords, args="{} {}".format(output_coords, chains))

'''
def extract_experiment_chains(experiment):
    """ Extract the selected chain of each molecule in an experiment
    and put the output files in the experiment's folder. """
    output_files = list()
    for molecule in experiment["Molecules"]:
        # Clean the molecule to remove any residues that have been excluded.
        # The molecule is at position 2 in the list.
        molecule_file = os.path.join(experiment["Folder"], molecule[2].
                generate_name(procedure="", fileFormat="PDB"))
        molecule[2].output(name=molecule_file)
        # The selected chain is item 1.
        chain = molecule[1]
        output_file = extracted_chains_name(molecule_file, chain)
        extract_chains(molecule_file, chain, output_file)
        output_files.append(output_file)
    return output_files
'''

def generate_PSF(experiment, molecule):
    """ Add missing atoms to the antigen coordinates and generate a PSF
    of the complete structure. """
    # Determine the directory in which to find and output the files.
    if os.path.isdir(os.path.join(experiment["Folder"], "structures")):
        struct_directory = os.path.join(experiment["Folder"], "structures")
    else:
        struct_directory = experiment["Folder"]
    if os.path.isdir(os.path.join(experiment["Folder"], "input_files")):
        inputs_directory = os.path.join(experiment["Folder"], "input_files")
    else:
        inputs_directory = experiment["Folder"]
    # Generate the names of the files.
    input_coords = os.path.join(struct_directory, molecule.generate_name(
            fileFormat="PDB", procedure=""))
    output_coords = os.path.join(struct_directory, molecule.generate_name(
            fileFormat="PDB", procedure=""))
    output_struct = os.path.join(struct_directory, psf_name(molecule.generate_name(
            fileFormat="PDB")))
    topology_files = [os.path.join(inputs_directory, f) for f in experiment[
            "CHARMM Topology Files"]]
    # Format the argument list.
    args = [input_coords, output_coords, output_struct] + topology_files
    # Run VMD with the make_antigen_psf script.
    run_vmd_script(os.path.join(ModulesFolder, "make_antigen_psf.tcl"), args=args)
    # Make sure the output files were generated.
    if not all(os.path.isfile(output) for output in [output_coords,
            output_struct]):
        raise Exception("VMD could not generate the structure and/or coordinate"
                " files.")
    # Return a molecule that contains the coordinates of added atoms.
    return MOLECULES.Molecule(open(output_coords).readlines())


def generate_PSFs(experiment):
    """ Generate a PSF of each molecule in the experiment. """
    for molecule in experiment["Molecules"]:
        # The molecule object is at position 2 of each item in the list of
        # molecules.
        generate_PSF(experiment, molecule[2])


def relax_molecule(experiment, molecule):
    """ Use VMD to perform an energy minimization on a molecule in an
    experiment. """
    # The directories of the experiment.
    struct_directory = os.path.join(experiment["Folder"], "structures")
    inputs_directory = os.path.join(experiment["Folder"], "input_files")
    # Generate the name of the molecule.
    molecule_coords = os.path.join(struct_directory, molecule.generate_name())
    molecule_struct = os.path.join(struct_directory, os.path.splitext(
            molecule.generate_name())[0] + ".psf")
    molecule_prefix = os.path.splitext(molecule_coords)[0]
    molecule_out = molecule_prefix + ".coor"
    # Get the names of the parameter file.
    parameter_files = [os.path.join(inputs_directory, f) for f in experiment[
            "CHARMM Parameter Files"]]
    parameters = {
	    "structure": molecule_struct,
	    "coordinates": molecule_coords,
	    "parameters": parameter_files,
	    "set outputname": molecule_prefix
    }
    # Get the name of the NAMD configuration file.
    base_conf = os.path.join(STANDARDS.InstallFolder, "input_files",
            "namd_relaxation_base.conf")
    namd_conf = os.path.join(inputs_directory, "namd_relaxation.conf")
    # Copy the base configuration file and add the details of this experiment.
    lines = open(base_conf).readlines()
    format_line = lambda param, val: "{:<20s}{}\n".format(parameter, value)
    with open(namd_conf, "w") as f:
	    for line in lines:
		    if line.strip() != "":
			    # If the line begins with a parameter in the dict of parameters,
			    # add the value of that parameter to the configuration file.
			    parameter = line[0: min(20, len(line))].strip()
			    if parameter in parameters:
			        if isinstance(parameters[parameter], list):
			            for i, value in enumerate(parameters[parameter]):
			                if i < len(parameters[parameter]) - 1:
			                    f.write(format_line(parameter, value))
			                else:
			                    line = format_line(parameter, value)
			        else:
			            line = "{:<20s}{}\n".format(parameter, parameters[
			                    parameter])
		    f.write(line)
    # Run the minimization in NAMD.
    run_namd(namd_conf)
    # Rename the output coordinates to clearly be a PDB.
    os.rename(molecule_out, molecule_coords)
    # Remove other NAMD output files.
    for ext in ("vel", "xsc", "xst"):
        try:
            os.remove("{}.{}".format(molecule_prefix, ext))
        except OSError:
            pass


def relax_molecules(experiment):
    """ Use VMD to perform an energy minimization on each molecule in
    an experiment. """
    for molecule in experiment["Molecules"]:
        relax_molecule(experiment, molecule[2])
        

def initial_antigen_position(experiment, molecule):
    """ Move an antigen to its initial position. """
    struct_directory = os.path.join(experiment["Folder"], "structures")
    molecule_file = os.path.join(struct_directory, molecule.generate_name())
    name = molecule.name
    details_file = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    run_vmd_script("initial_antigen_position.tcl", molecules=molecule_file,
            args=[molecule_file, details_file])
    

def initial_antigen_positions(experiment):
    """ Move all of the antigens in an experiment to their initial
    positions. """
    for molecule in experiment["Molecules"]:
        initial_antigen_position(experiment, molecule[2])
