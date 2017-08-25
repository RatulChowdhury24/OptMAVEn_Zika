__doc__ = """ This module is the interface between IPRO and Klaus Schulten
programs VMD and NAMD and contains a built-in function for performing each task
 with VMD and NAMD. """

from collections import OrderedDict 
import itertools
import os
import subprocess
import tempfile
import time

import numpy as np

import MOLECULES
import CHARMM
import STANDARDS
import SUBMITTER
import PERFORMANCE
import MAPs

# Define how to run VMD.
ModulesFolder = os.path.join(STANDARDS.InstallFolder, "modules")
InputsFolder = os.path.join(STANDARDS.InstallFolder, "input_files")
vmd_command = "vmd"
no_display_flag = "-dispdev none"
execute_flag = "-e"
molecules_flag = "-m"
frames_flag = "-f"
args_flag = "-args"
AgSeg = "ANTI"
MAPsSeg = "MAPP"
AbHSeg = "ABHC"
AbLSeg = "ABLC"

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


def make_vmd_command(script, molecules=None, frames=None, args=None,
        directory=None, time_=None):
    """ Create a VMD command. """
    command = "{} {} {} {}".format(vmd_command, no_display_flag, execute_flag,
            os.path.join(STANDARDS.InstallFolder, "modules", script))
    if isinstance(molecules, list) or isinstance(molecules, tuple):
        molecules = " ".join(map(str, molecules))
    if isinstance(molecules, str):
        command += " {} {}".format(molecules_flag, molecules)
    if isinstance(frames, list) or isinstance(frames, tuple):
        frames = " ".join(map(str, frames))
    if isinstance(frames, str):
        command += " {} {}".format(frames_flag, frames)
    if isinstance(args, list) or isinstance(args, tuple):
        args = " ".join(map(str, args))
    if isinstance(args, str):
        command += " {} {}".format(args_flag, args)
    if isinstance(time_, tuple):
        time_file, process = time_
        command = "(time {c}) 2> {t}\npython {i}/modules/PERFORMANCE.py {t} {p}"\
                .format(c=command, t=time_file, i=STANDARDS.InstallFolder,
                p=process)
    if isinstance(directory, str):
        command = "cd {}\n{}".format(directory, command)
    return command


def run_vmd_script(script, molecules=None, frames=None, args=None, time_=None):
    """ Run VMD using a script and optional arguments. """
    command = make_vmd_command(script, molecules, frames, args, time_=time_)
    i = os.system(command) #FIXME: use subprocess.Popen instead
    if i != 0:
        raise Exception("Running VMD with this command has failed:\n{}".format(
                command))

'''
def queue_vmd_script(script, molecules=None, frames=None, args=None, directory=
        None, time_=None):
    """ Run VMD using a script and optional arguments. """
    command = make_vmd_command(script, molecules, frames, args, time_=time_)
    SUBMITTER.experiment_script(command, time_=time_)


def queue_vmd_scripts(argument_list):
    """ Run a series of VMD scripts given by the argument list. """
    command = "\n".join([make_vmd_command(**entry) for entry in argument_list])
    SUBMITTER.experiment_script(command, time_=time_)
'''
 
    
def run_namd(configuration_file, output_prefix=None, run=True):
    """ Run NAMD using a configuration file. """
    # FIXME: use subprocess.Popen instead
    command = "{} {}".format(namd_command, configuration_file)
    if output_prefix is not None:
        command += " > {} 2> {}".format(output_prefix + ".out", output_prefix +
                ".err")
    if run:
        i = os.system(command)
        if i != 0:
            raise Exception("Running NAMD with this command has failed:\n{}".format(
                    command))
    else:
        return command


def read_namd_energies(energy_file, key="TOTAL", expected=None):
    """ Read the latest energy from a NAMD energy file. """
    index = None
    energies = list()
    with open(energy_file + ".out") as f:
        for line in f:
            # Get the index of the energy term.
            if line.startswith("ETITLE:"):
                index = line.split().index(key)
            elif line.startswith("ENERGY:"):
                energies.append(float(line.split()[index]))
    if expected is not None and len(energies) != expected:
        raise Exception("NAMD generated {} energies (expected {})".format(
                len(energies), expected))
    return energies


'''
def extract_chains(input_coords, chains, output_coords):
    """ Extract specific chains from a coordinate file. """
    # The molecule's existence should have been verified, but make sure.
    if not os.path.isfile(input_coords):
        raise IOError("File does not exist: {}".format(input_coords))
    if isinstance(chains, (list, tuple)):
        chains = " ".join(chains)
    run_vmd_script(os.path.join(ModulesFolder, "extract_chains.tcl"), molecules=
            input_coords, args="{} {}".format(output_coords, chains))


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
    output_struct = os.path.join(struct_directory, psf_name(molecule.
            generate_name(fileFormat="PDB")))
    topology_files = [os.path.join(inputs_directory, f) for f in experiment[
            "CHARMM Topology Files"]]
    # Format the argument list.
    args = [input_coords, output_coords, output_struct, AgSeg] + topology_files
    # Run VMD with the make_antigen_psf script.
    run_vmd_script(os.path.join(ModulesFolder, "make_antigen_psf.tcl"),
            args=args)
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


def relax_molecule(experiment, molecule_coords, molecule_struct):
    """ Use NAMD to perform an energy minimization on a molecule in an
    experiment. """
    timer = PERFORMANCE.Timer()
    timer.start()
    performance_file = PERFORMANCE.performance_file(experiment)
    # Generate the name of the molecule.
    head, tail = os.path.split(os.path.splitext(molecule_coords)[0])
    molecule_prefix = os.path.join(head, "relaxed_" + tail)
    # Copy the original coordinates of the molecule to a new file on which the
    # relaxation will be performed.
    molecule_out = molecule_prefix + ".coor"
    os.system("cp {} {}".format(molecule_coords, molecule_out))
    # Define some parameters for controlling the length of the minimization.
    tolerance = 0.0001  # Once the energy changes by less than the tolerance,
    # the minimization will stop.
    max_iterations = 2000  # In case the tolerance is not met, the maximum
    # number of iterations the minimization can run.
    steps_per_cycle = 25
    cycles_per_check = 2
    check_each = steps_per_cycle * cycles_per_check  # The number of iterations
    # run between energy checks.
    struct_directory = os.path.dirname(molecule_struct)
    namd_out = os.path.join(struct_directory, "namd")  # NAMD output
    # Get the names of the parameter file.
    inputs_directory = os.path.join(STANDARDS.InstallFolder, "input_files")
    parameter_files = [os.path.join(inputs_directory, f) for f in experiment[
            "CHARMM Parameter Files"]]
    parameters = {
	    "structure": molecule_struct,
	    "coordinates": molecule_out,
	    "parameters": parameter_files,
	    "set outputname": molecule_prefix,
            "outputEnergies": 1,
            "outputPressure": 1,
            "stepspercycle": steps_per_cycle,
            "fullElectFrequency": steps_per_cycle,
            "minimize": check_each
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
    running_energy = list()
    running_time = list()
    iteration = 0
    total_time = 0.0
    while (len(running_energy) < 2 or running_energy[-1] - running_energy[-2]
                < -abs(tolerance)) and iteration < max_iterations:
        start_time = time.time()
        run_namd(namd_conf, namd_out)
        iteration += check_each
        total_time += time.time() - start_time
        running_time.append(total_time)
        running_energy.append(read_namd_energies(namd_out, "TOTAL",
                check_each + 1)[-1])
        status = "Iteration {:>5} Time {:>3f} Energy {:>5f}".format(iteration,
                total_time, running_energy[-1])
        # Use this to benchmark against CHARMM.
        benchmark = False
        if benchmark:
            for infile in experiment["CHARMM Topology Files"] + experiment["CHARMM Parameter Files"]:
                os.system("ln {} {}".format(os.path.join(experiment["folder"], "input_files", infile), infile))
            CHARMM_energy = CHARMM.Energy(MOLECULES.MoleculeFile(os.path.basename(molecule_out))["E"])
            status += " CHARMM Energy {:>5f}".format(CHARMM_energy)
        print status
        open(os.path.join(experiment["Folder"], "structures", molecule_prefix +
            "_relaxation.dat"), "a").write(status + "\n")
    #run_namd(namd_conf, namd_out)  ## Use this w/ subprocess.Popen when that works.
    performance_file.sub_time("NAMD relaxation of molecule {}".format(molecule_coords), timer)
    # Rename the output coordinates to clearly be a PDB.
    os.rename(molecule_out, molecule_prefix + ".pdb")
    # Remove other NAMD output files.
    for ext in ("vel", "xsc", "xst"):
        try:
            os.remove("{}.{}".format(molecule_prefix, ext))
        except OSError:
            pass


def relax_molecules(experiment):
    """ Use VMD to perform an energy minimization on each molecule in
    an experiment. """
    struct_directory = os.path.join(experiment["Folder"], "structures")
    inputs_directory = os.path.join(experiment["Folder"], "input_files")
    for molecule_list in experiment["Molecules"]:
        molecule = molecule_list[2]
        molecule_coords = os.path.join(struct_directory, molecule.generate_name())
        molecule_struct = os.path.join(struct_directory, os.path.splitext(
                molecule.generate_name())[0] + ".psf")
        relax_molecule(experiment, molecule_coords, molecule_struct)
        

def initial_antigen_position(experiment, molecule):
    """ Move an antigen to its initial position. """
    timer = PERFORMANCE.Timer()
    timer.start()
    struct_directory = os.path.join(experiment["Folder"], "structures")
    input_molecule = os.path.join(struct_directory, "relaxed_" +
            molecule.generate_name())
    output_molecule = os.path.join(struct_directory, "mounted_" +
            molecule.generate_name()) 
    name = molecule.name
    details_file = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    final_centers = os.path.join(experiment["Folder"], "final_centers.txt")
    run_vmd_script("initial_antigen_position.tcl", molecules=input_molecule,
            args=[output_molecule, AgSeg, details_file, final_centers])
    # Record the final epitope z coordinate.
    with open(final_centers) as f:
        f.readline()
        antigen_center = np.array(map(float, f.readline().split()))
        f.readline()
        epitope_center = np.array(map(float, f.readline().split()))
    # The epitope should be at the origin.
    if not np.all(np.isclose(epitope_center, 0, atol=0.001)):
        raise Exception("VMD has failed to move the epitope to the origin.")
    # If the z coordinates have been minimized, the antigen will be above the
    # epitope.
    if not np.all(np.isclose(antigen_center[0: 2], 0, atol=0.001)):
        raise Exception("VMD has failed to minimize the z coordinates of the epitope.")
    try:
    	os.remove(final_centers)
    except OSError:
        pass
    performance_file = PERFORMANCE.performance_file(experiment)
    performance_file.sub_time("Mounting antigen {}".format(name), timer)


def initial_antigen_positions(experiment):
    """ Move all of the antigens in an experiment to their initial
    positions. """
    for molecule in experiment["Molecules"]:
        initial_antigen_position(experiment, molecule[2])


def cull_clash(experiment, molecule):
    """ Find the antigen positions that do not cause clashes between
    the antigen and the antibodies. """
    timer = PERFORMANCE.Timer()
    timer.start()
    # Get the location of the antigen.
    Ag = os.path.join(experiment["Folder"], "structures",
            "mounted_" + molecule.generate_name())
    # Define the locations of the prototype heavy and kappa chain antibodies.
    chains = ("H", "K")
    Ig = {chain: os.path.join(STANDARDS.InstallFolder, "input_files",
            "Molecule{}.pdb".format(chain)) for chain in chains}
    # Make sure that all of the files exist.
    missing = [chain for chain, path in Ig.items() if not os.path.exists(path)]
    if len(missing) > 0:
        raise IOError("Cannot locate files for prototype antibody chains: {}."
                .format(", ".join(missing)))
    # Define the output file.
    posFile = os.path.join(experiment["Folder"], "input_files", "positions.dat")
    # Define the experiment details file.
    detFile = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    detFile = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    # Define the clash cutoff (Angstroms).
    clashCutoff = 1.25
    # Define the maximum number of clashes permitted.
    clashesPermitted = 0
    # Run VMD.
    run_vmd_script("cull_clashes.tcl", molecules=[Ag, Ig["H"], Ig["K"]], args=[
            detFile, posFile, clashCutoff, clashesPermitted])
    # Ensure that the position file exists.
    if not os.path.isfile(posFile):
        raise Exception("VMD failed to generate a file of the non-clashing "
                "positions.")
    # Ensure that there is at least one non-clashing position.
    num_positions = len(open(posFile).readlines())
    if num_positions == 0:
        raise Exception("There were no positions that did not clash. Please "
                "re-examine the structure of the antigen and adjust your grid "
                "settings.")
    performance_file = PERFORMANCE.performance_file(experiment)
    performance_file.sub_time("culling_clashes", timer)


def cull_clashes(experiment):
    """ Find the antigen positions that do not cause clashes between
    the antigen and the antibodies. """
    for molecule in experiment["Molecules"]:
        cull_clash(experiment, molecule[2])


def prepare_antigen_part(experiment, antigen, part_file, prefix, time_=None):
    """ Combine the structures and coordinates of the antigen and a
    part in the MAPs database. Write, but do not run, the command. """
    # Define the files.
    details = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    antigen = os.path.join(experiment["Folder"], "structures", "mounted_" +
            antigen.generate_name())
    args = [experiment["Folder"], antigen, part_file, prefix, AgSeg, MAPsSeg] +\
            [os.path.join(STANDARDS.InstallFolder, "input_files", f) for f in
            experiment["CHARMM Topology Files"]]
    return make_vmd_command("merge_antigen_part.tcl", args=args, time_=time_)


def MAPs_interaction_energy(structure, coordinates, positions, output, details,
    parameters):
    """ Write, but do not run a command to calculate the interaction
    energy between the antigen and a MAPs part. """
    if isinstance(parameters, (list, tuple)):
        parameters = " ".join(parameters)
    return make_vmd_command("interaction_energies.tcl", args=[structure,
            coordinates, AgSeg, MAPsSeg, positions, output, details, parameters])


def MAPs_interaction_energies(experiment, batch_size=1):
    """ Calculate the interaction energy between a MAPs part and the
    antigen in each antigen position. """
    """ Create a structure and coordinate file of the antigen combined
    with each MAPs part. """
    # Make a directory to store the interaction energies.
    ie_directory = os.path.join(experiment["Folder"], "energies")
    if not os.path.isdir(ie_directory):
        os.mkdir(ie_directory)
    # Get some information from the experiment.
    details = os.path.join(experiment["Folder"], "Experiment_Details.txt")
    parameter_files = [os.path.join(STANDARDS.InstallFolder, "input_files", f)
            for f in experiment["CHARMM Parameter Files"]]
    # Loop through each antigen.
    for antigen in [mol[2] for mol in experiment["Molecules"]]:
        name = os.path.splitext(antigen.generate_name())[0]
        mol_directory = os.path.join(ie_directory, name)
        # Make a sub-directory for each antigen.
        if not os.path.exists(mol_directory):
            os.mkdir(mol_directory)
        # Make a sub-directory for each MAPs part.
        MAPs_parts = MAPs.list_parts()
        n_parts = len(MAPs_parts)
        for first in range(1, n_parts + 1, batch_size):
            commands = ["date"]
            last = min(first + batch_size - 1, n_parts)
            part_range = "parts_{}_to_{}".format(first, last)
            for part in MAPs_parts[first - 1: last]:
                part_directory = os.path.join(mol_directory, part)
                if not os.path.isdir(part_directory):
                    os.mkdir(part_directory)
                commands.append("cd {}".format(part_directory))
                # Combine the antigen and MAPs part.
                part_file = MAPs.get_part_file(part)
                output_prefix = os.path.join(part_directory, name)
                commands.append(prepare_antigen_part(experiment, antigen,
			            part_file, output_prefix))
                # Calculate the interaction energies.
                struct_file = output_prefix + ".psf"
                coords_file = output_prefix + ".pdb"
                positions_file = os.path.join(experiment["Folder"],
                        "input_files", "positions.dat")
                energy_file = os.path.join(part_directory, "energies.dat")
                commands.append(MAPs_interaction_energy(struct_file,
                        coords_file, positions_file, energy_file, details,
                        parameter_files))
            # Call back OptMAVEn once finished.
            commands.append("cd {}".format(experiment["Folder"]))
            commands.append("date")
            commands.append("python {}".format(os.path.join(
                    STANDARDS.InstallFolder, "programs", "Optmaven.py")))
            # Run the calculations on the queue.
            command = "\n".join(commands)
            script = os.path.join(part_directory, "{}.sh".format(part_range))
            SUBMITTER.experiment_script(experiment, script, command,
                    directory=experiment["Folder"], time_=
                    "parts_{}_to_{}_interaction_energies".format(first, last))


def assemble_antibody_antigen_complex(experiment, antigen, assembly, d, calculate_relaxed_energy=True):
    folder = experiment["Folder"]
    topology_files = [os.path.join(STANDARDS.InstallFolder, "input_files", f)
            for f in experiment["CHARMM Topology Files"]]
    args = [folder, os.path.join(folder, "Experiment_Details.txt"), d["zRot"],
            d["x"], d["y"], d["z"], antigen]
    args.extend(map(MAPs.get_part_file, [d["HV"], d["HJ"], d["HCDR3"], d["LV"],
            d["LJ"], d["LCDR3"]]))
    args.extend([assembly, AgSeg, AbHSeg, AbLSeg] + topology_files)
    run_vmd_script(os.path.join(ModulesFolder,
            "assemble_antibody_antigen_complex.tcl"), args=args)
    coords = assembly + ".pdb"
    structure = assembly + ".psf"
    energies = {"Pre-Relaxed": antibody_complex_energy(experiment, coords, structure)}
    if calculate_relaxed_energy:
        relaxed_coords = os.path.join(os.path.dirname(coords), "relaxed_" + os.path.basename(coords))
        if not os.path.isfile(relaxed_coords):
            relax_molecule(experiment, coords, structure)
        energies["Relaxed"] = antibody_complex_energy(experiment, relaxed_coords, structure)
    open(os.path.join(folder, "results", os.path.splitext(os.path.basename(assembly))[0] + "_energy.dat"), "w").write("\n".join(["{}: {}".format(key, value) for key, value in energies.items()]))
    return energies


def assemble_antibody_antigen_complexes(experiment, antibodies, calculate_relaxed_energies=True):
    folder = experiment["Folder"]
    for antigen in [mol[2] for mol in experiment["Molecules"]]:
        basename = os.path.splitext(os.path.basename(antigen.generate_name()))[0]
        antigen_file = os.path.join(folder, "structures", "mounted_" + basename + ".pdb")
        count_file = os.path.join(folder, "results", basename + "_assembly_info.txt")
        # Copy the antibody information to a new list of new dicts tbecause the dicts will be modified.
        antibodies_ = [OrderedDict(ab.items()) for ab in antibodies]
        index = "assembly"
        pre_relaxed = "Pre-Relaxed"
        relaxed = "Relaxed"
        for i, d in enumerate(antibodies_):
            d[index] = i + 1
            assembly_file = os.path.join(folder, "results", "assembly_{}_{}".format(d[index], basename))
            energies = assemble_antibody_antigen_complex(experiment, antigen_file, assembly_file, d, calculate_relaxed_energy=(calculate_relaxed_energies is True))
            d[pre_relaxed] = energies[pre_relaxed]
            if calculate_relaxed_energies is True:
                d[relaxed] = energies[relaxed]
            else:
                d[relaxed] = np.nan
        # Sort the antibodies by their pre-relaxation interaction energies.
        antibodies_.sort(key=lambda x: x[pre_relaxed])
        # Relax the number of specified complexes.
        if calculate_relaxed_energies is not True:
            if calculate_relaxed_energies is False:
                calculate_relaxed_energies = 0
            for i in range(min(calculate_relaxed_energies, len(antibodies_))):
                d = antibodies_[i]
                assembly_file = os.path.join(folder, "results", "assembly_{}_{}".format(d[index], basename))
                energies = assemble_antibody_antigen_complex(experiment, antigen_file, assembly_file, d, calculate_relaxed_energy=True)
                d[relaxed] = energies[relaxed]
        # Write the results to a file.            
        open(count_file, "w").write("\n".join([", ".join(["{}: {}".format(field, value) for field, value in d.items()]) for d in antibodies_]))
            

def antibody_complex_energy(experiment, coords, struct):
    energy = 0.0
    results = os.path.join(experiment["Folder"], "results")
    for AbSeg in (AbHSeg, AbLSeg):
        handle, energy_file = tempfile.mkstemp(prefix="energy_", suffix=".txt")
        args = [coords, struct, AgSeg, AbSeg, energy_file] + [os.path.join(STANDARDS.InstallFolder, "input_files", f) for f in experiment["CHARMM Parameter Files"]]
        run_vmd_script("interaction_energy.tcl", args=args)
        try:
            energy += float(open(energy_file).read())
            os.remove(energy_file)
        except OSError:
            raise Exception("Calculating the antigen-antibody interaction energy has failed.")
    return energy


def MAPs_part_clashes(cut_file, clash_cutoff=1.0):
    """ This function is only used for the pre-processing step of
    determining which pairs of MAPs parts clash sterically and
    developing integer cuts. """
    parts = MAPs.get_parts()
    # Choose all combinations of categories.
    # Do not check for clashes between parts of the same category.
    finished = True
    text = ""
    for cat1, cat2 in itertools.combinations(parts, 2):
        # Also, do not check for clashes between one kappa and one lambda part.
        if "".join(sorted(cat1[0] + cat2[0])) != "KL":
            result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
            if os.path.isfile(result):
                with open(result) as f:
                    text += f.read().strip() + "\n"
            else:
                script = "{}_{}_{}.sh".format(cut_file, cat1, cat2)
                command = "python {}/databases/MAPs/Make_Integer_Cuts.py {} {} {} {}".format(STANDARDS.InstallFolder, cat1, cat2, result, clash_cutoff)
                SUBMITTER.experiment_script(None, script, command, self_destruct=True)
                finished = False
    if finished:
        with open(cut_file + ".txt", "w") as f:
            f.write(text)
        for cat1, cat2 in itertools.combinations(parts, 2):
            if "".join(sorted(cat1[0] + cat2[0])) != "KL":
                result = "{}_{}_{}.txt".format(cut_file, cat1, cat2)
                try:
                    os.remove(result)
                except OSError:
                    pass
    

def MAPs_part_category_clashes(cat1, cat2, cut_file, clash_cutoff=1.0):
    parts = MAPs.get_parts()
    f = open(cut_file, "w")
    f.close()
    for part1, part2 in itertools.product(parts[cat1], parts[cat2]):
        pdbs = map(MAPs.get_part_file, [part1, part2])
        handle, clash_file = tempfile.mkstemp()
        run_vmd_script("count_clashes.tcl", molecules=pdbs,
                args=[clash_file, clash_cutoff])
        with open(clash_file) as f2:
            clashes = int(f2.read())
        if clashes > 0:
            with open(cut_file, "a") as f:
                f.write("{} {}\n".format(part1, part2).replace("_", " "))
        os.remove(clash_file)
