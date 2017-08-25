import os
import sys
import tempfile

import numpy as np

start_dir = os.path.dirname(os.path.realpath(__file__))
main_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
mods_dir = os.path.join(main_dir, "modules")
exps_dir = os.path.join(main_dir, "experiments")

if len(sys.argv) > 1:
    exp_spec = sys.argv[1]
else:
    exp_spec = ""

sys.path.append(mods_dir)
import PERFORMANCE
import CHARMM
import EXPERIMENT
import MOLECULES

output = os.path.join(start_dir, "results.csv")
open(output, "w").write("Antigen,Atoms,Residues,Epitope,Positions,Disk Usage (kb),CPU Time (s),Real Time (s),Initialization Time (s),Relax Time (s),Clash Time (s),Mount Time (s),Energy Time (s),MILP Time (s),Clustering Time (s),Pre-Relaxed NAMD Energy (kcal/mol),Relaxed NAMD Energy (kcal/mol),Pre-Relaxed CHARMM Energy (kcal/mol),Relaxed CHARMM Energy (kcal/mol)")

for exp in [e for e in os.listdir(exps_dir) if exp_spec in e]:
    exp_dir = os.path.join(exps_dir, exp)
    details = os.path.join(exp_dir, "Experiment_Details.txt")
    positions = os.path.join(exp_dir, "input_files/positions.dat")
    performance = os.path.join(exp_dir, "performance.txt")
    results = os.path.join(exp_dir, "results")
    complexes = os.path.join(results, "MoleculeA_assembly_info.txt")
    if all(map(os.path.exists, [performance, positions, details, results, complexes])):
        print exp
        os.chdir(exp_dir)
        experiment = EXPERIMENT.Experiment()
        directory = experiment["Folder"]
        # Get disk usage
        handle, du_file = tempfile.mkstemp()
        os.system("du -sk {} > {}".format(directory, du_file))
        du = int(open(du_file).read().split()[0])
        try:
            os.remove(du_file)
        except OSError:
            pass
        p_file = PERFORMANCE.File(performance)
        with open(os.path.join(directory, "structures", "mounted_MoleculeA.pdb")) as f:
            residue_numbers = [line[22: 26] for line in f if line.startswith("ATOM")]
            n_atoms = len(residue_numbers)
            n_residues = len(set(residue_numbers))
        n_epitope = len([line for line in open(details) if line.startswith(
                "Epitope Position")])
        n_pos = len(open(positions).readlines())
        cpu_time = p_file.clock_time()
        try:
            real_time = p_file.duration()
        except TypeError:
            real_time = np.nan
        init_time = p_file.time_process_category("antigen_relaxation_and_positioning", "time")
        relax_time = p_file.time_process_category("NAMD relaxation of molecule", "sub_time", indexes=[0])
        clash_time = p_file.get_sub_time("culling_clashes")
        mount_time = p_file.get_sub_time("Mounting antigen A") 
        energy_time = p_file.time_process_category(
                "interaction_energies", "time")
        milp_time = p_file.time_process_category("MILP", "time")
        cluster_time = p_file.time_process_category("clustering", "sub_time")
        charmm = [np.nan] * 2
        namd = dict()
        try:
            min_assembly = [int(datum.split(":")[1]) for datum in open(complexes).readline().split(", ") if datum.startswith("assembly")][0]
        except IndexError:
            pass
        else:
            """
            pre_relaxed_file = os.path.join(results, "assembly_{}_MoleculeA.pdb".format(min_assembly))
            relaxed_file = os.path.join(results, "relaxed_assembly_{}_MoleculeA.pdb".format(min_assembly))
            for i, f in enumerate([pre_relaxed_file, relaxed_file]):
                remove = list()
                pdb = "molecule.pdb"
                os.system("cp {} {}".format(f, pdb))
                remove.append(pdb)
                for inf in experiment["CHARMM Topology Files"] + experiment["CHARMM Parameter Files"]:
                    os.system("cp {} {}".format(os.path.join(main_dir, "input_files", inf), exp_dir))
                    remove.append(os.path.basename(inf))
                mf = MOLECULES.MoleculeFile(pdb)
                ag = mf[0]
                ab1 = mf[1]
                ab2 = mf[2]
                # Solvation is only used on the relaxed structure
                experiment["Use Solvation"] = (f == relaxed_file)
                complex_energy = CHARMM.Energy([ag, ab1, ab2], experiment)
                ag_energy = CHARMM.Energy(ag, experiment)
                ab_energy = CHARMM.Energy([ab1, ab2], experiment)
                print complex_energy, ag_energy, ab_energy
                charmm[i] = complex_energy - ag_energy - ab_energy
                for inf in remove:
                    os.remove(inf)
            """
            namd = {line.split(":")[0]: float(line.split(":")[1]) for line in open(os.path.join(results, "assembly_{}_MoleculeA_energy.dat".format(min_assembly)))}
        text = "\n" + ",".join(map(str, [exp, n_atoms, n_residues, n_epitope, n_pos, du, cpu_time, real_time, init_time, relax_time, clash_time, mount_time, energy_time, milp_time, cluster_time, namd.get("Pre-Relaxed", np.nan), namd.get("Relaxed", np.nan), charmm[0], charmm[1]]))
        open(output, "a").write(text)

os.chdir(start_dir) 
