#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf

set -u
echo " "
echo " "
echo "Job started on $(hostname -s) at $(date)"
cd /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/energies/MoleculeE/KCDR3_70
vmd -dispdev none -e /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/modules/merge_antigen_part.tcl -args /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/ /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/structures/MoleculeE.pdb /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/databases/MAPs/KCDR3/KCDR3_70.pdb /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/energies/MoleculeE/KCDR3_70/MoleculeE /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/input_files/top_all27_prot_na.rtf
vmd -dispdev none -e /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/modules/interaction_energies.tcl -args /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/energies/MoleculeE/KCDR3_70/MoleculeE.psf /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/energies/MoleculeE/KCDR3_70/MoleculeE.pdb /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/input_files/positions.dat /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/energies/MoleculeE/KCDR3_70/energies.dat /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/experiments/e4/Experiment_Details.txt /home/matthew/Maranas_Lab/IPRO_Suite_OptMAVEn2.0/input_files/par_all27_prot_na.prm
echo " "
echo "Job ended at $(date)"
echo " "
