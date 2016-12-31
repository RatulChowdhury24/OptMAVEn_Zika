#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf

set -u
echo " "
echo " "
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/HCDR3_405
vmd -dispdev none -e /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/modules/merge_antigen_part.tcl -args /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/ /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/structures/MoleculeE.pdb /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/databases/MAPs/HCDR3/HCDR3_405.pdb /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/HCDR3_405/MoleculeE /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/input_files/top_all27_prot_na.rtf
vmd -dispdev none -e /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/modules/interaction_energies.tcl -args /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/HCDR3_405/MoleculeE.psf /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/HCDR3_405/MoleculeE.pdb /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/input_files/positions.dat /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/HCDR3_405/energies.dat /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/Experiment_Details.txt /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/input_files/par_all27_prot_na.prm
cd /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/
python /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/programs/Optmaven2.0.py
echo " "
echo "Job ended at $(date)"
echo " "
