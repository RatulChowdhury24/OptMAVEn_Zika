#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf

set -u
echo " "
echo " "
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/energies/MoleculeE/LV_3
vmd -dispdev none -e /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/modules/merge_antigen_part.tcl -args /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/ /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/structures/MoleculeE.pdb /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/databases/MAPs/LV/LV_3.pdb /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/energies/MoleculeE/LV_3/MoleculeE ANTI MAPP /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/input_files/top_all27_prot_na.rtf
vmd -dispdev none -e /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/modules/interaction_energies.tcl -args /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/energies/MoleculeE/LV_3/MoleculeE.psf /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/energies/MoleculeE/LV_3/MoleculeE.pdb ANTI MAPP /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/input_files/positions.dat /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/energies/MoleculeE/LV_3/energies.dat /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/Experiment_Details.txt /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/input_files/par_all27_prot_na.prm
cd /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_orig1/
python /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/programs/Optmaven2.0.py
echo " "
echo "Job ended at $(date)"
echo " "
