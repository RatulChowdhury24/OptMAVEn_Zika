#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/energies/MoleculeA/KCDR3_79/parts_701_to_750.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/energies/MoleculeA/KCDR3_79/parts_701_to_750.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/energies/MoleculeA/KCDR3_79/parts_701_to_750.time.txt parts_701_to_750_interaction_energies
echo "Job ended at $(date)"
