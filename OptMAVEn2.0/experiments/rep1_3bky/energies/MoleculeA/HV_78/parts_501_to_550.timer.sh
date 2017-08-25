#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky/
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky/energies/MoleculeA/HV_78/parts_501_to_550.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky/energies/MoleculeA/HV_78/parts_501_to_550.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky/energies/MoleculeA/HV_78/parts_501_to_550.time.txt parts_501_to_550_interaction_energies
echo "Job ended at $(date)"
