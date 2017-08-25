#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_2igf_20/
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_2igf_20/energies/MoleculeA/KV_50/parts_801_to_820.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_2igf_20/energies/MoleculeA/KV_50/parts_801_to_820.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_2igf_20/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_2igf_20/energies/MoleculeA/KV_50/parts_801_to_820.time.txt parts_801_to_820_interaction_energies
echo "Job ended at $(date)"
