#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -q lionxf

set -u
echo " "
echo " "
echo "Job started on $(hostname -s) at $(date)"
python /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/modules/select_parts.py /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/energies/MoleculeE/zr0.0x0.0y-5.0z16.0.dat
cd /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/experiments/ZIKV_mini1/
python /gpfs/scratch/mfa5147/GitHub/OptMAVEn2.0/programs/Optmaven2.0.py
echo " "
echo "Job ended at $(date)"
echo " "
