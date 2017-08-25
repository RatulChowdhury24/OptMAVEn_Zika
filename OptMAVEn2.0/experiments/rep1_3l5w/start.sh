#PBS -l walltime=2:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/
echo "Job started on $(hostname -s) at $(date)"
(time -p python /gpfs/work/m/mfa5147/OptMAVEn2.0/programs/Optmaven.py) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/time.out
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3l5w/time.out antigen_relaxation_and_positioning
echo "Job ended at $(date)"
