#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/positions_1200_to_1249_MILP.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/positions_1200_to_1249_MILP.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/positions_1200_to_1249_MILP.time.txt MILP_for_positions_1200_to_1249
echo "Job ended at $(date)"
