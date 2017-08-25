#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/
(time -p positions_1500_to_1531_MILP.sh) 2> positions_1500_to_1531_MILP.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/bench_2g5b_batch/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py positions_1500_to_1531_MILP.time.txt MILP_for_positions_1500_to_1531
echo "Job ended at $(date)"
