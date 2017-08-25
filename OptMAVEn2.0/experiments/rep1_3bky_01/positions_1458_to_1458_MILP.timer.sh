#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.time.txt MILP_for_positions_1458_to_1458
rm /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.sh /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.time.txt
rm /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1458_to_1458_MILP.timer.sh
echo "Job ended at $(date)"
